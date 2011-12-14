// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Volume averaged quantities required by the Richards model.
 */
#ifndef DUMUX_RICHARDS_VOLUME_VARIABLES_HH
#define DUMUX_RICHARDS_VOLUME_VARIABLES_HH

#include "richardsproperties.hh"

#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/boxmodels/common/boxvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \ingroup BoxVolumeVariables
 * \brief Volume averaged quantities required by the Richards model.
 *
 * This contains the quantities which are are constant within a finite
 * volume in the Richards model
 */
template <class TypeTag>
class RichardsVolumeVariables : public BoxVolumeVariables<TypeTag>
{
    typedef BoxVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, RichardsIndices) Indices;
    enum {
        pwIdx = Indices::pwIdx,
        numPhases = FluidSystem::numPhases,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    //! The type returned by the fluidState() method
    typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> FluidState;

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The primary variables as a vector for the finite
     *                volume.
     * \param problem The physical problem which needs to be solved.
     * \param element The DUNE Codim<0> enitity which intersects
     *                the control volume of the box method
     * \param elemGeom The element's finite volume geometry
     * \param scvIdx The local index of the sub control volume inside the element
     * \param isOldSol Specifies whether the solution is from
     *                 the previous time step or from the current one
     */
    void update(const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        assert(!FluidSystem::isLiquid(nPhaseIdx));

        ParentType::update(elemCtx, scvIdx, timeIdx);

        completeFluidState(fluidState_, elemCtx, scvIdx, timeIdx);

        //////////
        // specify the other parameters
        //////////
        const auto &problem = elemCtx.problem();
        const MaterialLawParams &matParams =
            problem.materialLawParams(elemCtx, scvIdx, timeIdx);
        MaterialLaw::relativePermeabilities(relativePermeability_, matParams, fluidState_);
        porosity_ = problem.porosity(elemCtx, scvIdx, timeIdx);

        // energy related quantities not belonging to the fluid state
        asImp_().updateEnergy_(elemCtx, scvIdx, timeIdx);
    }

    /*!
     * \copydoc BoxModel::completeFluidState
     */
    static void completeFluidState(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int timeIdx)
    {
        Implementation::updateTemperature_(fluidState, elemCtx, scvIdx, timeIdx);

        // material law parameters
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        const auto &problem = elemCtx.problem();
        const typename MaterialLaw::Params &materialParams =
            problem.materialLawParams(elemCtx, scvIdx, timeIdx);
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);

        /////////
        // calculate the pressures
        /////////
                    
        // first, we have to find the minimum capillary pressure (i.e. Sw = 0)
        fluidState.setSaturation(wPhaseIdx, 1.0);
        fluidState.setSaturation(nPhaseIdx, 0.0);
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState);
                    
        // non-wetting pressure can be larger than the
        // reference pressure if the medium is fully
        // saturated by the wetting phase
        Scalar pW = priVars[pwIdx];
        Scalar pN = std::max(elemCtx.problem().referencePressure(elemCtx, scvIdx, /*timeIdx=*/0),
                             pW + (pC[nPhaseIdx] - pC[wPhaseIdx]));
                    
        /////////
        // calculate the saturations
        /////////
        fluidState.setPressure(wPhaseIdx, pW);
        fluidState.setPressure(nPhaseIdx, pN);

        PhaseVector sat;
        MaterialLaw::saturations(sat, materialParams, fluidState);
        fluidState.setSaturation(wPhaseIdx, sat[wPhaseIdx]);
        fluidState.setSaturation(nPhaseIdx, 1.0 - sat[wPhaseIdx]);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        // compute and set the wetting phase viscosity
        Scalar mu = FluidSystem::viscosity(fluidState, paramCache, wPhaseIdx);
        fluidState.setViscosity(wPhaseIdx, mu);
        fluidState.setViscosity(nPhaseIdx, 1e-20);

        // compute and set the wetting phase density
        Scalar rho = FluidSystem::density(fluidState, paramCache, wPhaseIdx);
        fluidState.setDensity(wPhaseIdx, rho);
        fluidState.setDensity(nPhaseIdx, 1e-20);

        Implementation::updateEnthalpy_(fluidState,
                                        paramCache,
                                        elemCtx,
                                        scvIdx,
                                        timeIdx);
    }

    /*!
     * \brief Returns a reference to the fluid state for the volume
     */
    const FluidState &fluidState() const
    { return fluidState_; }


    /*!
     * \brief Returns the average porosity [] within the control volume.
     *
     * The porosity is defined as the ratio of the pore space to the
     * total volume, i.e. \f[ \Phi := \frac{V_{pore}}{V_{pore} + V_{rock}} \f]
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns relative permeability [-] of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return relativePermeability(phaseIdx)/fluidState().viscosity(phaseIdx); }

    /*!
     * \brief Returns the effective capillary pressure \f$\mathrm{[Pa]}\f$ within the
     *        control volume.
     *
     * The capillary pressure is defined as the difference in
     * pressures of the non-wetting and the wetting phase, i.e.
     * \f[ p_c = p_n - p_w \f]
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(nPhaseIdx) - fluidState_.pressure(wPhaseIdx); }

    /*!
     * \brief Given a fluid state, set the temperature in the primary variables
     */
    template <class FluidState>
    static void setPriVarTemperatures(PrimaryVariables &priVars, const FluidState &fs)
    { }

    /*!
     * \brief Set the enthalpy rate per second of a rate vector, .
     */
    static void setEnthalpyRate(RateVector &rateVec, Scalar rate)
    {}
    
    /*!
     * \brief Given a fluid state, set the enthalpy rate which emerges
     *        from a volumetric rate.
     */
    template <class FluidState>
    static void setEnthalpyRate(RateVector &rateVec,
                                const FluidState &fluidState, 
                                int phaseIdx, 
                                Scalar volume)
    { };

    static void updateTemperature_(FluidState &fluidState,
                                   const ElementContext &elemCtx,
                                   int scvIdx,
                                   int timeIdx)
    {
        fluidState.setTemperature(elemCtx.problem().temperature(elemCtx, scvIdx, timeIdx));
    };

    template<class ParameterCache>
    static void updateEnthalpy_(FluidState &fluidState,
                                const ParameterCache &paramCache,
                                const ElementContext &elemCtx,
                                int scvIdx,
                                int timeIdx)
    { }

protected:
    /*!
     * \brief Called by update() to compute the energy related quantities
     */
    void updateEnergy_(const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx)
    { }

    FluidState fluidState_;
    Scalar relativePermeability_[numPhases];
    Scalar porosity_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
