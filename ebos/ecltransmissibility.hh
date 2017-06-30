// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::EclTransmissibility
 */
#ifndef EWOMS_ECL_TRANSMISSIBILITY_HH
#define EWOMS_ECL_TRANSMISSIBILITY_HH

#include <ewoms/common/propertysystem.hh>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/GridProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/TransMult.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>


#include <opm/grid/CpGrid.hpp>

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/ConditionalStorage.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

BEGIN_PROPERTIES

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Vanguard);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ElementMapper);
NEW_PROP_TAG(EnableEnergy);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This class calculates the transmissibilites for grid faces according to the
 *        Eclipse Technical Description.
 */
template <class TypeTag>
class EclTransmissibility
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GridView::Intersection Intersection;

    static const bool enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy);

    // Grid and world dimension
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    static const unsigned elemIdxShift = 32; // bits

public:

    EclTransmissibility(const Vanguard& vanguard)
        : vanguard_(vanguard)
    {
        const Opm::UnitSystem& unitSystem = vanguard_.deck().getActiveUnitSystem();
        transmissibilityThreshold_  = unitSystem.parse("Transmissibility").getSIScaling() * 1e-6;
    }

    /*!
     * \brief Actually compute the transmissibilty over a face as a pre-compute step.
     *
     * This code actually uses the direction specific "centroids" of
     * each element. These "centroids" are _not_ the identical
     * barycenter of the element, but the middle of the centers of the
     * faces of the logical Cartesian cells, i.e., the centers of the
     * faces of the reference elements. We do things this way because
     * the barycenter of the element can be located outside of the
     * element for sufficiently "ugly" (i.e., thin and "non-flat")
     * elements which in turn leads to quite wrong
     * permeabilities. This approach is probably not always correct
     * either but at least it seems to be much better.
     */
    void finishInit()
    { update(); }


    /*!
     * \brief Compute all transmissibilities
     *
     * Also, this updates the "thermal half transmissibilities" if energy is enabled.
     */
    void update()
    {
        const auto& gridView = vanguard_.gridView();

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        const ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
#else
        const ElementMapper elemMapper(gridView);
#endif

        const unsigned numElements = elemMapper.size();

        const std::vector<double> ntg = this->prepareTranCalc_(numElements);

        if (! this->gridIsRadial_())
            this->computeTrans_NewTran_(gridView, elemMapper, numElements, ntg);
        else
            this->computeTrans_OldTran_Radial_(gridView, elemMapper, numElements, ntg);
    }

    /*!
     * \brief Return the permeability for an element.
     */
    const DimMatrix& permeability(unsigned elemIdx) const
    { return permeability_[elemIdx]; }

    /*!
     * \brief Return the transmissibility for the intersection between two elements.
     */
    Scalar transmissibility(unsigned elemIdx1, unsigned elemIdx2) const
    { return trans_.at(isId_(elemIdx1, elemIdx2)); }

    /*!
     * \brief Return the transmissibility for a given boundary segment.
     */
    Scalar transmissibilityBoundary(unsigned elemIdx, unsigned boundaryFaceIdx) const
    { return transBoundary_.at(std::make_pair(elemIdx, boundaryFaceIdx)); }

    /*!
     * \brief Return the thermal "half transmissibility" for the intersection between two
     *        elements.
     *
     * The "half transmissibility" features all sub-expressions of the "thermal
     * transmissibility" which can be precomputed, i.e. they are not dependent on the
     * current solution:
     *
     * H_t = A * (n*d)/(d*d);
     *
     * where A is the area of the intersection between the inside and outside elements, n
     * is the outer unit normal, and d is the distance between the center of the inside
     * cell and the center of the intersection.
     */
    Scalar thermalHalfTrans(unsigned insideElemIdx, unsigned outsideElemIdx) const
    { return thermalHalfTrans_->at(directionalIsId_(insideElemIdx, outsideElemIdx)); }

    Scalar thermalHalfTransBoundary(unsigned insideElemIdx, unsigned boundaryFaceIdx) const
    { return thermalHalfTransBoundary_.at(std::make_pair(insideElemIdx, boundaryFaceIdx)); }

private:
    struct EntityID
    {
        template <class Entity, class CartMapper>
        EntityID(const Entity&        entity,
                 const ElementMapper& elemMapper,
                 const CartMapper&    cartMapper)
            : elemIdx(elemMapper.index(entity))
            , cartElemIdx(cartMapper.cartesianIndex(elemIdx))
        {}

        unsigned int elemIdx;
        unsigned int cartElemIdx;
    };

    /*!
     * Compute radial transmissibilities assuming fully matching RADIAL
     * grids.
     */
    class OldTran_Radial
    {
    public:
        explicit OldTran_Radial(EclTransmissibility& host)
            : host_(host)
        {}

        void computeTrans(const GridView& gridView,
                          const ElementMapper& elemMapper,
                          const unsigned int numElements,
                          const std::vector<double>& ntg)
        {
            const auto hTrans =
                this->computeHalfTrans_(gridView, elemMapper, numElements);

            this->computeFullTrans_(gridView, elemMapper, hTrans, ntg);
        }

    private:
        template <typename T>
        struct RemoveCVR
        {
            using type = typename std::remove_cv<
                typename std::remove_reference<T>::type
                >::type;
        };

        template <typename T>
        using RemoveCVR_t = typename RemoveCVR<T>::type;

        using ECLGridType = RemoveCVR_t<
            decltype(std::declval<Vanguard>().eclState().getInputGrid())>;

        using CartMapper = RemoveCVR_t<
            decltype(std::declval<Vanguard>().cartesianIndexMapper())>;

        using TransMultiplier = RemoveCVR_t<
            decltype(std::declval<Vanguard>().eclState().getTransMult())>;

        EclTransmissibility& host_;

        std::vector<double> angleChangeVector_() const
        {
            auto dtheta = std::vector<double>{};

            const auto& eclGrid = this->eclGrid_();

            auto pointAngleDeg =
                [&eclGrid](const std::size_t j, const std::size_t c) -> double
            {
                const double r2d =
                    180.0 / 3.14159265358979323846264338327950288;

                const auto& p = eclGrid.getCornerPos(0, j, 0, c);

                return r2d * std::atan2(p[1], p[0]);
            };

            const auto ny = eclGrid.getNY();

            dtheta.reserve(ny);

            // Angle increases in the clockwise direction.
            double t2 = pointAngleDeg(0, 0);
            for (auto j = 0*ny; j < ny; ++j) {
                const double t1 = t2;
                t2 = pointAngleDeg(j, 1);

                if (t2 > t1)
                    // We crossed (-\infty, 0) real line; unwrap.
                    t2 -= 360.0;

                dtheta.push_back(- (t2 - t1));
            }

            return dtheta;
        }

        std::vector<double> radiusVector_() const
        {
            auto radii = std::vector<double>{};

            const auto& eclGrid = this->eclGrid_();

            auto pointRadius =
                [&eclGrid](const std::size_t i, const std::size_t c) -> double
            {
                const auto& p = eclGrid.getCornerPos(i, 0, 0, c);
                return std::hypot(p[0], p[1]);
            };

            // This code relies on the fact that local corner 0 is on the
            // inner plane (minimum radius) while local corner 1 is on the
            // outer plane (maximum radius).

            const auto nx = eclGrid.getNX();

            radii.reserve(nx + 1);
            radii.push_back(pointRadius(0, 0));

            for (auto i = 0*nx; i < nx; ++i)
                radii.push_back(pointRadius(i, 1));

            return radii;
        }

        std::vector<double> thicknessVector_() const
        {
            auto dz = std::vector<double>{};

            const auto& eclGrid = this->eclGrid_();

            auto pointDepth =
                [&eclGrid](const std::size_t k, const std::size_t c) -> double
            {
                const auto& p = eclGrid.getCornerPos(0, 0, k, c);
                return p[2];
            };

            const auto nz = eclGrid.getNY();

            dz.reserve(nz);

            double z2 = pointDepth(0, 0);
            for (auto k = 0*nz; k < nz; ++k) {
                const double z1 = z2;

                z2 = pointDepth(k, 4);

                dz.push_back(z2 - z1);
            }

            return dz;
        }

        const ECLGridType& eclGrid_() const
        {
            return this->host_.vanguard_.eclState().getInputGrid();
        }

        const CartMapper& cartMapper_() const
        {
            return this->host_.vanguard_.cartesianIndexMapper();
        }

        const TransMultiplier& transMult_() const
        {
            return this->host_.vanguard_.eclState().getTransMult();
        }

        std::vector<Scalar>
        computeHalfTrans_(const GridView& gridView,
                          const ElementMapper& elemMapper,
                          const unsigned int numElements) const
        {
            // Fully matching radial grid => each element (cell) has exactly
            // six faces/connections/intersections.
            const auto facesPerElm = 6 + 0*numElements;
            auto hTrans = std::vector<Scalar>(numElements * facesPerElm);

            const auto& cartMapper = this->cartMapper_();

            const auto R      = this->radiusVector_();
            const auto dTheta = this->angleChangeVector_();
            const auto dZ     = this->thicknessVector_();

            for (auto elemIt    = gridView.template begin</*codim=*/ 0>(),
                      elemEndIt = gridView.template end</*codim=*/ 0>();
                 elemIt != elemEndIt; ++elemIt)
            {
                const auto E = EntityID(*elemIt, elemMapper, cartMapper);

                auto ijk = std::array<int,3>{};
                cartMapper.cartesianCoordinate(E.elemIdx, ijk);

                const auto& K = this->host_.permeability_[E.elemIdx];

                const auto r1 = R[ijk[0] + 0];  // Inner
                const auto r2 = R[ijk[0] + 1];  // Outer
                const auto dt = dTheta[ijk[1]]; // Degrees
                const auto dz = dZ[ijk[2]];

                const auto log_r2r1 = std::log(r2 / r1);
                const auto r12 = r1 * r1;
                const auto r22 = r2 * r2;

                auto* hT = &hTrans[E.elemIdx*facesPerElm + 0];

                // Radial direction.
                {
                    const auto half = 1.0 / 2.0;

                    const auto tR_num = K[0][0] * dt * dz;

                    const auto D_inner =
                        (r22 / (r22 - r12))*log_r2r1 - half;

                    const auto D_outer =
                        half - (r12 / (r22 - r12))*log_r2r1;

                    // R- (centre towards inner radius).
                    hT[0] = tR_num / D_inner;

                    // R+ (centre towards outer radius).
                    hT[1] = tR_num / D_outer;
                }

                // Azimuthal direction (symmetric).
                hT[2] = hT[3] = 2 * K[1][1] * dz * log_r2r1 / dt;

                // Vertical direction (symmetric).
                hT[4] = hT[5] = K[2][2] * dt * (r22 - r12) / dz;
            }

            return hTrans;
        }

        void computeFullTrans_(const GridView& gridView,
                               const ElementMapper& elemMapper,
                               const std::vector<Scalar>& hTrans,
                               const std::vector<double>& ntg)
        {
            const auto facesPerElem = 0*hTrans.size() + 6;

            const auto& cartMapper = this->cartMapper_();
            const auto& transMult = this->transMult_();

            for (auto elemIt    = gridView.template begin</*codim=*/ 0>(),
                      elemEndIt = gridView.template end</*codim=*/ 0>();
                 elemIt != elemEndIt; ++elemIt)
            {
                const auto& elem = *elemIt;
                const auto inside = EntityID(elem, elemMapper, cartMapper);

                const auto* hTi = &hTrans[inside.elemIdx*facesPerElem + 0];

                for (auto isIt    = gridView.ibegin(elem),
                          isEndIt = gridView.iend(elem);
                     isIt != isEndIt; ++isIt)
                {
                    const auto& intersection = *isIt;

                    if (! intersection.neighbor())
                        continue;

                    const auto outside = EntityID(intersection.outside(),
                                                  elemMapper, cartMapper);

                    if (inside.elemIdx > outside.elemIdx)
                        continue;

                    const auto* hTo = &hTrans[outside.elemIdx*facesPerElem + 0];
                    const unsigned insideFaceIdx = intersection.indexInInside();
                    const unsigned outsideFaceIdx = intersection.indexInOutside();

                    auto halfTrans1 = hTi[insideFaceIdx];
                    auto halfTrans2 = hTo[outsideFaceIdx];

                    this->host_.computeFullTrans_(std::move(halfTrans1),
                                                  std::move(halfTrans2),
                                                  inside, insideFaceIdx,
                                                  outside, outsideFaceIdx,
                                                  transMult, ntg);
                }
            }
        }
    };

    std::vector<double> prepareTranCalc_(const unsigned numElements)
    {
        // get the ntg values, the ntg values are modified for the cells merged with minpv
        std::vector<double> ntg;
        minPvFillNtg_(ntg);

        extractPermeability_();

        // reserving some space in the hashmap upfront saves quite a bit of time because
        // resizes are costly for hashmaps and there would be quite a few of them if we
        // would not have a rough idea of how large the final map will be (the rough idea
        // is a conforming Cartesian grid).
        trans_.clear();
        trans_.reserve(numElements*3*1.05);

        transBoundary_.clear();

        // if energy is enabled, let's do the same for the "thermal half transmissibilities"
        if (enableEnergy) {
            thermalHalfTrans_->clear();
            thermalHalfTrans_->reserve(numElements*6*1.05);

            thermalHalfTransBoundary_.clear();
        }

        return ntg;
    }

    void computeTrans_NewTran_(const GridView& gridView,
                               const ElementMapper& elemMapper,
                               const unsigned int numElements,
                               const std::vector<double>& ntg)
    {
        const auto& eclState = vanguard_.eclState();
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        auto& transMult = eclState.getTransMult();

        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();

        // calculate the axis specific centroids of all elements
        std::array<std::vector<DimVector>, dimWorld> axisCentroids;
        {
            const auto& eclGrid = eclState.getInputGrid();

            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                axisCentroids[dimIdx].resize(numElements);

            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;
                unsigned elemIdx = elemMapper.index(elem);

                // compute the axis specific "centroids" used for the transmissibilities. for
                // consistency with the flow simulator, we use the element centers as
                // computed by opm-parser's Opm::EclipseGrid class for all axes.
                unsigned cartesianCellIdx = cartMapper.cartesianIndex(elemIdx);
                const auto& centroid = eclGrid.getCellCenter(cartesianCellIdx);
                for (unsigned axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
                    for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                        axisCentroids[axisIdx][elemIdx][dimIdx] = centroid[dimIdx];
            }
        }

        // compute the transmissibilities for all intersections
        elemIt = gridView.template begin</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            unsigned elemIdx = elemMapper.index(elem);

            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
            unsigned boundaryIsIdx = 0;
            for (; isIt != isEndIt; ++ isIt) {
                // store intersection, this might be costly
                const auto& intersection = *isIt;

                // deal with grid boundaries
                if (intersection.boundary()) {
                    // compute the transmissibilty for the boundary intersection
                    const auto& geometry = intersection.geometry();
                    const auto& faceCenterInside = geometry.center();

                    auto faceAreaNormal = intersection.centerUnitOuterNormal();
                    faceAreaNormal *= geometry.volume();

                    Scalar transBoundaryIs;
                    computeHalfTrans_(transBoundaryIs,
                                      faceAreaNormal,
                                      intersection.indexInInside(),
                                      distanceVector_(faceCenterInside,
                                                      intersection.indexInInside(),
                                                      elemIdx,
                                                      axisCentroids),
                                      permeability_[elemIdx]);

                    // normally there would be two half-transmissibilities that would be
                    // averaged. on the grid boundary there only is the half
                    // transmissibility of the interior element.
                    transBoundary_[std::make_pair(elemIdx, boundaryIsIdx)] = transBoundaryIs;

                    // for boundary intersections we also need to compute the thermal
                    // half transmissibilities
                    if (enableEnergy) {
                        const auto& n = intersection.centerUnitOuterNormal();
                        const auto& inPos = elem.geometry().center();
                        const auto& outPos = intersection.geometry().center();
                        const auto& d = outPos - inPos;

                        // eWoms expects fluxes to be area specific, i.e. we must *not*
                        // the transmissibility with the face area here
                        Scalar thermalHalfTrans = std::abs(n*d)/(d*d);

                        thermalHalfTransBoundary_[std::make_pair(elemIdx, boundaryIsIdx)] =
                            thermalHalfTrans;
                    }

                    ++ boundaryIsIdx;
                    continue;
                }

                if (!intersection.neighbor())
                    // elements can be on process boundaries, i.e. they are not on the
                    // domain boundary yet they don't have neighbors.
                    continue;

                const auto inside = EntityID(intersection.inside(), elemMapper, cartMapper);
                const auto outside = EntityID(intersection.outside(), elemMapper, cartMapper);

                // update the "thermal half transmissibility" for the intersection
                if (enableEnergy) {
                    const auto& n = intersection.centerUnitOuterNormal();
                    Scalar A = intersection.geometry().volume();

                    const auto& inPos = elem.geometry().center();
                    const auto& outPos = intersection.geometry().center();
                    const auto& d = outPos - inPos;

                    (*thermalHalfTrans_)[directionalIsId_(inside.elemIdx, outside.elemIdx)] =
                        A * (n*d)/(d*d);
                }

                // we only need to calculate a face's transmissibility
                // once...
                if (inside.elemIdx > outside.elemIdx)
                    continue;

                // local indices of the faces of the inside and
                // outside elements which contain the intersection
                int insideFaceIdx  = intersection.indexInInside();
                int outsideFaceIdx = intersection.indexInOutside();

                if (insideFaceIdx == -1) {
                    // NNC. Set zero transmissibility, as it will be
                    // *added to* by applyNncToGridTrans_() later.
                    assert(outsideFaceIdx == -1);
                    trans_[isId_(elemIdx, outsideElemIdx)] = 0.0;
                    continue;
                }

                DimVector faceCenterInside;
                DimVector faceCenterOutside;
                DimVector faceAreaNormal;

                typename std::is_same<Grid, Dune::CpGrid>::type isCpGrid;
                computeFaceProperties(intersection,
                                      inside.elemIdx,
                                      insideFaceIdx,
                                      outside.elemIdx,
                                      outsideFaceIdx,
                                      faceCenterInside,
                                      faceCenterOutside,
                                      faceAreaNormal,
                                      isCpGrid);

                Scalar halfTrans1;
                Scalar halfTrans2;

                computeHalfTrans_(halfTrans1,
                                  faceAreaNormal,
                                  insideFaceIdx,
                                  distanceVector_(faceCenterInside,
                                                  insideFaceIdx,
                                                  inside.elemIdx,
                                                  axisCentroids),
                                  permeability_[inside.elemIdx]);
                computeHalfTrans_(halfTrans2,
                                  faceAreaNormal,
                                  outsideFaceIdx,
                                  distanceVector_(faceCenterOutside,
                                                  outsideFaceIdx,
                                                  outside.elemIdx,
                                                  axisCentroids),
                                  permeability_[outside.elemIdx]);

                // convert half transmissibilities to full face
                // transmissibilities using the harmonic mean.  Assign
                // result to appropriate location in 'trans_'.
                computeFullTrans_(std::move(halfTrans1),
                                  std::move(halfTrans2),
                                  inside, insideFaceIdx,
                                  outside, outsideFaceIdx,
                                  transMult, ntg);
            }
        }

        // potentially overwrite and/or modify  transmissibilities based on input from deck
        updateFromEclState_();

        // Create mapping from global to local index
        const size_t cartesianSize = cartMapper.cartesianSize();
        // reserve memory
        std::vector<int> globalToLocal(cartesianSize, -1);

        // loop over all elements (global grid) and store Cartesian index
        elemIt = vanguard_.grid().leafGridView().template begin<0>();

        for (; elemIt != elemEndIt; ++elemIt) {
            int elemIdx = elemMapper.index(*elemIt);
            int cartElemIdx = vanguard_.cartesianIndexMapper().cartesianIndex(elemIdx);
            globalToLocal[cartElemIdx] = elemIdx;
        }
        applyEditNncToGridTrans_(elemMapper, globalToLocal);
        applyNncToGridTrans_(globalToLocal);

        //remove very small non-neighbouring transmissibilities
        removeSmallNonCartesianTransmissibilities_();
    }

    void computeTrans_OldTran_Radial_(const GridView& gridView,
                                      const ElementMapper& elemMapper,
                                      const unsigned int numElements,
                                      const std::vector<double>& ntg)
    {
        this->computeTrans_OldTran_Radial_(typename std::is_same<Grid, Dune::CpGrid>::type{},
                                           gridView, elemMapper, numElements, ntg);
    }

    void removeSmallNonCartesianTransmissibilities_()
    {
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        for (auto&& trans: trans_) {
            if (trans.second < transmissibilityThreshold_) {
                const auto& id = trans.first;
                const auto& elements = isIdReverse_(id);
                int gc1 = std::min(cartMapper.cartesianIndex(elements.first), cartMapper.cartesianIndex(elements.second));
                int gc2 = std::max(cartMapper.cartesianIndex(elements.first), cartMapper.cartesianIndex(elements.second));

                // only adjust the NNCs
                if (gc2 - gc1 == 1 || gc2 - gc1 == cartDims[0] || gc2 - gc1 == cartDims[0]*cartDims[1])
                    continue;

                //remove transmissibilities less than the threshold (by default 1e-6 in the deck's unit system)
                trans.second = 0.0;
            }
        }
    }

    void applyAllZMultipliers_(Scalar& trans,
                               unsigned insideFaceIdx,
                               unsigned insideCartElemIdx,
                               unsigned outsideCartElemIdx,
                               const Opm::TransMult& transMult,
                               const std::array<int, dimWorld>& cartDims)
    {
        if (insideFaceIdx > 3) { // top or or bottom
            Scalar mult = 1e20;
            unsigned cartElemIdx = insideCartElemIdx;
            // pick the smallest multiplier while looking down the pillar untill reaching the other end of the connection
            // for the inbetween cells we apply it from both sides
            while (cartElemIdx != outsideCartElemIdx) {
                if (insideFaceIdx == 4 || cartElemIdx !=insideCartElemIdx )
                    mult = std::min(mult, transMult.getMultiplier(cartElemIdx, Opm::FaceDir::ZMinus));
                if (insideFaceIdx == 5 || cartElemIdx !=insideCartElemIdx)
                    mult = std::min(mult, transMult.getMultiplier(cartElemIdx, Opm::FaceDir::ZPlus));

                cartElemIdx += cartDims[0]*cartDims[1];
            }
            trans *= mult;
        }
        else
            applyMultipliers_(trans, insideFaceIdx, insideCartElemIdx, transMult);
    }

    void updateFromEclState_()
    {
        const auto& gridView = vanguard_.gridView();
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        const auto& properties = vanguard_.eclState().get3DProperties();

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
#else
        ElementMapper elemMapper(gridView);
#endif

        const auto& inputTranx = properties.getDoubleGridProperty("TRANX");
        const auto& inputTrany = properties.getDoubleGridProperty("TRANY");
        const auto& inputTranz = properties.getDoubleGridProperty("TRANZ");

        // compute the transmissibilities for all intersections
        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();

        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                // store intersection, this might be costly
                const auto& intersection = *isIt;
                if (!intersection.neighbor())
                    continue; // intersection is on the domain boundary

                unsigned c1 = elemMapper.index(intersection.inside());
                unsigned c2 = elemMapper.index(intersection.outside());

                if (c1 > c2)
                    continue; // we only need to handle each connection once, thank you.

                auto isId = isId_(c1, c2);

                int gc1 = std::min(cartMapper.cartesianIndex(c1), cartMapper.cartesianIndex(c2));
                int gc2 = std::max(cartMapper.cartesianIndex(c1), cartMapper.cartesianIndex(c2));

                if (gc2 - gc1 == 1) {
                    if (inputTranx.deckAssigned())
                        // set simulator internal transmissibilities to values from inputTranx
                        trans_[isId] = inputTranx.iget(gc1);
                    else
                        // Scale transmissibilities with scale factor from inputTranx
                        trans_[isId] *= inputTranx.iget(gc1);
                }
                else if (gc2 - gc1 == cartDims[0]) {
                    if (inputTrany.deckAssigned())
                        // set simulator internal transmissibilities to values from inputTrany
                        trans_[isId] = inputTrany.iget(gc1);
                    else
                        // Scale transmissibilities with scale factor from inputTrany
                        trans_[isId] *= inputTrany.iget(gc1);
                }
                else if (gc2 - gc1 == cartDims[0]*cartDims[1]) {
                    if (inputTranz.deckAssigned())
                        // set simulator internal transmissibilities to values from inputTranz
                        trans_[isId] = inputTranz.iget(gc1);
                    else
                        // Scale transmissibilities with scale factor from inputTranz
                        trans_[isId] *= inputTranz.iget(gc1);
                }
                //else.. We don't support modification of NNC at the moment.
            }
        }
    }


    template <class Intersection>
    void computeFaceProperties(const Intersection& intersection,
                               const int insideElemIdx,
                               const int insideFaceIdx,
                               const int outsideElemIdx,
                               const int outsideFaceIdx,
                               DimVector& faceCenterInside,
                               DimVector& faceCenterOutside,
                               DimVector& faceAreaNormal,
                               /*isCpGrid=*/std::false_type) const
    {
        // default implementation for DUNE grids
        const auto& geometry = intersection.geometry();
        faceCenterInside = geometry.center();
        faceCenterOutside = faceCenterInside;

        faceAreaNormal = intersection.centerUnitOuterNormal();
        faceAreaNormal *= geometry.volume();
    }

    template <class Intersection>
    void computeFaceProperties(const Intersection& intersection,
                               const int insideElemIdx,
                               const int insideFaceIdx,
                               const int outsideElemIdx,
                               const int outsideFaceIdx,
                               DimVector& faceCenterInside,
                               DimVector& faceCenterOutside,
                               DimVector& faceAreaNormal,
                               /*isCpGrid=*/std::true_type) const
    {
        int faceIdx = intersection.id();
        faceCenterInside = vanguard_.grid().faceCenterEcl(insideElemIdx, insideFaceIdx);
        faceCenterOutside = vanguard_.grid().faceCenterEcl(outsideElemIdx, outsideFaceIdx);
        faceAreaNormal = vanguard_.grid().faceAreaNormalEcl(faceIdx);
    }

    /*
     * \brief Applies additional transmissibilities specified via NNC keyword.
     *
     * Applies only those NNC that are actually represented by the grid. These may
     * NNCs due to faults or NNCs that are actually neighbours. In both cases that
     * specified transmissibilities (scaled by EDITNNC) will be added to the already
     * existing models.
     *
     * \param cartesianToCompressed Vector containing the compressed index (or -1 for inactive
     *                              cells) as the element at the cartesian index.
     * \return Two vector of NNCs (scaled by EDITNNC). The first one are the NNCs that have been applied
     *         and the second the NNCs not resembled by faces of the grid. NNCs specified for
     *         inactive cells are omitted in these vectors.
     */
    std::tuple<std::vector<Opm::NNCdata>, std::vector<Opm::NNCdata> >
    applyNncToGridTrans_(const std::vector<int>& cartesianToCompressed)
    {
        // First scale NNCs with EDITNNC.
        std::vector<Opm::NNCdata> unprocessedNnc;
        std::vector<Opm::NNCdata> processedNnc;
        const auto& nnc = vanguard_.eclState().getInputNNC();
        if (!nnc.hasNNC())
            return make_tuple(processedNnc, unprocessedNnc);

        auto nncData = nnc.nncdata();
        auto editnncData = vanguard_.eclState().getInputEDITNNC().data();
        auto nncLess =
            [](const Opm::NNCdata& d1, const Opm::NNCdata& d2)
            {
                return
                    (d1.cell1 < d2.cell1)
                    || (d1.cell1 == d2.cell1 && d1.cell2 < d2.cell2);
            };
        std::sort(nncData.begin(), nncData.end(), nncLess);
        auto candidate = nncData.begin();
        for (const auto& edit: editnncData) {
            auto printNncWarning =
                [](int c1, int c2)
                {
                    std::ostringstream sstr;
                    sstr << "Cannot edit NNC from " << c1 << " to " << c2
                         << " as it does not exist";
                    Opm::OpmLog::warning(sstr.str());
                };
            if (candidate == nncData.end()) {
                // no more NNCs left
                printNncWarning(edit.cell1, edit.cell2);
                continue;
            }
            if (candidate->cell1 != edit.cell1 || candidate->cell2 != edit.cell2) {
                candidate = std::lower_bound(candidate, nncData.end(), Opm::NNCdata(edit.cell1, edit.cell2, 0), nncLess);
                if (candidate == nncData.end()) {
                    // no more NNCs left
                    printNncWarning(edit.cell1, edit.cell2);
                    continue;
                }
            }
            auto firstCandidate = candidate;
            while (candidate != nncData.end()
                   && candidate->cell1 == edit.cell1
                   && candidate->cell2 == edit.cell2)
            {
                candidate->trans *= edit.trans;
                ++candidate;
            }
            // start with first match in next iteration to catch case where next
            // EDITNNC is for same pair.
            candidate = firstCandidate;
        }

        for (const auto& nncEntry : nnc.nncdata()) {
            auto c1 = nncEntry.cell1;
            auto c2 = nncEntry.cell2;
            auto low = cartesianToCompressed[c1];
            auto high = cartesianToCompressed[c2];

            if (low > high)
                std::swap(low, high);

            if (low == -1 && high == -1)
                // Silently discard as it is not between active cells
                continue;

            if (low == -1 || high == -1)
                OPM_THROW(std::logic_error, "NNC between active and inactive cells (" <<
                          low << " -> " << high);

            auto candidate = trans_.find(isId_(low, high));

            if (candidate == trans_.end())
                // This NNC is not resembled by the grid. Save it for later
                // processing with local cell values
                unprocessedNnc.push_back({c1, c2, nncEntry.trans});
            else {
                // NNC is represented by the grid and might be a neighboring connection
                // In this case the transmissibilty is added to the value already
                // set or computed.
                candidate->second += nncEntry.trans;
                processedNnc.push_back({c1, c2, nncEntry.trans});
            }
        }
        return make_tuple(processedNnc, unprocessedNnc);
    }

    /// \brief Multiplies the grid transmissibilities according to EDITNNC.
    void applyEditNncToGridTrans_(const ElementMapper& elementMapper,
                                  const std::vector<int>& globalToLocal)
    {
        const auto& editNnc = vanguard_.eclState().getInputEDITNNC();
        if (editNnc.empty())
            return;

        // editNnc is supposed to only reference non-neighboring connections and not
        // neighboring connections. Use all entries for scaling if there is an NNC.
        // variable nnc incremented in loop body.
        auto nnc = editNnc.data().begin();
        auto end = editNnc.data().end();
        while (nnc != end) {
            auto c1 = nnc->cell1;
            auto c2 = nnc->cell2;
            auto low = globalToLocal[c1];
            auto high = globalToLocal[c2];
            if (low > high)
                std::swap(low, high);

            auto candidate = trans_.find(isId_(low, high));
            if (candidate == trans_.end()) {
                std::ostringstream sstr;
                sstr << "Cannot edit NNC from " << c1 << " to " << c2
                     << " as it does not exist";
                Opm::OpmLog::warning(sstr.str());
                ++nnc;
            }
            else {
                // NNC exists
                while (nnc!= end && c1==nnc->cell1 && c2==nnc->cell2) {
                    candidate->second *= nnc->trans;
                    ++nnc;
                }
            }
        }
    }

    void extractPermeability_()
    {
        const auto isRadial = this->gridIsRadial_();

        const std::string k00_kw = isRadial ? "PERMR"   : "PERMX";
        const std::string k11_kw = isRadial ? "PERMTHT" : "PERMY";
        const std::string k22_kw =            "PERMZ";

        const auto& props = vanguard_.eclState().get3DProperties();

        unsigned numElem = vanguard_.gridView().size(/*codim=*/0);
        permeability_.resize(numElem);

        // read the intrinsic permeabilities from the eclState. Note that all arrays
        // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
        // simulation grid might remove a few elements. (e.g. because it is distributed
        // over several processes.)
        if (props.hasDeckDoubleGridProperty(k00_kw)) {
            const std::vector<double>& K00_Data =
                props.getDoubleGridProperty(k00_kw).getData();
            const std::vector<double>& K11_Data = props.hasDeckDoubleGridProperty(k11_kw)
                ? props.getDoubleGridProperty(k11_kw).getData() : K00_Data;
            const std::vector<double>& K22_Data = props.hasDeckDoubleGridProperty(k22_kw)
                ? props.getDoubleGridProperty(k22_kw).getData() : K00_Data;

            for (size_t dofIdx = 0; dofIdx < numElem; ++ dofIdx) {
                unsigned cartesianElemIdx = vanguard_.cartesianIndex(dofIdx);
                permeability_[dofIdx] = 0.0;
                permeability_[dofIdx][0][0] = K00_Data[cartesianElemIdx];
                permeability_[dofIdx][1][1] = K11_Data[cartesianElemIdx];
                permeability_[dofIdx][2][2] = K22_Data[cartesianElemIdx];
            }

            // for now we don't care about non-diagonal entries
        }
        else
            throw std::logic_error("Can't read the intrinsic permeability from the ecl state. "
                                   "(The PERM{X,Y,Z} keywords are missing)");
    }

    std::uint64_t isId_(std::uint32_t elemIdx1, std::uint32_t elemIdx2) const
    {
        std::uint32_t elemAIdx = std::min(elemIdx1, elemIdx2);
        std::uint64_t elemBIdx = std::max(elemIdx1, elemIdx2);

        return (elemBIdx<<elemIdxShift) + elemAIdx;
    }

    std::pair<std::uint32_t, std::uint32_t> isIdReverse_(const std::uint64_t& id) const
    {
        // Assigning an unsigned integer to a narrower type discards the most significant bits.
        // See "The C programming language", section A.6.2.
        // NOTE that the ordering of element A and B may have changed
        std::uint32_t elemAIdx = id;
        std::uint32_t elemBIdx = (id - elemAIdx) >> elemIdxShift;

        return std::make_pair(elemAIdx, elemBIdx);
    }

    std::uint64_t directionalIsId_(std::uint32_t elemIdx1, std::uint32_t elemIdx2) const
    {
        return (std::uint64_t(elemIdx1)<<elemIdxShift) + elemIdx2;
    }

    void computeHalfTrans_(Scalar& halfTrans,
                           const DimVector& areaNormal,
                           int faceIdx, // in the reference element that contains the intersection
                           const DimVector& distance,
                           const DimMatrix& perm) const
    {
        assert(faceIdx >= 0);
        unsigned dimIdx = faceIdx/2;
        assert(dimIdx < dimWorld);
        halfTrans = perm[dimIdx][dimIdx];

        Scalar val = 0;
        for (unsigned i = 0; i < areaNormal.size(); ++i)
            val += areaNormal[i]*distance[i];

        halfTrans *= std::abs(val);
        halfTrans /= distance.two_norm2();
    }

    template <class TransMultipliers>
    void computeFullTrans_(Scalar&& halfTrans1,
                           Scalar&& halfTrans2,
                           const EntityID& inside,
                           const unsigned int insideFaceIdx,
                           const EntityID& outside,
                           const unsigned int outsideFaceIdx,
                           const TransMultipliers& transMult,
                           const std::vector<double>& ntg)
    {
        applyNtg_(halfTrans1, insideFaceIdx, inside.cartElemIdx, ntg);
        applyNtg_(halfTrans2, outsideFaceIdx, outside.cartElemIdx, ntg);

        Scalar trans;
        if (std::abs(halfTrans1) < 1e-30 || std::abs(halfTrans2) < 1e-30)
            // avoid division by zero
            trans = 0.0;
        else
            trans = 1.0 / (1.0/halfTrans1 + 1.0/halfTrans2);

        // apply the full face transmissibility multipliers
        // for the inside ...
        //
        // The MULTZ needs special case if the option is ALL
        // Then the smallest multiplier is applied.
        // Default is to apply the top and bottom multiplier
        bool useSmallestMultiplier = eclGrid.getMultzOption() == Opm::PinchMode::ModeEnum::ALL;
        if (useSmallestMultiplier)
            applyAllZMultipliers_(trans, insideFaceIdx, inside.cartElemIdx, outside.cartElemIdx, transMult, cartDims);
        else
            applyMultipliers_(trans, insideFaceIdx, inside.cartElemIdx, transMult);
        // ... and outside elements
        applyMultipliers_(trans, outsideFaceIdx, outside.cartElemIdx, transMult);

        // apply the region multipliers (cf. the MULTREGT keyword)
        Opm::FaceDir::DirEnum faceDir;
        switch (insideFaceIdx) {
        case 0:
        case 1:
            faceDir = Opm::FaceDir::XPlus;
            break;

        case 2:
        case 3:
            faceDir = Opm::FaceDir::YPlus;
            break;

        case 4:
        case 5:
            faceDir = Opm::FaceDir::ZPlus;
            break;

        default:
            OPM_THROW(std::logic_error, "Could not determine a face direction");
        }

        trans *= transMult.getRegionMultiplier(inside.cartElemIdx,
                                               outside.cartElemIdx,
                                               faceDir);

        trans_[isId_(inside.elemIdx, outside.elemIdx)] = trans;
    }

    DimVector distanceVector_(const DimVector& center,
                              int faceIdx, // in the reference element that contains the intersection
                              unsigned elemIdx,
                              const std::array<std::vector<DimVector>, dimWorld>& axisCentroids) const
    {
        assert(faceIdx >= 0);
        unsigned dimIdx = faceIdx/2;
        assert(dimIdx < dimWorld);
        DimVector x = center;
        x -= axisCentroids[dimIdx][elemIdx];

        return x;
    }

    void applyMultipliers_(Scalar& trans,
                           unsigned faceIdx,
                           unsigned cartElemIdx,
                           const Opm::TransMult& transMult) const
    {
        // apply multiplyer for the transmissibility of the face. (the
        // face index is the index of the reference-element face which
        // contains the intersection of interest.)
        switch (faceIdx) {
        case 0: // left
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::XMinus);
            break;
        case 1: // right
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::XPlus);
            break;

        case 2: // front
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::YMinus);
            break;
        case 3: // back
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::YPlus);
            break;

        case 4: // bottom
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::ZMinus);
            break;
        case 5: // top
            trans *= transMult.getMultiplier(cartElemIdx, Opm::FaceDir::ZPlus);
            break;
        }
    }

    void applyNtg_(Scalar& trans,
                   unsigned faceIdx,
                   unsigned cartElemIdx,
                   const std::vector<double>& ntg) const
    {
        // apply multiplyer for the transmissibility of the face. (the
        // face index is the index of the reference-element face which
        // contains the intersection of interest.)
        switch (faceIdx) {
        case 0: // left
            trans *= ntg[cartElemIdx];
            break;
        case 1: // right
            trans *= ntg[cartElemIdx];
            break;

        case 2: // front
            trans *= ntg[cartElemIdx];
            break;
        case 3: // back
            trans *= ntg[cartElemIdx];
            break;

            // NTG does not apply to top and bottom faces
        }
    }

    void minPvFillNtg_(std::vector<double>& averageNtg) const
    {
        // compute volume weighted arithmetic average of NTG for
        // cells merged as an result of minpv.
        const auto& eclState = vanguard_.eclState();
        const auto& eclGrid = eclState.getInputGrid();

        const std::vector<double>& ntg =
            eclState.get3DProperties().getDoubleGridProperty("NTG").getData();

        averageNtg = ntg;

        bool opmfil = eclGrid.getMinpvMode() == Opm::MinpvMode::OpmFIL;

        // just return the unmodified ntg if opmfil is not used
        if (!opmfil)
            return;

        const auto& porv = eclState.get3DProperties().getDoubleGridProperty("PORV").getData();
        const auto& actnum = eclState.get3DProperties().getIntGridProperty("ACTNUM").getData();
        const auto& cartMapper = vanguard_.cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        assert(dimWorld > 1);
        const size_t nxny = cartDims[0] * cartDims[1];
        for (size_t cartesianCellIdx = 0; cartesianCellIdx < ntg.size(); ++cartesianCellIdx) {
            // use the original ntg values for the inactive cells
            if (!actnum[cartesianCellIdx])
                continue;

            // Average properties as long as there exist cells above
            // that has pore volume less than the MINPV threshold
            const double cellVolume = eclGrid.getCellVolume(cartesianCellIdx);
            double ntgCellVolume = ntg[cartesianCellIdx] * cellVolume;
            double totalCellVolume = cellVolume;
            int cartesianCellIdxAbove = cartesianCellIdx - nxny;
            while (cartesianCellIdxAbove >= 0 &&
                   actnum[cartesianCellIdxAbove] > 0 &&
                   porv[cartesianCellIdxAbove] < eclGrid.getMinpvVector()[cartesianCellIdxAbove])
            {
                // Volume weighted arithmetic average of NTG
                const double cellAboveVolume = eclGrid.getCellVolume(cartesianCellIdxAbove);
                totalCellVolume += cellAboveVolume;
                ntgCellVolume += ntg[cartesianCellIdxAbove]*cellAboveVolume;
                cartesianCellIdxAbove -= nxny;
            }

            averageNtg[cartesianCellIdx] = ntgCellVolume / totalCellVolume;
        }
    }


    bool gridIsRadial_() const
    {
        return this->gridIsRadial_(typename std::is_same<Grid, Dune::CpGrid>::type{});
    }


    bool gridIsRadial_(std::false_type /* !CpGrid */) const
    {
        return false;
    }


    bool gridIsRadial_(std::true_type /* CpGrid */) const
    {
        return this->vanguard_.eclState().getInputGrid().isRadial();
    }


    void computeTrans_OldTran_Radial_(std::false_type /* !CpGrid */,
                                      const GridView& /* gridView */,
                                      const ElementMapper& /* elemMapper */,
                                      const unsigned int /* numElements*/,
                                      const std::vector<double>& /* ntg */)
    {
    }


    void computeTrans_OldTran_Radial_(std::true_type /* CpGrid */,
                                      const GridView& gridView,
                                      const ElementMapper& elemMapper,
                                      const unsigned int numElements,
                                      const std::vector<double>& ntg)
    {
        OldTran_Radial calc(*this);

        calc.computeTrans(gridView, elemMapper, numElements, ntg);
    }


    const Vanguard& vanguard_;
    Scalar transmissibilityThreshold_;
    std::vector<DimMatrix> permeability_;
    std::unordered_map<std::uint64_t, Scalar> trans_;
    std::map<std::pair<unsigned, unsigned>, Scalar> transBoundary_;
    std::map<std::pair<unsigned, unsigned>, Scalar> thermalHalfTransBoundary_;
    Opm::ConditionalStorage<enableEnergy,
                            std::unordered_map<std::uint64_t, Scalar> > thermalHalfTrans_;
};

} // namespace Ewoms

#endif
