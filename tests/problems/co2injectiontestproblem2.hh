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
 * \copydoc Opm::Co2InjectionProblem
 */
#ifndef EWOMS_CO2_INJECTION_PROBLEM_HH
#define EWOMS_CO2_INJECTION_PROBLEM_HH

#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/simulators/linalg/parallelamgbackend.hh>
#include <opm/models/io/structuredgridvanguard.hh>
#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>

#include <opm/material/fluidsystems/compositionalfluid/twophasefluidsystem.hh>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <iostream>
#include <string>

namespace Opm {
//! \cond SKIP_THIS
template <class TypeTag>
class Co2InjectionTestProblem;

namespace Co2InjectionTest {
#include <opm/material/components/co2tables.inc>
}
//! \endcond
}

namespace Opm::Properties {

namespace TTag {
struct Co2InjectionTestBaseProblem {};
}

// declare the CO2 injection problem specific property tags
template<class TypeTag, class MyTypeTag>
struct FluidSystemPressureLow { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemPressureHigh { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemNumPressure { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemTemperatureLow { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemTemperatureHigh { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemNumTemperature { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct MaxDepth { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct Temperature { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct SimulationName { using type = UndefinedProperty; };
template <class TypeTag, class MyTypeTag>
struct Inflowrate { using type = UndefinedProperty;};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Co2InjectionTestBaseProblem> 
{ using type = Dune::YaspGrid<3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Co2InjectionTestBaseProblem>
{ using type = Opm::Co2InjectionTestProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = Opm::Co2InjectionTest::CO2Tables;

public:
    using type = Opm::BrineCO2FluidSystem<Scalar, CO2Tables>;
    //using type = Opm::H2ON2FluidSystem<Scalar, /*useComplexRelations=*/false>;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                               /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx>;

    // define the material law which is parameterized by effective
    // saturations
    using EffMaterialLaw = Opm::NullMaterial<Traits>;

public:
    // define the material law parameterized by absolute saturations
    using type = EffMaterialLaw;
};

// Use the algebraic multi-grid linear solver for this problem
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::Co2InjectionTestBaseProblem> { using type = TTag::ParallelAmgLinearSolver; };

// Write the Newton convergence behavior to disk?
template<class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::Co2InjectionTestBaseProblem> { static constexpr bool value = false; };

// Enable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::Co2InjectionTestBaseProblem> { static constexpr bool value = true; };

// set the defaults for the problem specific properties
template<class TypeTag>
struct FluidSystemPressureLow<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3e7;
};
template<class TypeTag>
struct FluidSystemPressureHigh<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 4e7;
};
template<class TypeTag>
struct FluidSystemNumPressure<TypeTag, TTag::Co2InjectionTestBaseProblem> { static constexpr int value = 100; };
template<class TypeTag>
struct FluidSystemTemperatureLow<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 290;
};
template<class TypeTag>
struct FluidSystemTemperatureHigh<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 500;
};
template<class TypeTag>
struct FluidSystemNumTemperature<TypeTag, TTag::Co2InjectionTestBaseProblem> { static constexpr int value = 100; };

template<class TypeTag>
struct MaxDepth<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 2500;
};
template<class TypeTag>
struct Temperature<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 293.15;
};
template <class TypeTag>
struct Inflowrate<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = INFLOW_RATE;
};
template<class TypeTag>
struct SimulationName<TypeTag, TTag::Co2InjectionTestBaseProblem> { static constexpr auto value = "co2injection"; };

// The default for the end time of the simulation
template <class TypeTag>
struct EndTime<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = SIM_TIME * 24. * 60. * 60.;
};

// The default for the initial time step size of the simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::Co2InjectionTestBaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 30;
};

template <class TypeTag>
struct LinearSolverTolerance<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};

template <class TypeTag>
struct LinearSolverAbsTolerance<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.;
};
template <class TypeTag>
struct NewtonTolerance<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};

template <class TypeTag>
struct MaxTimeStepSize<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 60 * 60;
};

 template <class TypeTag>
 struct NewtonMaxIterations<TypeTag, TTag::Co2InjectionTestBaseProblem> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = 10;
 };

 template <class TypeTag>
 struct NewtonTargetIterations<TypeTag, TTag::Co2InjectionTestBaseProblem> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = 6;
 };

// The default DGF file to load
//template<class TypeTag>
//struct GridFile<TypeTag, TTag::Co2InjectionTestBaseProblem> { static constexpr auto value = "data/co2injection.dgf"; };
// mesh grid
template <class TypeTag>
struct Vanguard<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = Opm::StructuredGridVanguard<TypeTag>;
};

template <class TypeTag>
struct DomainSizeX<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = X_SIZE; // meter
};

template <class TypeTag>
struct CellsX<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    static constexpr int value = NX;
};

template <class TypeTag>
struct DomainSizeY<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = Y_SIZE; // meter
};

template <class TypeTag>
struct CellsY<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    static constexpr int value = NY;
};

 template <class TypeTag>
 struct DomainSizeZ<TypeTag, TTag::Co2InjectionTestBaseProblem> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = Z_SIZE; // meter
 };

 template <class TypeTag>
 struct CellsZ<TypeTag, TTag::Co2InjectionTestBaseProblem> {
     static constexpr int value = NZ;
 };

template <class TypeTag>
struct EnableEnergy<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    static constexpr bool value = false;
};

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief Problem where \f$CO_2\f$ is injected under a low permeable
 *        layer at a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, one
 * which is moderately permeable (\f$K = 10^{-12}\;m^2\f$) for \f$ y >
 * 22\; m\f$ and one with a lower intrinsic permeablility (\f$
 * K=10^{-13}\;m^2\f$) in the rest of the domain.
 *
 * \f$CO_2\f$ gets injected by means of a forced-flow boundary
 * condition into water-filled aquifer, which is situated 2700m below
 * sea level, at the lower-right boundary (\f$5m<y<15m\f$) and
 * migrates upwards due to buoyancy. It accumulates and eventually
 * enters the lower permeable aquitard.
 *
 * The boundary conditions applied by this problem are no-flow
 * conditions on the top bottom and right boundaries and a free-flow
 * boundary condition on the left. For the free-flow condition,
 * hydrostatic pressure is assumed.
 */
template <class TypeTag>
class Co2InjectionTestProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    enum { numPhases = FluidSystem::numPhases };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { BrineIdx = FluidSystem::BrineIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + CO2Idx };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;

    using H2O = typename Opm::H2O<Scalar>;
    using Brine = typename Opm::Brine<Scalar, H2O>;    
    using Toolbox = Opm::MathToolbox<Evaluation>;
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    const unsigned XDIM = 0;
    const unsigned YDIM = 1;
    const unsigned ZDIM = 2;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Co2InjectionTestProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1e-6;

        temperatureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow);
        temperatureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh);
        nTemperature_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumTemperature);

        pressureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureLow);
        pressureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureHigh);
        nPressure_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumPressure);

        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        K_ = this->toDimMatrix_(PERMEABILITY * 1.e-15);


        // initialize the tables of the fluid system
        // FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        fineLayerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-13);
        coarseK_ = this->toDimMatrix_(1e-12);

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow,
                             "The lower temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh,
                             "The upper temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumTemperature,
                             "The number of intervals between the lower and "
                             "upper temperature");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureLow,
                             "The lower pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureHigh,
                             "The upper pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumPressure,
                             "The number of intervals between the lower and "
                             "upper pressure");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Inflowrate, 
                            "The inflow rate [?] on the left boundary of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxDepth,
                             "The maximum depth [m] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName)
            << "_" << Model::name();
        if (getPropValue<TypeTag, Properties::EnableEnergy>())
            oss << "_ni";
        oss << "_" << Model::discretizationName();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        Scalar tol = this->model().newtonMethod().tolerance()*1e5;
        this->model().checkConservativeness(tol);

        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, /*phaseIdx=*/0);
        this->model().globalPhaseStorage(storageG, /*phaseIdx=*/1);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: liquid=[" << storageL << "]"
                      << " gas=[" << storageG << "]\n" << std::flush;
        }
#endif // NDEBUG
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (inHighTemperatureRegion_(pos))
            return temperature_ + 100;
        return temperature_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        return coarsePorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context OPM_UNUSED,
                                               unsigned spaceIdx OPM_UNUSED,
                                               unsigned timeIdx OPM_UNUSED) const
    {
        return this->mat_;
    }

    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (onLeftBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.checkDefined();

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onInlet_(pos)) {
            RateVector massRate(0.0);
            massRate[contiCO2EqIdx] = -1e-3; // [kg/(m^3 s)]

            using FluidState = Opm::ImmiscibleFluidState<Scalar, FluidSystem>;
            FluidState fs;
            fs.setSaturation(gasPhaseIdx, 1.0);
            const auto& pg =
                context.intensiveQuantities(spaceIdx, timeIdx).fluidState().pressure(gasPhaseIdx);
            fs.setPressure(gasPhaseIdx, Toolbox::value(pg));
            fs.setTemperature(temperature(context, spaceIdx, timeIdx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            paramCache.updatePhase(fs, gasPhaseIdx);
            Scalar h = FluidSystem::template enthalpy<FluidState, Scalar>(fs, paramCache, gasPhaseIdx);

            // impose an forced inflow boundary condition for pure CO2
            values.setMassRate(massRate);
            values.setEnthalpyRate(massRate[contiCO2EqIdx] * h);
        }
        else
            // no flow on top and bottom
            values.setNoFlow();
    }



    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        // const auto& matParams = this->materialLawParams(context, spaceIdx,
        // timeIdx);
        // values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        values.assignNaive(fs);
    }


   // No source terms
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    { 
        rate = Scalar(0.0); 
    }

private:
    bool onLeftBoundary_(const GlobalPosition& pos) const
    {
        return pos[XDIM] < 1e-6;
    }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    { return onRightBoundary_(pos) && (5 < pos[1]) && (pos[1] < 15); }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    {
        return pos[ZDIM] < 1e-6;
    }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    {
        return pos[ZDIM] > this->boundingBoxMax()[ZDIM] - 1e-6;
    }

    bool aboveMiddle_(const GlobalPosition& pos) const
    {
        return pos[ZDIM] >= (this->boundingBoxMax()[ZDIM] + this->boundingBoxMin()[ZDIM]) / 2;
    }

    bool leftMiddle(const GlobalPosition& pos) const
    {
        return pos[XDIM] > (this->boundingBoxMax()[XDIM] + this->boundingBoxMin()[XDIM]) / 2;
    }

   /*!
     * \copydoc FvBaseProblem::initial
     */    
    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        //////
        // set temperature
        //////
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));

        //////
        // set saturations
        //////
        fs.setSaturation(FluidSystem::liquidPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

        //////
        // set pressures
        //////
        Scalar densityL = FluidSystem::Brine::liquidDensity(temperature_, Scalar(1e5));
        Scalar depth = maxDepth_ - pos[dim - 1];
        Scalar pl = 1e5 - densityL * this->gravity()[dim - 1] * depth;

        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(liquidPhaseIdx, pl + (pC[liquidPhaseIdx] - pC[liquidPhaseIdx]));
        fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

        //////
        // set composition of the liquid phase
        //////
        fs.setMoleFraction(liquidPhaseIdx, CO2Idx, 0.005);
        fs.setMoleFraction(liquidPhaseIdx, BrineIdx,
                           1.0 - fs.moleFraction(liquidPhaseIdx, CO2Idx));

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        using CFRP = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/liquidPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/true);
    }


    bool inHighTemperatureRegion_(const GlobalPosition& pos) const
    { return (pos[0] > 20) && (pos[0] < 30) && (pos[1] > 5) && (pos[1] < 35); }

    bool isFineMaterial_(const GlobalPosition& pos) const
    { return pos[dim - 1] > fineLayerBottom_; }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar fineLayerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

        DimMatrix K_;
    Scalar porosity_;
    Scalar temperature_;
    MaterialLawParams mat_;
    DimVector gravity_;
    Scalar maxDepth_;
    Scalar eps_;

    unsigned nTemperature_;
    unsigned nPressure_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
};
} // namespace Opm

#endif
