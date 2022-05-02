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
 * \copydoc Opm::Co2InjectionTestProblem
 */
#ifndef EWOMS_CO2_INJECTION_TEST_PROBLEM_HH
#define EWOMS_CO2_INJECTION_TEST_PROBLEM_HH

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
template <class TypeTag>
class Co2InjectionTestProblem;

namespace Co2InjectionTest {
#include <opm/material/components/co2tables.inc>
}
} // namespace Opm

namespace Opm::Properties {

namespace TTag {
struct Co2InjectionTestBaseProblem {};
}

// declare the CO2 injection problem specific property tags
template <class TypeTag, class MyTypeTag>
struct Temperature { using type = UndefinedProperty; };
template <class TypeTag, class MyTypeTag>
struct SimulationName { using type = UndefinedProperty; };
template <class TypeTag, class MyTypeTag>
struct Inflowrate { using type = UndefinedProperty;};

template <class TypeTag, class MyTypeTag>
struct Initialpressure { using type = UndefinedProperty;};

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
};

// Set the material Law
template <class TypeTag>
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
    //using type = Opm::EffToAbsLaw<EffMaterialLaw>;
     using type = EffMaterialLaw;
};

// Write the Newton convergence behavior to disk?
template<class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::Co2InjectionTestBaseProblem> { 
static constexpr bool value = false; };

// Enable gravity
template <class TypeTag>
struct EnableGravity<TypeTag, TTag::Co2InjectionTestBaseProblem> 
{ static constexpr bool value = true;
};

template <class TypeTag>
struct Temperature<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 273.15 + TEMPERATURE;
};

template <class TypeTag>
struct Inflowrate<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = INFLOW_RATE;
};
template <class TypeTag>
struct Initialpressure<TypeTag, TTag::Co2InjectionTestBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = MIN_PRES;
};

template<class TypeTag>
struct SimulationName<TypeTag, TTag::Co2InjectionTestBaseProblem> { 
    static constexpr auto value = "co2injectiontest"; 
};

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
//struct GridFile<TypeTag, TTag::Co2InjectionTestBaseProblem> 
//{ static constexpr auto value = "data/co2injection.dgf"; };
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
 * \brief Problem that mimics the problem in compositional branch
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
    using Indices = GetPropType<TypeTag, Properties::Indices>;
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
    
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { BrineIdx = FluidSystem::BrineIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + CO2Idx };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    //enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableGravity = getPropValue<TypeTag, Properties::EnableGravity>() };
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
    void initGravity()
    {
        gravity_ = 0.0;
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dimWorld - 1] =  (-9.81);
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();
        // initialize fixed parameters; temperature, permeability, porosity
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        K_ = this->toDimMatrix_(PERMEABILITY * 1.e-15);
        porosity_ = POROSITY;

        // initialize gravity
        initGravity();
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Inflowrate, 
                            "The inflow rate [?] on the left boundary of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Initialpressure, 
                            "The initial pressure [Pa s] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, 
                            std::string, 
                            SimulationName,
                             "The name of the simulation used for the output "
                             "files");
    }

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
        Scalar tol = this->model().newtonMethod().tolerance() * 1e5;
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
    }

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Evaluation, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);
        values.assignNaive(fs);
    }

    // Constant temperature
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED, unsigned spaceIdx OPM_UNUSED, unsigned timeIdx OPM_UNUSED) const
    {
        return temperature_;
    }

    // Constant permeability
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context OPM_UNUSED,
                                           unsigned spaceIdx OPM_UNUSED,
                                           unsigned timeIdx OPM_UNUSED) const
    {
        return K_;
    }

    // Constant porosity
    template <class Context>
    Scalar porosity(const Context& context OPM_UNUSED, unsigned spaceIdx OPM_UNUSED, unsigned timeIdx OPM_UNUSED) const
    {
        return porosity_;
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
    void boundary(BoundaryRateVector& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        // const Scalar eps = std::numeric_limits<double>::epsilon();
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos)) {
            Scalar inflowrate = EWOMS_GET_PARAM(TypeTag, Scalar, Inflowrate);
            RateVector massRate(0.0);
            massRate = 0.0;
            massRate[contiCO2EqIdx] = inflowrate; // kg / (m^2 * s)

            values.setMassRate(massRate);
        } else if (onRightBoundary_(pos)) {
            Opm::CompositionalFluidState<Evaluation, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        } else
            values.setNoFlow(); // closed on top and bottom
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
    {
        return pos[XDIM] > this->boundingBoxMax()[XDIM] - 1e-6;
    }

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
        // get capillary pressure
        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);        
        //////
        // set temperature
        //////
        fs.setTemperature(temperature_);
        //////
        // set saturations
        //////
        fs.setSaturation(FluidSystem::liquidPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);
        //////
        // set pressures
        //////
        // pressure; set simple hydrostatic pressure initially.
        // OBS: If horizontal (NZ = 1), then h = Z_SIZE/2 since boundingBoxMax is cell edge
        Scalar init_pressure = EWOMS_GET_PARAM(TypeTag, Scalar, Initialpressure);
        bool enable_gravity = EWOMS_GET_PARAM(TypeTag, bool, EnableGravity);
        Scalar p_init;        
         if (enable_gravity == true) {
             Scalar densityW = Brine::liquidDensity(temperature_, Scalar(init_pressure));
             const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
             Scalar h = this->boundingBoxMax()[ZDIM] - pos[ZDIM];
             p_init = (init_pressure * 1e5) + densityW * h * (-this->gravity_[1]);
         } else
            p_init = init_pressure * 1e5;
        fs.setPressure(liquidPhaseIdx, p_init);
        fs.setPressure(gasPhaseIdx, p_init);


        //////
        // set composition of the liquid phase
        //////
        fs.setMoleFraction(liquidPhaseIdx, CO2Idx, 0);
        fs.setMoleFraction(liquidPhaseIdx, BrineIdx,1);

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        using CFRP = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/liquidPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/true);
    }


    DimMatrix K_;
    Scalar porosity_;
    Scalar temperature_;
    MaterialLawParams mat_;
    DimVector gravity_;
};
} // namespace Opm

#endif