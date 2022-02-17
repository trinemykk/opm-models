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
 * \copydoc Opm::SmeaheiaProblem
 */
#ifndef SMEAHEIA_PROBLEM_HH
#define SMEAHEIA_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <iostream>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/immiscible/immisciblemodel.hh>
//#include <opm/models/common/transfluxmodule.hh>
//#include <opm/multidomain/multidimproperties.hh>
#include <opm/simulators/linalg/parallelamgbackend.hh>
//#include <opm/paper/pressureDependentIntensiveQuantities.hh>//
#include <sstream>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>

namespace Opm {
//! \cond SKIP_THIS
template <class TypeTag>
class SmeaheiaProblem;

namespace Co2Injection {
#include <opm/material/components/co2tables.inc>
}
//! \endcond
}  // namespace Opm

namespace Opm::Properties {

namespace TTag {
struct SmeaheiaBaseProblem {};
}  // namespace TTag

template <class TypeTag, class MyTypeTag>
struct InjectionRate {
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct ReservoirLayerTop {
    using type = UndefinedProperty;
};

// template <class TypeTag, class MyTypeTag>
// struct FracturePermeability {
//     using type = UndefinedProperty;
// };

// template <class TypeTag, class MyTypeTag>
// struct EffectiveDeformationBandPermeability {
//     using type = UndefinedProperty;
// };

template <class TypeTag, class MyTypeTag>
struct PressureDependentPermeability {
    using type = UndefinedProperty;
};

// template <class TypeTag, class MyTypeTag>
// struct IrreversibleDeformation {
//     using type = UndefinedProperty;
// };

// template <class TypeTag, class MyTypeTag>
// struct DilationAngle {
//     using type = UndefinedProperty;
// };

template <class TypeTag, class MyTypeTag>
struct CsvFilePath {
    using type = UndefinedProperty;
};

// declare the CO2 injection problem specific property tags
template <class TypeTag, class MyTypeTag>
struct FluidSystemPressureLow {
    using type = UndefinedProperty;
};
template <class TypeTag, class MyTypeTag>
struct FluidSystemPressureHigh {
    using type = UndefinedProperty;
};
template <class TypeTag, class MyTypeTag>
struct FluidSystemNumPressure {
    using type = UndefinedProperty;
};
template <class TypeTag, class MyTypeTag>
struct FluidSystemTemperatureLow {
    using type = UndefinedProperty;
};
template <class TypeTag, class MyTypeTag>
struct FluidSystemTemperatureHigh {
    using type = UndefinedProperty;
};
template <class TypeTag, class MyTypeTag>
struct FluidSystemNumTemperature {
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct MaxDepth {
    using type = UndefinedProperty;
};
template <class TypeTag, class MyTypeTag>
struct Temperature {
    using type = UndefinedProperty;
};
template <class TypeTag, class MyTypeTag>
struct SimulationName {
    using type = UndefinedProperty;
};

// template <class TypeTag, class MyTypeTag>
// struct EnableDeformationBand {
//     using type = UndefinedProperty;
// };

// template <class TypeTag, class MyTypeTag>
// struct EnableFaultBoundary {
//     using type = UndefinedProperty;
// };

// Enable adaptivity
// template <class TypeTag>
// struct EnableGridAdaptation<TypeTag, TTag::SmeaheiaBaseProblem> {
//     static constexpr bool value = false;
// };
// No async vtk with adaptivity
// template <class TypeTag>
// struct EnableAsyncVtkOutput<TypeTag, TTag::SmeaheiaBaseProblem> {
//     static constexpr bool value = true;
// };

// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = Opm::SmeaheiaProblem<TypeTag>;
};

// // Set the problem property
// template <class TypeTag>
// struct FluxModule<TypeTag, TTag::SmeaheiaBaseProblem> {
//     using type = TransFluxModule<TypeTag>;
// };

// // Set the problem property
// template <class TypeTag>
// struct IntensiveQuantities<TypeTag, TTag::SmeaheiaBaseProblem> {
//     using type = Opm::PressureDependentIntensiveQuantities<TypeTag>;
// };


// Set fluid configuration
template <class TypeTag>
struct FluidSystem<TypeTag, TTag::SmeaheiaBaseProblem> {
   private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = Opm::Co2Injection::CO2Tables;

   public:
    using type = Opm::BrineCO2FluidSystem<Scalar, CO2Tables>;
    // using type = Opm::H2ON2FluidSystem<Scalar, /*useComplexRelations=*/false>;
};

// Set the material Law
template <class TypeTag>
struct MaterialLaw<TypeTag, TTag::SmeaheiaBaseProblem> {
   private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/liquidPhaseIdx,
                                               /*nonWettingPhaseIdx=*/gasPhaseIdx>;

    // define the material law which is parameterized by effective
    // saturations
    using EffMaterialLaw = Opm::RegularizedBrooksCorey<Traits>;

   public:
    // define the material law parameterized by absolute saturations
    using type = Opm::EffToAbsLaw<EffMaterialLaw>;
};

template <class TypeTag>
struct PressureDependentPermeability<TypeTag, TTag::SmeaheiaBaseProblem> {
    static constexpr bool value = false;
};

// template <class TypeTag>
// struct IrreversibleDeformation<TypeTag, TTag::SmeaheiaBaseProblem> {
//     static constexpr bool value = false;
// };

// template <class TypeTag>
// struct DilationAngle<TypeTag, TTag::SmeaheiaBaseProblem> {
//     using type = GetPropType<TypeTag, Scalar>;
//     static constexpr type value = 5.0 / 180.0 * 3.1415;
// };

template <class TypeTag>
struct InjectionRate<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.0;
};

template <class TypeTag>
struct ReservoirLayerTop<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 120.0;
};

// template <class TypeTag>
// struct FracturePermeability<TypeTag, TTag::SmeaheiaBaseProblem> {
//     using type = GetPropType<TypeTag, Scalar>;
//     static constexpr type value = 1e-16;
// };

// template <class TypeTag>
// struct EffectiveDeformationBandPermeability<TypeTag, TTag::SmeaheiaBaseProblem> {
//     using type = GetPropType<TypeTag, Scalar>;
//     static constexpr type value = 1e-1;
// };

// Set the thermal conduction law
template <class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::SmeaheiaBaseProblem> {
   private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

   public:
    // define the material law parameterized by absolute saturations
    using type = Opm::SomertonThermalConductionLaw<FluidSystem, Scalar>;
};

// set the energy storage law for the solid phase
template <class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = Opm::ConstantSolidHeatCapLaw<GetPropType<TypeTag, Properties::Scalar>>;
};

// Use the algebraic multi-grid linear solver for this problem
template <class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = TTag::ParallelAmgLinearSolver;
};

// Write the Newton convergence behavior to disk?
template <class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::SmeaheiaBaseProblem> {
    static constexpr bool value = false;
};

// By default, include the intrinsic permeability tensor to the VTK output files
template <class TypeTag>
struct VtkWriteIntrinsicPermeabilities<TypeTag, TTag::SmeaheiaBaseProblem> {
    static constexpr bool value = true;
};

// Enable gravity
template <class TypeTag>
struct EnableGravity<TypeTag, TTag::SmeaheiaBaseProblem> {
    static constexpr bool value = true;
};

template <class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = TTag::EcfvDiscretization;
};

// use automatic differentiation for this simulator
template <class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = TTag::AutoDiffLocalLinearizer;
};

// set the defaults for the problem specific properties
template <class TypeTag>
struct FluidSystemPressureLow<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e5;
};
template <class TypeTag>
struct FluidSystemPressureHigh<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e8;
};
template <class TypeTag>
struct FluidSystemNumPressure<TypeTag, TTag::SmeaheiaBaseProblem> {
    static constexpr int value = 100;
};
template <class TypeTag>
struct FluidSystemTemperatureLow<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 280;
};
template <class TypeTag>
struct FluidSystemTemperatureHigh<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 400;
};
template <class TypeTag>
struct FluidSystemNumTemperature<TypeTag, TTag::SmeaheiaBaseProblem> {
    static constexpr int value = 50;
};

template <class TypeTag>
struct MaxDepth<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.0;
};
template <class TypeTag>
struct Temperature<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 363.15;
};

template <class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::SmeaheiaBaseProblem> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct SimulationName<TypeTag, TTag::SmeaheiaBaseProblem> {
    static constexpr auto value = "Smeaheia";
};


// The default for the end time of the simulation
template <class TypeTag>
struct EndTime<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3 * 365 * 24 * 60 * 60;
};

// The default for the initial time step size of the simulation
template <class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0 * 24 * 60 * 60;
};

// The default for the initial time step size of the simulation
template <class TypeTag>
struct MaxTimeStepSize<TypeTag, TTag::SmeaheiaBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 30 * 24 * 60 * 60;
};

// template <class TypeTag>
// struct EnableDeformationBand<TypeTag, TTag::SmeaheiaBaseProblem> {
//     static constexpr bool value = true;
// };

// template <class TypeTag>
// struct EnableFaultBoundary<TypeTag, TTag::SmeaheiaBaseProblem> {
//     static constexpr bool value = false;
// };

}  // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems TODO
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
class SmeaheiaProblem 
    : public GetPropType<TypeTag, Properties::BaseProblem>{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using BoundaryContext = GetPropType<TypeTag, Properties::BoundaryContext>;

    //enum { dim = getPropValue<TypeTag, Properties::DomainDim>() };
    //enum { dimWorld = getPropValue<TypeTag, Properties::WorldDim>() };
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
    using ThermalConductionLaw = GetPropType<TypeTag, Properties::ThermalConductionLaw>;
    using SolidEnergyLawParams = GetPropType<TypeTag, Properties::SolidEnergyLawParams>;
    using ThermalConductionLawParams = typename ThermalConductionLaw::Params;

    using Toolbox = Opm::MathToolbox<Evaluation>;
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;
    using FluxExtensiveQuantities = typename FluxModule::FluxExtensiveQuantities;

    //constexpr Scalar M_PI = 3.14159265358979323846;

   public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    SmeaheiaProblem(Simulator& simulator) : ParentType(simulator) {}
    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit() {
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
        
        csvFilePath_ = EWOMS_GET_PARAM(TypeTag, std::string, CsvFilePath);
        // initialize the tables of the fluid system
        // FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        reservoirLayerBottom_ = 0.0;
        reservoirLayerTop_ = EWOMS_GET_PARAM(TypeTag, Scalar, ReservoirLayerTop);

        //frac_density_ = .2; // Number of frac  / domain length

        // intrinsic permeabilities
        reservoirK_ = this->toDimMatrix_(1e-13);
		reservoirK_[dimWorld - 1][dimWorld - 1] *= 0.1;
		//fractureK_ = this->toDimMatrix_(EWOMS_GET_PARAM(TypeTag, Scalar, FracturePermeability));
		//fractureK_ = this->toDimMatrix_(1e-16);

        //aperture0_ = 1.0;

        // porosities
        //fracturePorosity_ = 0.1; //1.0;//1e-1 / 19.607843137254903;  // 1.5e-8;
        //reservoirPorosity_ = 0.17;

        loadPermeability_();
        loadPorosity_();
        // Dilation angle
        //dilationAngle_ = EWOMS_GET_PARAM(TypeTag, Scalar, DilationAngle);

        // residual saturations
        reservoirMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.2);
        reservoirMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);
        //fractureMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.2);
       // fractureMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the Brooks-Corey law
        reservoirMaterialParams_.setEntryPressure(2.5e3);
        reservoirMaterialParams_.setLambda(2.8);
       // fractureMaterialParams_.setEntryPressure(2.5e6);

       // fractureMaterialParams_.setLambda(2.8);
        reservoirMaterialParams_.finalize();
       // fractureMaterialParams_.finalize();

        // parameters for the somerton law thermal conduction
      //  computeThermalCondParams_(fineThermalCondParams_, fracturePorosity_);
        computeThermalCondParams_(coarseThermalCondParams_, reservoirPorosity_[0]);

        // assume constant heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(
            790.0       // specific heat capacity of granite [J / (kg K)]
            * 2700.0);  // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();
        injectionRate_ = EWOMS_GET_PARAM(TypeTag, Scalar, InjectionRate);

        // const auto& idxSet = this->simulator().gridView().indexSet();
        // irreversibleDilation_.resize(this->simulator().gridView().size(0));
        // for (auto& e : elements(this->simulator().gridView())) {
        //     irreversibleDilation_[idxSet.index(e)] = 0.0;
        // }

        std::ostringstream fileName;
        fileName << csvFilePath_ << getPropValue<TypeTag, Properties::SimulationName>();
        if (getPropValue<TypeTag, Properties::PressureDependentPermeability>())
            fileName << "PressureDependent";
        fileName << ".csv";
        csvFileName_ = fileName.str();

        // const char *cstr = csvFileName_.c_str();
        // if( remove( cstr ) != 0 )
        //     perror( "Error deleting file" );
        // calculateWellVolume_();

        shutin_ = 60 * 60 * 24 * 365 * 50;
		openinj_ = 0; //60 * 60 * 24 * 1 * 1; //wait one day
		max_time_step_size_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxTimeStepSize);
		boundaryFluxOut_.open (name() + "_flux.txt", std::ios::trunc);
		boundaryFluxOut_ << "Time TimeStep Total Boundary flux.\n";
		boundaryFluxOut_.close();
		total_boundary_flux_ = 0.0;
        total_boundary_flux2_ = 0.0;
		time_since_last_output_ = 1e9; // we want to output first timestep


    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters() {
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
        EWOMS_REGISTER_PARAM(TypeTag, bool, PressureDependentPermeability,
                             "If true, the permeability depends on pressure");
        // EWOMS_REGISTER_PARAM(TypeTag, bool, IrreversibleDeformation,
        //                      "If true, irreversible shear deformation is added");
		// EWOMS_REGISTER_PARAM(TypeTag, bool, EnableDeformationBand,
        //                      "If true, effect of defomation band is added");
        // EWOMS_REGISTER_PARAM(TypeTag, bool, EnableFaultBoundary,
        //                     "If true, effect of fault is added to the boundary");
        //EWOMS_REGISTER_PARAM(TypeTag, Scalar, DilationAngle, "Dilation angle of shear");
        //EWOMS_REGISTER_PARAM(TypeTag, std::string, CsvFilePath, "File path for mass leak");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InjectionRate,
                             "The injection rate [kg / s]");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, ReservoirLayerTop,
                             "Position of upper reservoir layer [m]");
        // EWOMS_REGISTER_PARAM(TypeTag, Scalar, FracturePermeability,
        //                      "Permeability in fracture [m2]");
		// EWOMS_REGISTER_PARAM(TypeTag, Scalar, EffectiveDeformationBandPermeability,
        //                      "Effective Deformationband Permeability Kf_aKm [1/m]");
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
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxDepth,
                             "The maximum depth [m] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
       // EWOMS_REGISTER_PARAM(TypeTag, int, DomainDim, "The domain dimension.");
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const {
        std::ostringstream oss;
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName) << "_"
            << Model::name();
        if (getPropValue<TypeTag, Properties::EnableEnergy>()) oss << "_ni";
        oss << "_rate_" << EWOMS_GET_PARAM(TypeTag, Scalar, InjectionRate);
        // oss << "_enabledeformationband_"
        //     << EWOMS_GET_PARAM(TypeTag, bool, EnableDeformationBand);
		//  oss << "_fracturePerm_"
        //     << EWOMS_GET_PARAM(TypeTag, Scalar, FracturePermeability);
		// oss << "_Kf_aKm_"
        //     << EWOMS_GET_PARAM(TypeTag, Scalar, EffectiveDeformationBandPermeability);
     //   oss << "_dim_" << EWOMS_GET_PARAM(TypeTag, int, DomainDim);
        return oss.str();
    }
    const std::string csvFileName() const {
        return csvFileName_;
    }
    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep() {
		boundaryFluxOut_.open (name() + "_flux.txt", std::ios::app);
		boundaryFluxOut_ << this->simulator().time() << " " << this->simulator().timeStepSize() << " " << -total_boundary_flux_ << std::endl;
		total_boundary_flux_ = 0.0;
        total_boundary_flux2_ = 0.0;
		boundaryFluxOut_.close();
    }

    bool shouldWriteOutput() {
		time_since_last_output_ += this->simulator().timeStepSize();
		if (time_since_last_output_ > (60*60*24*365)) //every 365th days
		{
			time_since_last_output_ = 0;
			return true;
		}
		return false;
	}

    static void appendLineToFile(std::string filepath, std::string line)
    {
        std::ofstream file;
        //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
        //file.exceptions(file.exceptions() | std::ios::failbit);
        file.open(filepath, std::ios::out | std::ios::app);
        if (file.fail())
            throw std::ios_base::failure(std::strerror(errno));

        //make sure write fails with exception if something is wrong
        file.exceptions(file.exceptions() | std::ios::failbit | std::ifstream::badbit);

        file << line;
    }

    Scalar nextTimeStepSize() const
    {
        Scalar dt = ParentType::nextTimeStepSize();
	   if (asImp_().simulator().time() < openinj_){ 
            return std::min(dt, openinj_ - asImp_().simulator().time());}
       else if (asImp_().simulator().time() < shutin_){
            return std::min(dt, shutin_ - asImp_().simulator().time());}
       else
			return std::min(dt, max_time_step_size_);
    }
    template <class Context>
    Scalar extrusionFactor(const Context& context OPM_UNUSED,
                           unsigned spaceIdx OPM_UNUSED,
                           unsigned timeIdx OPM_UNUSED) const {
        //return std::pow(asImp_().initialAperture(), dimWorld - dim);
		if (dim == dimWorld)
			return 1;
			
		return aperture0_;

    }

    Scalar extrusionFactor() const {
        throw;
    }
    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx,
                       unsigned timeIdx) const {
        return temperature_;
    }

    /*!
     * \brief Returns the intrinsic permeability of an intersection.
     *
     * This method is specific to the finite volume discretizations. If left unspecified,
     * it calls the intrinsicPermeability() method for the intersection's interior and
     * exterior finite volumes and averages them harmonically. Note that if this function
     * is defined, the intrinsicPermeability() method does not need to be defined by the
     * problem (if a finite-volume discretization is used).
     */
    template <class Context>
    void intersectionIntrinsicPermeability(DimMatrix& result,
                                           const Context& context,
                                           unsigned intersectionIdx,
                                           unsigned timeIdx,
                                           short upstreamDofIdx) const 
    {
        const auto& scvf = context.stencil(timeIdx).interiorFace(intersectionIdx);

        const DimMatrix& K1 = asImp_().intrinsicPermeability(context, scvf.interiorIndex(), timeIdx, scvf.interiorIndex());
        const DimMatrix& K2 = asImp_().intrinsicPermeability(context, scvf.exteriorIndex(), timeIdx, scvf.exteriorIndex());
        
		//const DimVector& normal = scvf.normal();
		//std::cout << normal << std::endl;
		//Scalar K1_tmp = 0.0;
		//Scalar K2_tmp = 0.0;

		//for (unsigned i = 0; i < dimWorld; ++i) {
		//	K1_tmp = std::max(K1_tmp, K1[i][i]*normal[i]);
		//	K2_tmp = std::max(K2_tmp, K2[i][i]*normal[i]);
		//}


        // entry-wise harmonic mean. this is almost certainly wrong if
        // you have off-main diagonal entries in your permeabilities!
        for (unsigned i = 0; i < dimWorld; ++i)
            for (unsigned j = 0; j < dimWorld; ++j) {
				//results[i][j] = 
                result[i][j] = Opm::harmonicMean(K1[i][j], K2[i][j]);
			}
    }
    template <class Context>
    void intersectionIntrinsicPermeability(DimMatrix& result,
                                           const Context& context,
                                           unsigned intersectionIdx,
                                           unsigned timeIdx) const 
    {
        ParentType::intersectionIntrinsicPermeability(result, context, intersectionIdx, timeIdx);
    }
    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                          unsigned timeIdx, short upstreamDofIdx) const {
        std::runtime_error(
            "not implemented: intrinsicPermeability(const Context& context, unsigned spaceIdx, usigned timeIdx, short upstreamDofIdx)"
            );
    }

    // template <class Context>
    // double deformationBandTransmissibilityMultiplicator(const Context& context, unsigned spaceIdx,
	// 													unsigned timeIdx) const {

	// if (!EWOMS_GET_PARAM(TypeTag, bool, EnableDeformationBand))
	// 	return 1.0;
		
	// //const auto& pos = context.template pos<0>(spaceIdx, timeIdx);
    // const auto& pos = context.pos(spaceIdx, timeIdx);
    // double t = throw_(pos[1]); //pass the y dim;
    // if (t <= 0)
    //     return 1.0;

    // double dz_width = damage_zone_width_full_(t);
    // auto perms_eff = damage_zone_permeability_(t);
	
	// const auto& stencil = context.stencil(/*timeIdx=*/0);
    // //const auto& face = stencil.template interiorFace<0>(0);
    // const auto& face = stencil.interiorFace(0);

	// const auto& faceCenter = face.integrationPos();
	// //const auto& interiorPos = stencil.template subControlVolume<0>(0).globalPos();
	// const auto& interiorPos = stencil.subControlVolume(0).globalPos();

    // // Calculate vector form cells to face
    // auto distVec0 = faceCenter;
    // for (unsigned dimIdx = 0; dimIdx < distVec0.size(); ++dimIdx) {
    //     distVec0[dimIdx] -= interiorPos[dimIdx];
    // }
	// Scalar dist = std::sqrt ( distVec0 * distVec0);
	// Scalar X = dz_width / dist;
    // //std::cout << "DEF " << t << " " << X << " " << perms_eff.first << " " << perms_eff.second <<std::endl;

	// // return multiplier
	// // T_DB * T_CELL_WO_DB / (T_CELL_WO_DB + T_DB) / T_CELL
	// // = perm_eff / (perm_eff ( 1 - X) + X )
    // return perms_eff.first / ( perms_eff.first * ( 1.0 - X ) + X);
	// }

    // template <class Context>
    // double boundaryTransMultiplier( const Context& context, unsigned scvfIdx,
	// 					        	unsigned timeIdx) const {

    // if (!EWOMS_GET_PARAM(TypeTag, bool, EnableFaultBoundary))
	// 	return 1.0;

	// const auto& stencil = context.stencil(/*timeIdx=*/0);
    // const auto& face = stencil.boundaryFace(scvfIdx);
    // unsigned spaceIdx = face.interiorIndex();
	// const auto& faceCenter = face.integrationPos();
    // const auto& interiorPos = stencil.subControlVolume(spaceIdx).globalPos();

    // // Calculate vector form cells to face
    // auto distVec0 = faceCenter;
    // for (unsigned dimIdx = 0; dimIdx < distVec0.size(); ++dimIdx) {
    //     distVec0[dimIdx] -= interiorPos[dimIdx];
    // }
	// Scalar dist = std::sqrt ( distVec0 * distVec0);
	// Scalar X = aperture0_ / (dist + aperture0_);
    // Scalar perm_eff = fractureK_[0][0]; // / kxx_[context.globalSpaceIndex(spaceIdx, timeIdx)][0][0]; 
	// // return multiplier
	// // T_FRAC * T_CELL / (T_CELL + T_FRAC) / T_CELL
	// // = perm_eff / (perm_eff ( 1 - X) + X )
    // Scalar mult = perm_eff / ( perm_eff * ( 1.0 - X ) + X);
    // //std::cout << mult << " " << deformationBandTransmissibilityMultiplicator(context, spaceIdx, timeIdx) << std::endl;
    // return mult * deformationBandTransmissibilityMultiplicator(context, spaceIdx, timeIdx);
    // }

	

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const {
		const auto& pos = context.pos(spaceIdx, timeIdx);
       // if (dim==dimWorld) {
			//return reservoirK_;

		// std::cout << "space idx: " <<context.globalSpaceIndex(spaceIdx, timeIdx)<< std::endl;
         //std::cout <<"resperm: " <<kxx_[context.globalSpaceIndex(spaceIdx, timeIdx)]<< " " << pos << std::endl;
        // std::cout << "finished"<<std::endl;
			return kxx_[context.globalSpaceIndex(spaceIdx, timeIdx)];
	//	}
       
	    //std::cout << "fracp: " << fractureK_ << " " << pos << std::endl;
       // return fractureK_;
    }

    // template <class Context>
    // const DimMatrix fracPermeability(const Context& context, unsigned spaceIdx,
    //                                       unsigned timeIdx, short upstreamDofIdx) const {
    //     Evaluation ap = asImp_().aperture(context, spaceIdx, timeIdx, upstreamDofIdx);
    //     Evaluation kr = Opm::pow(ap, 2) / 12.0;
    //     DimMatrix ki = this->toDimMatrix_(kr);
    //     if (dim == dimWorld) {
    //         ki[0][0] *= fracVolumeDensity();
    //         for (int i=1; i<dimWorld; i++)
    //             ki[i][i] *= fracVolumeDensity();//1e-20;
    //     }
    //     return ki;
    // }

    // template <class Context>
    // Evaluation aperture(const Context& context, unsigned spaceIdx,
    //                     unsigned timeIdx, short upstreamDofIdx) const {
    //     return asImp_().initialAperture() + apertureChange(context, spaceIdx, timeIdx, upstreamDofIdx);
    // }

    // template <class Context>
    // Evaluation apertureChange(const Context& context, unsigned spaceIdx,
    //                           unsigned timeIdx, short upstreamDofIdx) const {
    //     const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
    //     Evaluation p0;
    //     Evaluation p1;

    //     if (EWOMS_GET_PARAM(TypeTag, bool, PressureDependentPermeability)) {
    //         unsigned focusDofIdx = context.focusDofIndex();
    //         // we only carry the derivatives along if the upstream DOF is the one which
    //         // we currently focus on
    //         const auto& up = context.intensiveQuantities(upstreamDofIdx, timeIdx);
    //         const auto& down = context.intensiveQuantities(spaceIdx, timeIdx);
    //         if (upstreamDofIdx == static_cast<int>(focusDofIdx))
    //             p0 = up.fluidState().pressure(0);
    //         else
    //             p0 = Toolbox::value(up.fluidState().pressure(0));
    //         if (spaceIdx == static_cast<int>(focusDofIdx))
    //             p1 = down.fluidState().pressure(0);
    //         else
    //             p1 = Toolbox::value(down.fluidState().pressure(0));
    //         p0 = (p0 + p1) / 2.0;
    //     }
    //     else{
    //         Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
    //         initialFluidState_(fs, context, spaceIdx, timeIdx);
    //         p0 = 1e8; // fs.pressure(0);//
    //     }
    //     const auto& element = context.element();
    //     if (dim == 0) {
    //         return asImp_().initialAperture() * 10.0;
    //     }
    //     const auto sigma = calculateSigma(element);
    //     if (EWOMS_GET_PARAM(TypeTag, bool, IrreversibleDeformation))
    //         return reversibleApertureChange(sigma, p0) + irreversibleApertureChange(element, sigma, p0);
    //     else
    //         return  reversibleApertureChange(sigma, p0);
    // }

    // template <class Element>
    // std::vector<Scalar> calculateSigma(Element element) const {
    //     Scalar sigma_xx = 25.1 * 1e6 * maxDepth_ / 1000.0;
    //     Scalar sigma_yy = 15.8 * 1e6 * maxDepth_ / 1000.0;
    //     Scalar sigma_zz = 22.2 * 1e6 * maxDepth_ / 1000.0;

    //     Scalar theta;

    //     if (dim == dimWorld || dimWorld==3) {
    //         theta = 3.1415 / 8.0;//1.34591;
    //     }
    //     else if (dim == 1) {
    //         GlobalPosition faceNormal;
    //         const auto& face = element.template subEntity<1>(0);
    //         faceNormal = face.geometry().center() - element.geometry().center();
    //         Scalar weight = 0;
    //         for (auto& val : faceNormal) weight += val * val;
    //         faceNormal *= 1.0 / std::sqrt(weight);

    //         Dune::FieldVector<Scalar, dimWorld> fracNormal;
    //         fracNormal[0] = faceNormal[1];
    //         fracNormal[1] = -faceNormal[0];
    //         theta = std::atan(fracNormal[1] / fracNormal[0]);
    //     }
    //     else {
    //         throw std::runtime_error("Sigma only defined for dimension >=1");
    //     }
    //     Scalar sigma_n = std::abs((sigma_yy + sigma_zz) / 2.0 +
    //                               (sigma_yy - sigma_zz) / 2.0 * std::cos(2 * theta));
    //     Scalar sigma_t = std::abs(-(sigma_yy - sigma_zz) / 2.0 * std::sin(2 * theta));
    //     std::vector<Scalar> sigma{sigma_n, sigma_t};
    //     return sigma;
    // }

    // template <class Element>
    // Evaluation irreversibleApertureChange(const Element& element,
    //                                       const std::vector<Scalar>& sigma,
    //                                       Evaluation p) const {
    //     const auto& idxSet = this->simulator().gridView().indexSet();
    //     auto sigma_en = sigma[0];
    //     auto sigma_t = sigma[1];
    //     Scalar nu = 0.2;
    //     Scalar E = 10 * 1e9;
    //     Scalar mu = 0.5;
    //     Scalar amax = (asImp_().initialAperture() * 1.0);
    //     Scalar C = Opm::tan(dilationAngle_) * (2.0 * (1 - Opm::pow(nu, 2)) * 100.0 / E);
    //     // Evaluation sigma_c = sigma_t - mu * sigma_n;
    //     // Evaluation sigma_en = sigma_n + p;
    //     Evaluation pref = (2 * sigma_en - 2 * sigma_t / mu + amax / (C * mu)) / 2.0;
    //     Scalar k = Opm::min(1e-6, C * 3.1415 * mu / amax);
    //     Evaluation dilation = amax * (0.5 + 1.0 / 3.1415 * Opm::atan(k * (p - pref)));
    //     return dilation;//Opm::max(dilation, irreversibleDilation_[idxSet.index(element)]);
    // }

    // Evaluation reversibleApertureChange(const std::vector<Scalar>& sigma,
    //                                     Evaluation p0) const {
    //     Evaluation sigma_n = sigma[0];
    //     sigma_n -= p0;
    //     Scalar Kn = 100 * 1e6;
    //     Scalar Em = 0.9 * asImp_().initialAperture();
    //     // sigma_n = Opm::max(sigma_n, 0.0);
    //     // return - Opm::abs(sigma_n * Em / (Em * Kn + sigma_n));
    //     Scalar k = 2e-7;
    //     return -Em / (1 + Opm::exp(-2.0 * k * sigma_n));
    // }

    Scalar initialAperture() const {
        return aperture0_;
    }

    Scalar fracVolumeDensity() const {
        return asImp_().initialAperture() * asImp_().frac_density();
    }

    DimMatrix reservoirK() const {
        return reservoirK_;
    }

    // Scalar frac_density() const {
    //     return frac_density_;
    // }

    Scalar reservoirPorosity() const {
        return reservoirPorosity_;
    }
    // Scalar fracPorosity() const {
    //     return fracturePorosity_;
    // }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const {
        const auto& pos = context.pos(spaceIdx, timeIdx);
    //    if (dim==dimWorld) {
			//std::cout << "resphi: " << reservoirPorosity_[context.globalSpaceIndex(spaceIdx, timeIdx)] << " " << pos << std::endl;
            return reservoirPorosity_[context.globalSpaceIndex(spaceIdx, timeIdx)];
	//	}
       // else {
        //    Scalar phi = fracturePorosity_;
			//std::cout << "fracphi: " << phi << " " << pos << std::endl;
            //if (dim == dimWorld)
            //    phi *= fracVolumeDensity();
       //     return phi;
        //}
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context, unsigned spaceIdx,
                                               unsigned timeIdx) const {
        const auto& pos = context.pos(spaceIdx, timeIdx); 
       // if (dim == dimWorld)// && dim==dimWorld)
            return asImp_().reservoirMaterialParams();
      //  else            
       //     return asImp_().fractureMaterialParams();
    }

        template <class Context>
    const MaterialLawParams& materialLawParamsBoundary(const Context& context, unsigned spaceIdx,
                                                       unsigned timeIdx) const {      
            return asImp_().fractureMaterialParams();
    }

    const MaterialLawParams& reservoirMaterialParams() const{
        return reservoirMaterialParams_;
    }

    // const MaterialLawParams& fractureMaterialParams() const{
    //     return fractureMaterialParams_;
    // }

    /*!
     * \brief Return the parameters for the heat storage law of the rock
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    const SolidEnergyLawParams& solidEnergyLawParams(const Context& context OPM_UNUSED,
                                                     unsigned spaceIdx OPM_UNUSED,
                                                     unsigned timeIdx OPM_UNUSED) const {
        return solidEnergyLawParams_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams& thermalConductionLawParams(const Context& context,
                                                                 unsigned spaceIdx,
                                                                 unsigned timeIdx) const {
        return coarseThermalCondParams_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context, unsigned spaceIdx,
                  unsigned timeIdx) const {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (asImp_().onDirichletBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx,false);
            fs.checkDefined();
            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
			total_boundary_flux_ += Opm::getValue(values[conti0EqIdx]);
        }
        else
            values.setNoFlow();
    }

    // \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        values.assignNaive(fs);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate, const Context& context,
                unsigned spaceIdx, unsigned timeIdx) const {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (asImp_().inWell_(pos) && asImp_().simulator().time() < (shutin_ - eps_) && asImp_().simulator().time() > (openinj_ - eps_)) {
			//std::cout << pos << " " << injectionRate_ << " " << asImp_().getWellVolume() << std::endl;
            rate[contiCO2EqIdx] = injectionRate_ / asImp_().getWellVolume();
        }
        else if (false && dim == 2) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.checkDefined();
            BoundaryContext boundaryCtx(context);
            BoundaryRateVector values;
            //FluxExtensiveQuantities::calculateBoundaryGradients_(boundaryCtx,
            //                                                spaceIdx,
            //                                                 timeIdx,
            //                                                 fs);
            // impose an freeflow boundary condition as a source term 

            values.setFreeFlow(boundaryCtx, spaceIdx, timeIdx, fs);
            //rate = values;
			total_boundary_flux2_ += Opm::getValue(values[conti0EqIdx]);
        }
        else
            rate = Scalar(0.0);
    }

    Scalar getWellVolume() const{
        return wellVolume_;
    }
    //! \}

    template <class Vec>
    GlobalPosition dofPos(const Vec& pos) const{
        return pos;
    }

   protected:
    
    void calculateWellVolume_(){
        wellVolume_ = 0.0;
        for (auto& ele : elements(this->simulator().gridView())) {
            if (asImp_().inWell_(ele.geometry().center()))
                wellVolume_ += ele.geometry().volume();
        }
        if (dimWorld==2)
            wellVolume_ *= 1000.0;
    }

    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs, const Context& context, unsigned spaceIdx,
                            unsigned timeIdx, bool boundary = false) const {
        const GlobalPosition pos = asImp_().dofPos(context.pos(spaceIdx, timeIdx));

        //////
        // set temperature
        //////
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));

        //////
        // set saturations
        //////
        Scalar Sw = 1.0;
        fs.setSaturation(FluidSystem::liquidPhaseIdx, Sw);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 1 - Sw);

        //////
        // set pressures
        //////
        Scalar densityL = FluidSystem::Brine::liquidDensity(temperature_, Scalar(1e5));
        Scalar depth = maxDepth_ - pos[dimWorld - 1];
        Scalar pl = 1e5 - densityL * this->gravity()[dimWorld - 1] * depth;
        if (false && dim == 2){
            std::cout << "pos = " << context.pos(spaceIdx, timeIdx);
            std::cout << ", dim dim world = " << dim << " " << dimWorld << std::endl;
			std::cout << "grav = " << this->gravity()[dimWorld - 1] << std::endl;
			std::cout << pl << std::endl;
        }
        Scalar pC[numPhases];
        if (boundary) {
            const auto& matParams = this->fractureMaterialParams();
            MaterialLaw::capillaryPressures(pC, matParams, fs);
        } else {
            const auto& matParams = this->reservoirMaterialParams();
            MaterialLaw::capillaryPressures(pC, matParams, fs);
        }

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

    bool inReservoir_(const GlobalPosition& pos) const {
        return (pos[dimWorld - 1] <= 20 || pos[dimWorld - 1] >= 20 + 100);
    }

    bool onLeftBoundary_(const GlobalPosition& pos) const { return pos[1] < 6.69e6; }
    bool onRightBoundary_(const GlobalPosition& pos) const { return pos[0] > 2e3 - eps_; }
    bool onTopBoundary_(const GlobalPosition& pos) const { return (pos[dimWorld - 1] > -840 - eps_); }
    bool onLowerBoundary_(const GlobalPosition& pos) const { return pos[dimWorld - 1] < eps_; }
    bool onDomainBoundary_(const GlobalPosition& pos) const {
        return asImp_().onLeftBoundary_(pos) || asImp_().onRightBoundary_(pos) || asImp_().onLowerBoundary_(pos) ||
               asImp_().onTopBoundary_(pos);
    }
    bool onDirichletBoundary_(const GlobalPosition& pos) const {
        return onDirichletBoundaryDim3_(pos) || onDirichletBoundaryDim2_(pos);
    }
    bool onDirichletBoundaryDim2_(const GlobalPosition& pos) const {
        return (dim == 2 && (pos[dimWorld - 1] > -840 - eps_) && (pos[1] > 6.75677e6 - 1000));
    }
    // bool onDirichletBoundaryDim3_(const GlobalPosition& pos) const {
    //     return (dim == 3 && (pos[0] < 560237) && EWOMS_GET_PARAM(TypeTag, bool, EnableFaultBoundary) );
    // }
    bool onInlet_(const GlobalPosition& pos) const { return false; }

    bool inWell_(const GlobalPosition& pos) const {
        return (pos[0] > 559600) && (pos[1] > 6.7237e6) && (pos[0] < 560001) && (pos[1] < 6.72415e6) && (pos[2] < -1488) && (pos[2] > -1545);
    }

    bool inHighTemperatureRegion_(const GlobalPosition& pos) const { return false; }

    void computeThermalCondParams_(ThermalConductionLawParams& params, Scalar poro) {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet =
            std::pow(lambdaGranite, (1 - poro)) * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1 - poro));

        params.setFullySaturatedLambda(gasPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(liquidPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    private:
    void loadPermeability_(){
		
		if (dim < dimWorld)
		return;
        std::ifstream file;
        //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
        //file.exceptions(file.exceptions() | std::ios::failbit);


        file.open("/home/AD.NORCERESEARCH.NO/tmyk/OPM/opm-data/smeaheia/smeaheia_perm.txt", std::ios::in);
        if (file.fail())
            throw std::ios_base::failure(std::strerror(errno));

        std::string line;
        int i = 0;
        kxx_.resize(this->simulator().gridView().size(0), reservoirK_);
        while ( getline (file, line) )
        {
            kxx_[i] = this->toDimMatrix_(std::stod(line));
            kxx_[i] *= 1e-3; // from d to md
			kxx_[i][dimWorld - 1][dimWorld - 1] *= 0.1;
			//kxx_[i][0][0] *= 0.1;
			//std::cout << i << " " << kxx_[i] << std::endl;
            i++;
        }
        //std::cout << "size " << kxx_.size() << std::endl;
        file.close();
    }
    void loadPorosity_(){
		if (dim < dimWorld) {
		    reservoirPorosity_.resize(1);
			reservoirPorosity_[0] = 1.0;
			return;
		}
        std::ifstream file;
        //can't enable exception now because of gcc bug that raises ios_base::failure with useless message
        //file.exceptions(file.exceptions() | std::ios::failbit);
        file.open("/home/AD.NORCERESEARCH.NO/tmyk/OPM/opm-data/smeaheia/smeaheia_poro.txt", std::ios::in);
        if (file.fail())
            throw std::ios_base::failure(std::strerror(errno));

        std::string line;
        int i = 0;
        reservoirPorosity_.resize(this->simulator().gridView().size(0));
        while ( getline (file, line) )
        {
            reservoirPorosity_[i] = std::stod(line);
            i++;
        }
        file.close();
    }

// Calculate the effective permeability over the whole damage zone
// std::pair<double,double> damage_zone_permeability_(double t) const {
//         double W5 = damage_zone_width_(t);
//         double Kf_aKm = EWOMS_GET_PARAM(TypeTag, Scalar, EffectiveDeformationBandPermeability);
//         double band_length=2;
//         double sigma=M_PI/12;
//         std::pair<double,double> density_params = band_density_params_(W5);
//         double A = density_params.first;
//         double L = density_params.second;
//         double X = std::exp(-A / L);
//         std::pair<double,double> Kxy = interval_effective_perm_(0, X, t, W5, Kf_aKm, band_length, sigma);
//         return Kxy;
// }

double throw_(double y) const {
  //std::cout << y << std::endl;
  double maxPoint = 6739819.56;
  double minPoint = 6667544.54;
  //std::cout << 1000 * (1.0 - (abs(maxPoint - y) / (maxPoint - minPoint)) ) << std::endl;
  return std::max(0.0, 1000 * (1.0 - (abs(maxPoint - y) / (maxPoint - minPoint)) ));
}



// double damage_zone_width_(double t) const  {
//     //Estimate damage zone width from throw. Average value from Schueller (2013), Fig 6
//     return 1.74 * std::pow(t, 0.43);
// }

// double damage_zone_width_full_(double t)const  {
//     double W5 = damage_zone_width_(t);
//     std::pair<double,double> density_params = band_density_params_(W5);
//     double A = density_params.first;
//     double L = density_params.second;
//     return std::exp(-A / L);
// }

// std::pair<double,double> interval_effective_perm_(double a, double b, double t, double W5, double Kf_aKm, double band_length, double sigma) const {
//     //Calculate the effective permeability at an interval from a fault.  
//     //Params:
//     //a [float]: The starting point of the interval
//     //b [float]: The end point of the interval
//     //throw [float]: The throw of the fault. Either throw or W5 must be set, but not both.
//     //W5 [float] : The damage zone width of the fault.Either throw or W5 must be set, but not both.
//     //Kf_aKm [float]: The permeability contrast between the bands and the rock matrix. The
//     //                permeability contrast is defined as perm_band / (aperture * perm_matrix)
//     //band_length [float]: The length of individual bands.
//     //sigma [float]: The standard deviation of the rotation of the bands
//     //Returns:
//     //Kex_km [float]: The effective permeability in x-direction divided by matrix permeability
//     //Key_km [float]: The effective permeability in y-direction divided by matrix permeability
    
//     // discretize segment
//     double segment_length = b - a;
//     double dx = segment_length / 100;
//     // use midtpoint (mostly to avoid infinite value at 0)
//     double x = dx/2;
//     double Kx = 0.0;
//     double Ky = 0.0;
//     for (int i = 0; i < 100; i++) {
//        auto Kp = point_wise_perm_(x, W5, Kf_aKm, band_length, sigma);
//        Kx += dx / Kp.first;
//        Ky += Kp.second;
//        x+=dx;
//     }
//     // Harmonic average in x-direction
//     Kx = segment_length / Kx;

//     // Arithmetic average in y-direction
//     Ky = Ky / 100;

//     return std::make_pair(Kx, Ky);
// }

// std::pair<double,double> band_density_params_(double W5) const {
//     // Estimate band density:
//     double D5 = 13.33;  //Averaged band density across the whole damage zone
//     // Assuming D5 is independent of W5:
//     double L = 5 - D5;
//     double A = D5 - L * (std::log(W5) - 1);
//     //Band density as function of distance
//     return std::make_pair(A, L);
// }


// std::pair<double,double> point_wise_perm_(double distance, double W5, double Kf_aKm, double band_length, double sigma) const {
//     //Calculate the point-wise effective permeability at a given distance from a fault.
//     ///
//     ///Params:
//     ///distance [float or array]: The distance from the fault
//     ///W5 [float] : The damage zone width of the fault.
//     //Kf_aKm [float]: The permeability contrast between the bands and the rock matrix. The
//     //                permeability contrast is defined as perm_band / (aperture * perm_matrix)
//     //band_length [float]: The length of individual bands.
//     //sigma [float]: The standard deviation of the rotation of the bands

//     //Returns:
//     //Kx_Km : The effective permeability in x-direction divided by matrix permeability (Kex / Km)
//     //Ky_Km : The effective permeability in y-direction divided by matrix permeability (Key / Km)
    
//     std::pair<double,double> density_params = band_density_params_(W5);
//     double A = density_params.first;
//     double L = density_params.second;
//     double rho = A + L * std::log(distance);
//     return effective_permeability_(rho, Kf_aKm, band_length, sigma);
// }


std::pair<double,double> effective_permeability_(double rho, double Kf_aKm, double band_length, double sigma) const  {

    // Avoid dividing by 0 (rho=0 just gives us rock matrix permeability anyway)
    if (rho < 1e-8)
        return std::make_pair(1.0,1.0);

    // The expected value of the rotation
    double theta = sigma * std::sqrt(2) / std::sqrt(M_PI);

    // Chain length:
    std::pair<double,double> cl = chain_length_(rho, band_length, sigma);
    double Lx = cl.first;
    double Ly = cl.second;

    // permeability x-direction
    double Yy = Ly + 1 / (2 * rho);
    double Kfy = 1 / (1 + std::pow(Kf_aKm,-1) * rho);
    double Kx = (Kfy * Ly + 1.0 / (2 * rho)) / Yy;

    // permeability y-direction
    double num_cross = rho * band_length * std::sin(theta);
    double Yx = Lx + 1 / rho;
    double Kfx = 1 / (1 + std::pow(Kf_aKm,-1) * num_cross * std::sin(theta) / (band_length * std::cos(theta)));
    double Ky = (Kfx * Lx + 1 / rho) / Yx;
    return std::make_pair(Kx, Ky);
}



std::pair<double,double> chain_length_(double rho, double band_length, double sigma) const {
    double theta = sigma * std::sqrt(2) / std::sqrt(M_PI);
    // probability that a segment contains a crossing:
    double p_cross = 1 - std::exp(-band_length * rho * std::sin(theta));

    // Expected number of crossings in a chain:
    double Ecross = p_cross / (1 - p_cross);

    // Expected lenght of chain in x- and y-direction
    double Lx = std::sin(theta) * (band_length + (0.5 * band_length * Ecross));
    double Ly = std::cos(theta) * (band_length + (0.5 * band_length * Ecross));
    return std::make_pair(Lx, Ly);
}

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    std::vector<DimMatrix> kxx_;
    DimMatrix reservoirK_;
	DimMatrix fractureK_;
    DimMatrix fineK_;
    Scalar aperture0_;
    Scalar dilationAngle_;
    std::vector<Scalar> irreversibleDilation_;

    std::string csvFilePath_;
    std::string csvFileName_;

    Scalar reservoirLayerBottom_;
    Scalar reservoirLayerTop_;

    Scalar fracturePorosity_;
    std::vector<Scalar> reservoirPorosity_;
    Scalar frac_density_;

    MaterialLawParams fractureMaterialParams_;
    MaterialLawParams reservoirMaterialParams_;

    ThermalConductionLawParams fineThermalCondParams_;
    ThermalConductionLawParams coarseThermalCondParams_;
    SolidEnergyLawParams solidEnergyLawParams_;

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;

    Scalar injectionRate_;
    Scalar wellVolume_;

    unsigned nTemperature_;
    unsigned nPressure_;
    Scalar shutin_;
	Scalar openinj_;
	Scalar max_time_step_size_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
	Scalar time_since_last_output_;
	std::ofstream boundaryFluxOut_;
	mutable Scalar total_boundary_flux_;
    mutable Scalar total_boundary_flux2_;
    

};
}  // namespace Opm

#endif
