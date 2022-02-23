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
 * \copydoc Opm::Co2SmeaheiaInjectionProblem
 */
#ifndef CO2_SMEAHEIA_INJECTION_PROBLEM_HH
#define CO2_SMEAHEIA_INJECTION_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
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
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/cpgrid/dgfparser.hh>
#include <opm/grid/polyhedralgrid/dgfparser.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/models/common/transfluxmodule.hh>
#include <opm/models/io/unstructuredgridvanguard.hh>
#include <opm/simulators/linalg/parallelamgbackend.hh>
#include <sstream>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>

namespace Opm {
//! \cond SKIP_THIS
template <class TypeTag>
class Co2SmeaheiaInjectionProblem;

namespace Co2SmeaheiaInjection {
#include <opm/material/components/co2tables.inc>
}
//! \endcond
}  // namespace Opm

namespace Opm::Properties {

namespace TTag {
struct Co2SmeaheiaInjectionBaseProblem {};
}  // namespace TTag
template <class TypeTag, class MyTypeTag>
struct InjectionRate { 
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct ReservoirLayerTop {
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct DomainDim { 
    using type = UndefinedProperty; 
};
template <class TypeTag, class MyTypeTag>
struct WorldDim { 
    using type = UndefinedProperty; 
};

template <class TypeTag, class MyTypeTag>
struct GridDim { 
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


// Enable adaptivity
template <class TypeTag>
struct EnableGridAdaptation<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    static constexpr bool value = false;
};
// No async vtk with adaptivity
template <class TypeTag>
struct EnableAsyncVtkOutput<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    static constexpr bool value = true;
};

// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = Opm::Co2SmeaheiaInjectionProblem<TypeTag>;
};

// Set the problem property
template <class TypeTag>
struct FluxModule<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = TransFluxModule<TypeTag>;
};


// Set fluid configuration
template <class TypeTag>
struct FluidSystem<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
   private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = Opm::Co2SmeaheiaInjection::CO2Tables;

   public:
    using type = Opm::BrineCO2FluidSystem<Scalar, CO2Tables>;
    // using type = Opm::H2ON2FluidSystem<Scalar, /*useComplexRelations=*/false>;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem>
{
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
struct ReservoirLayerTop<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 120.0;
};
// Set the thermal conduction law
template <class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    // define the material law parameterized by absolute saturations
    using type = Opm::SomertonThermalConductionLaw<FluidSystem, Scalar>;
};

// set the energy storage law for the solid phase
template <class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    using type = Opm::ConstantSolidHeatCapLaw<GetPropType<TypeTag, Properties::Scalar>>; 
};

// Use the algebraic multi-grid linear solver for this problem
template <class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    using type = TTag::ParallelAmgLinearSolver; 
};

// Write the Newton convergence behavior to disk?
template <class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    static constexpr bool value = false; 
};

// By default, include the intrinsic permeability tensor to the VTK output files
template <class TypeTag>
struct VtkWriteIntrinsicPermeabilities<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    static constexpr bool value = true;
};

// Enable gravity
template <class TypeTag>
struct EnableGravity<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    static constexpr bool value = true; 
};

template <class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = TTag::EcfvDiscretization;
};

// use automatic differentiation for this simulator
template <class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = TTag::AutoDiffLocalLinearizer;
};

// set the defaults for the problem specific properties
template <class TypeTag>
struct FluidSystemPressureLow<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e5;
};
template <class TypeTag>
struct FluidSystemPressureHigh<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e8;
};
template <class TypeTag>
struct FluidSystemNumPressure<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    static constexpr int value = 100; 
};
template <class TypeTag>
struct FluidSystemTemperatureLow<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 280;
};
template <class TypeTag>
struct FluidSystemTemperatureHigh<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 400;
};
template <class TypeTag>
struct FluidSystemNumTemperature<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    static constexpr int value = 50; 
};

template <class TypeTag>
struct MaxDepth<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.0;
};
template <class TypeTag>
struct Temperature<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 363.15;
};

template <class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    static constexpr bool value = false;
};

template <class TypeTag>
struct SimulationName<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    static constexpr auto value = "Smeaheia"; 
};


// The default for the end time of the simulation
template <class TypeTag>
struct EndTime<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3 * 365 * 24 * 60 * 60;
};

// The default for the initial time step size of the simulation
template <class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.0 * 24 * 60 * 60;
};

// The default for the initial time step size of the simulation
template <class TypeTag>
struct MaxTimeStepSize<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 30 * 24 * 60 * 60;
};


template <class TypeTag>
struct WorldDim<TypeTag,  TTag::Co2SmeaheiaInjectionBaseProblem> {
    static constexpr auto value = 3; };
template <class TypeTag>
struct DomainDim<TypeTag,  TTag::Co2SmeaheiaInjectionBaseProblem> {
    static constexpr auto value = 3; 
};
template <class TypeTag>
struct GridDim<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    static constexpr auto value = 3; 
};
template <class TypeTag>
struct Grid<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> {
    private:
        enum{ worldDim = getPropValue<TypeTag, Properties::WorldDim>() };
    public:
    using type = Dune::PolyhedralGrid<3, worldDim>;
};

template <class TypeTag>
struct Vanguard<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    using type = Opm::UnstructuredGridVanguard<TypeTag>; };

template <class TypeTag>
struct GridFile<TypeTag, TTag::Co2SmeaheiaInjectionBaseProblem> { 
    static constexpr auto value = "./data/smeaheia_3.txt"; 
};

template <class TypeTag>
struct InjectionRate<TypeTag,  TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 54.0;
};
template <class TypeTag>
struct NewtonTolerance<TypeTag,  TTag::Co2SmeaheiaInjectionBaseProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};






}  // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief CO2 injection in Smeaheia
 */
template <class TypeTag>
class Co2SmeaheiaInjectionProblem 
    : public GetPropType<TypeTag, Properties::BaseProblem>{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;
    
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using BoundaryContext = GetPropType<TypeTag, Properties::BoundaryContext>;

    enum { dim = getPropValue<TypeTag, Properties::DomainDim>() };
    enum { dimWorld = getPropValue<TypeTag, Properties::WorldDim>() };

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

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Co2SmeaheiaInjectionProblem(Simulator& simulator)
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

        // intrinsic permeabilities
        reservoirK_ = this->toDimMatrix_(1e-13);
        reservoirK_[dimWorld - 1][dimWorld - 1] *= 0.1;
        coarseK_ = this->toDimMatrix_(1e-12);

        // porosities
        
        coarsePorosity_ = 0.3;
        loadPorosity_();//reservoirPorosity_ = 0.3;
        loadPermeability_();

        // residual saturations
        reservoirMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.2);
        reservoirMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);
        coarseMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.2);
        coarseMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the Brooks-Corey law
        reservoirMaterialParams_.setEntryPressure(1e4);
        coarseMaterialParams_.setEntryPressure(5e3);
        reservoirMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);

        reservoirMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        // parameters for the somerton law thermal conduction
        computeThermalCondParams_(coarseThermalCondParams_, reservoirPorosity_[0]);

        // assume constant heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(
            790.0       // specific heat capacity of granite [J / (kg K)]
            * 2700.0);  // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();
        injectionRate_ = EWOMS_GET_PARAM(TypeTag, Scalar, InjectionRate);


        calculateWellVolume_();
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
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InjectionRate,
                             "The injection rate [kg / s]");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, ReservoirLayerTop,
                             "Position of upper reservoir layer [m]");
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
        EWOMS_REGISTER_PARAM(TypeTag, int, DomainDim, "The domain dimension.");
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
        return oss.str();
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

        /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx,
                       unsigned timeIdx) const {
        return temperature_;
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
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const {
        return kxx_[context.globalSpaceIndex(spaceIdx, timeIdx)];
    }
    DimMatrix reservoirK() const {
        return reservoirK_;
    }
    Scalar reservoirPorosity() const {
        return reservoirPorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        return reservoirPorosity_[context.globalSpaceIndex(spaceIdx, timeIdx)];
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context, unsigned spaceIdx,
                                               unsigned timeIdx) const {
        return asImp_().reservoirMaterialParams();
    }

    const MaterialLawParams& reservoirMaterialParams() const{
        return reservoirMaterialParams_;
    }

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
    void boundary(BoundaryRateVector& values, const Context& context,unsigned spaceIdx, 
                unsigned timeIdx) const {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (asImp_().onDirichletBoundary_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.checkDefined();
            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else
            // no flow on top and bottom
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
     */
    template <class Context>
    void source(RateVector& rate, const Context& context,
                unsigned spaceIdx, unsigned timeIdx) const {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (asImp_().inWell_(pos) && asImp_().simulator().time() < (shutin_ - eps_) && asImp_().simulator().time() > (openinj_ - eps_)) {
			//std::cout << pos << " " << injectionRate_ << " " << asImp_().getWellVolume() << std::endl;
            rate[contiCO2EqIdx] = injectionRate_ / asImp_().getWellVolume();
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
                            unsigned timeIdx) const {
        //const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
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
        Scalar pC[numPhases];
        const auto& matParams = this->reservoirMaterialParams();
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
    bool onDirichletBoundaryDim3_(const GlobalPosition& pos) const {
        return false;
    }
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


        file.open("data/smeaheia_perm.txt", std::ios::in);
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
        file.open("data/smeaheia_poro.txt", std::ios::in);
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

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    std::vector<DimMatrix> kxx_;
    DimMatrix reservoirK_;
    DimMatrix coarseK_;
    Scalar reservoirLayerBottom_;
    Scalar reservoirLayerTop_;

    std::vector<Scalar> reservoirPorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams reservoirMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    ThermalConductionLawParams reservoirThermalCondParams_;
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
