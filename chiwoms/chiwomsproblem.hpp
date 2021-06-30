#ifndef CHIWOMS_PROBLEM_HPP
#define CHIWOMS_PROBLEM_HPP

#include <iostream>
#include <cassert>
#include <stdexcept>  // invalid_argument
#include <sstream>
#include <iostream>
#include <string>
#include <random>    // mt19937, normal_distribution
#include <limits>    // epsilon
#include <boost/format.hpp>  // boost::format

// note that this header must be included before anything that directly
// or indirectly may reference DUNE headers
#include "threephasefluidsystem.hh"
#include "ChiFlash.hpp"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/start.hh>

#include <opm/models/io/structuredgridvanguard.hh>

//#include <ewoms/models/ncp/ncpmodel.hh>
#include <opm/models/flash/flashmodel.hh>
#include <opm/models/immiscible/immisciblemodel.hh>

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>

#include  <opm/simulators/linalg/parallelamgbackend.hh>

#include <opm/common/Exceptions.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Unused.hpp>

//#include <dune/alugrid/grid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>


namespace Opm {
template <class TypeTag>
class ChiwomsProblem;
} // namespace Opm

namespace Opm::Properties {

// Create new type tags
namespace TTag {
struct ChiwomsProblem {};
} // end namespace TTag


template<class TypeTag>
struct Grid<TypeTag, TTag::ChiwomsProblem>
{ using type = Dune::YaspGrid<2>; };

//SET_TYPE_PROP(ChiwomsProblem, Grid, Dune::YaspGrid<2>);

// declare the CO2 finger problem specific property tags
template<class TypeTag, class MyTypeTag>
struct TopResvPres{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct Temperature{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct SimulationName{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct EpisodeLength{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct WaveLength{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct ConnateWater{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct ResidualOil{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct EntryPressure{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct PoreSizeDist{ using type = UndefinedProperty; };


// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ChiwomsProblem> { using type = Opm::ChiwomsProblem<TypeTag>; };

template<class TypeTag>
struct FlashSolver<TypeTag, TTag::ChiwomsProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    using type = Opm::ChiFlash<Scalar, Evaluation, FluidSystem >;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChiwomsProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using type = Opm::ThreePhaseThreeComponentFluidSystem<Scalar>;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::ChiwomsProblem>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum {
        oilPhaseIdx = FluidSystem::oilPhaseIdx,
        waterPhaseIdx = FluidSystem::waterPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::ThreePhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                               /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                               /*gasPhaseIdx=*/ FluidSystem::gasPhaseIdx>;

    // define the material law which is parameterized by effective
    // saturations

    using EffMaterialLaw = Opm::NullMaterial<Traits>;

public:
    // define the material law parameterized by absolute saturations
    using type = EffMaterialLaw;
};

// Write the Newton convergence behavior to disk?
template<class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = false; };

// Enable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = false; };

// set the defaults for the problem specific properties
template<class TypeTag>
struct TopResvPres<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 200 * 1.e5;
};

template<class TypeTag>
struct Temperature<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 273.15 + TEMPERATURE;
};


template<class TypeTag>
struct SimulationName<TypeTag, TTag::ChiwomsProblem>
{ static constexpr auto value = "chiwoms"; };


// The default for the end time of the simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value =  SIM_TIME * 24. * 60. * 60.;
};


// convergence control
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 30;
};

template<class TypeTag>
struct LinearSolverTolerance<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};

template<class TypeTag>
struct LinearSolverAbsTolerance<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.;
};
template<class TypeTag>
struct NewtonTolerance<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};
template<class TypeTag>
struct MaxTimeStepSize<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 60 * 60;
};
template<class TypeTag>
struct NewtonMaxIterations<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 10;
};
template<class TypeTag>
struct NewtonTargetIterations<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 6;
};

// output
template<class TypeTag>
struct VtkWriteFilterVelocities<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWritePotentialGradients<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteTotalMassFractions<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteTotalMoleFractions<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteFugacityCoeffs<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = true; };

// write restart for every hour
template<class TypeTag>
struct EpisodeLength<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 60. * 60.;
};

// mesh grid
template<class TypeTag>
struct Vanguard<TypeTag, TTag::ChiwomsProblem> { using type = Opm::StructuredGridVanguard<TypeTag>; };

template<class TypeTag>
struct DomainSizeX<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = (X_SIZE / 100.);
};
template<class TypeTag>
struct CellsX<TypeTag, TTag::ChiwomsProblem> { static constexpr int value = NX; };
template<class TypeTag>
struct DomainSizeY<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = (Y_SIZE / 100.);
};
template<class TypeTag>
struct CellsY<TypeTag, TTag::ChiwomsProblem> { static constexpr int value = NY; };
template<class TypeTag>
struct DomainSizeZ<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.;
};
template<class TypeTag>
struct CellsZ<TypeTag, TTag::ChiwomsProblem> { static constexpr int value = 1; };


// compositional, with diffusion
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = false; };

// injection rate parameter
template<class TypeTag>
struct WaveLength<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = WAVE_LENGTH;
};

// default hydrology: almost (but not quite) Snohvit-like Brooks-Corey
template<class TypeTag>
struct ConnateWater<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.15;
};
template<class TypeTag>
struct ResidualOil<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.10;
};
template<class TypeTag>
struct EntryPressure<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1.5e4 /* [Pa] */;
};
template<class TypeTag>
struct PoreSizeDist<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 2.0;
};

}// namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 */
template <class TypeTag>
class ChiwomsProblem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;


    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    enum { numPhases = FluidSystem::numPhases };
    enum {
        oilPhaseIdx = FluidSystem::oilPhaseIdx,
        waterPhaseIdx = FluidSystem::waterPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
    };
    enum { Comp1Idx = FluidSystem::Comp1Idx };//change to comp1
    enum { Comp0Idx = FluidSystem::Comp0Idx };//change to comp0
    enum { Comp2Idx = FluidSystem::Comp2Idx };//change to comp2
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + Comp1Idx };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;

    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;

    using Toolbox = Opm::MathToolbox<Evaluation>;
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using FullField = Dune::FieldMatrix<Scalar, NX, NY>;//    typedef Dune::FieldMatrix<Scalar, NX, NY> FullField;
    using FieldColumn = Dune::FieldVector<Scalar, NY>;//typedef Dune::FieldVector<Scalar, NY> FieldColumn;

    const unsigned XDIM = 0;
    const unsigned YDIM = 1;

    // we need to layout the initial state from top to down, so fill
    // a full matrix with that now and access it (possibly) randomly later
    // note that these matrices are logical to fill from the top, so index
    // 0 is the *top-most* element, and the index increase as we move *down*
    // the column. this is opposite of the orientation of the grid in eWoms
    FieldColumn init_pres;  // oleic phase pressure; ref. pres. w/o cap. pres

    // influx on the left boundary
    Scalar rate;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    ChiwomsProblem(Simulator& simulator)
        : ParentType(simulator)
    {
	    Scalar epi_len = EWOMS_GET_PARAM(TypeTag, Scalar, EpisodeLength);
	    simulator.setEpisodeLength(epi_len);
    }

    void initPetrophysics() {
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        K_ = this->toDimMatrix_(PERMEABILITY * 1.e-15);
        porosity_ = POROSITY;
    }

    void initHydrology() {
	    this->mat_.finalize();
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        // set seed so that any random number generated is reused next run
        rand_gen.seed(SEED);

        // initialize fixed parameters; temperature, permeability, porosity
        initPetrophysics();

        // initialize two-phase unsaturated zone functions
        initHydrology();

    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, TopResvPres,
                             "Pressure [Pa] at the top of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EpisodeLength,
                             "Time interval [s] for episode length");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, WaveLength,
                             "Instill critical wavelength for this fraction of domain height");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, ConnateWater,
                             "Saturation of water that cannot be removed once there");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, ResidualOil,
                             "Saturation of oil that cannot be removed once there");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EntryPressure,
                             "Extra pressure required for oil to enter water-filled pore");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, PoreSizeDist,
                             "Distribution parameter for pore size");
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
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName);
        return oss.str();
    }

    // This method must be overridden for the simulator to continue with
    // a new episode. We just start a new episode with the same length as
    // the old one.
    void endEpisode() {
	    Scalar epi_len = EWOMS_GET_PARAM(TypeTag, Scalar, EpisodeLength);
	    this->simulator().startNextEpisode(epi_len);
    }

    // only write output when episodes change, aka. report steps, and
    // include the initial timestep too
    bool shouldWriteOutput() {
        return true;
	    return this->simulator().episodeWillBeOver()
		    || (this->simulator().timeStepIndex() == -1);
    }

    // we don't care about doing restarts from every fifth timestep, it
    // will just slow us down
    bool shouldWriteRestartFile() {
	    return false;
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
        Scalar tol = this->model().newtonMethod().tolerance()*1e5;
        this->model().checkConservativeness(tol);

        // Calculate storage terms
        PrimaryVariables storageO, storageW;
        this->model().globalPhaseStorage(storageO, oilPhaseIdx);
        this->model().globalPhaseStorage(storageW, waterPhaseIdx);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
	        std::cout << "Storage: oleic = [" << storageO << "], "
	                  << "aqueous = [" << storageW << "]" << std::endl << std::flush;
        }
    }

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFs(fs, context, spaceIdx, timeIdx);
        values.assignNaive(fs);

    //std::cout << "primary variables for cell " << context.globalSpaceIndex(spaceIdx, timeIdx) << ": " << values << "\n";
    }

    // Constant temperature
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
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
    Scalar porosity(const Context& context  OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    {
        return porosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        return this->mat_;
    }

    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        //values.setNoFlow();
        //return;

        const Scalar eps = std::numeric_limits<double>::epsilon();
	    const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        // left side has a fixed inflow rate
	    if((pos[XDIM] < this->boundingBoxMin()[XDIM] + eps) )
        {
		    // assign rate to the CO2 component of the inflow
		    RateVector massRate(0.);
            massRate[contiCO2EqIdx] = -1e-3;// -1e-7;
		    values.setMassRate(massRate);
        } 
        else if((pos[XDIM] > this->boundingBoxMax()[XDIM] - eps))
        {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            initialFs(fs, context, spaceIdx, timeIdx);
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        } 
	    else {
            values.setNoFlow();// closed on top and bottom
	    }
    }

    // No source terms
    template <class Context>
    void source(RateVector& rate, const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    {
        rate = Scalar(0.0);
    }

private:
    DimMatrix K_;
    Scalar porosity_;
    Scalar temperature_;
    MaterialLawParams mat_;

    // initialize a random sequence that will return some normal distributed
    // number with up to 99% probability within the pertubation interval
    // random sources can be mutable; we expect them to behave erratically!
    mutable std::mt19937 rand_gen;
    mutable std::normal_distribution<double> norm_dist{0., PERTUBATION / 3.};
    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class FluidState, class Context>
    void initialFs(FluidState& fs, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        // get capillary pressure
        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // pressure; oleic phase is the reference
        fs.setPressure(oilPhaseIdx, 150*1e5);
        fs.setPressure(waterPhaseIdx, 150*1e5);
        fs.setPressure(gasPhaseIdx,150*1e5);

        // composition
        fs.setMoleFraction(oilPhaseIdx, Comp1Idx, 0.01);
        fs.setMoleFraction(oilPhaseIdx, Comp0Idx, 0.99);
        fs.setMoleFraction(oilPhaseIdx, Comp2Idx, 0.0);

        fs.setMoleFraction(gasPhaseIdx, Comp1Idx, 0.01);
        fs.setMoleFraction(gasPhaseIdx, Comp0Idx, 0.99);
        fs.setMoleFraction(gasPhaseIdx, Comp2Idx, 0.0);

        fs.setMoleFraction(waterPhaseIdx, Comp1Idx, 0.0);
        fs.setMoleFraction(waterPhaseIdx, Comp0Idx, 0.0);
        fs.setMoleFraction(waterPhaseIdx, Comp2Idx, 1.0);

        // temperature
        fs.setTemperature(temperature_);

        // saturation, oil-filled
        fs.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::waterPhaseIdx, 0.0);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);


        // fill in viscosity and enthalpy based on the state set above
        // and the fluid system defined in this class
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
                using CFRP = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/oilPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/true);

    }
};
} // namespace Opm

#endif
