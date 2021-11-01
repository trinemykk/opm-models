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
#include "twophasefluidsystem.hh"
#include "ChiFlash.hpp"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/start.hh>

#include <opm/models/io/structuredgridvanguard.hh>

#include <opm/models/flash/flashmodel.hh>
#include <opm/models/immiscible/immisciblemodel.hh>

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>

// #include  <opm/simulators/linalg/parallelamgbackend.hh>
#include  <opm/simulators/linalg/parallelistlbackend.hh>

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
{ using type = Dune::YaspGrid<3>; };

template<class TypeTag, class MyTypeTag>
struct Temperature{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct Gravityfactor{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct SimulationName{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct EpisodeLength{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct Inflowrate{ using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct Initialpressure{ using type = UndefinedProperty; };

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
    using type = Opm::TwoPhaseThreeComponentFluidSystem<Scalar>;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::ChiwomsProblem>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum {
        oilPhaseIdx = FluidSystem::oilPhaseIdx,
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                           //    /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx, TODO
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
struct EnableGravity<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = true; };

template<class TypeTag>
struct Temperature<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 273.15 + TEMPERATURE;
};


template<class TypeTag>
struct Gravityfactor<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = GRAVITYFACTOR;
};

template<class TypeTag>
struct Inflowrate<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = INFLOW_RATE;
};

template<class TypeTag>
struct Initialpressure<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = MIN_PRES;
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

// template<class TypeTag>
// struct LinearSolverSplice<TypeTag, TTag::ChiwomsProblem>
// { using type = TTag::ParallelIstlLinearSolver; };

// template<class TypeTag>
// struct LinearSolverWrapper<TypeTag, TTag::ChiwomsProblem>
// { using type = Opm::Linear::SolverWrapperRestartedGMRes<TypeTag>; };

// template<class TypeTag>
// struct PreconditionerWrapper<TypeTag, TTag::ChiwomsProblem>
// { using type = Opm::Linear::PreconditionerWrapperILU<TypeTag>; };

// template<class TypeTag>
// struct PreconditionerOrder<TypeTag, TTag::ChiwomsProblem> { static constexpr int value = 2; };

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
template<class TypeTag>
struct VtkWriteLiquidMoleFractions<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteEquilibriumConstants<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = true; };

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
    static constexpr type value = X_SIZE;  // meter
};
template<class TypeTag>
struct CellsX<TypeTag, TTag::ChiwomsProblem> { static constexpr int value = NX; };
template<class TypeTag>
struct DomainSizeY<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = Y_SIZE;  // meter
};
template<class TypeTag>
struct CellsY<TypeTag, TTag::ChiwomsProblem> { static constexpr int value = NY; };
template<class TypeTag>
struct DomainSizeZ<TypeTag, TTag::ChiwomsProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = Z_SIZE;  //meter
};
template<class TypeTag>
struct CellsZ<TypeTag, TTag::ChiwomsProblem> { static constexpr int value = NZ; };

// compositional, with diffusion
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::ChiwomsProblem> { static constexpr bool value = false; };
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
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using Toolbox = Opm::MathToolbox<Evaluation>;
    using CoordScalar = typename GridView::ctype;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using H2O = typename Opm::H2O<Scalar>;
    using Brine = typename Opm::Brine<Scalar, H2O>;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };    
    enum { Comp1Idx = FluidSystem::Comp1Idx };
    enum { Comp0Idx = FluidSystem::Comp0Idx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + Comp1Idx };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableGravity = getPropValue<TypeTag, Properties::EnableGravity>() };

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using FullField = Dune::FieldMatrix<Scalar, NX, NY>;
    using FieldColumn = Dune::FieldVector<Scalar, NY>;
    using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;
    using FlashSolver = GetPropType<TypeTag, Properties::FlashSolver>;
    
    const unsigned XDIM = 0;
    const unsigned YDIM = 1;
    const unsigned ZDIM = 2;

    // influx on the left boundary
    Scalar rate;

public:
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem, enableEnergy>;
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

    void initGravity() {
        gravity_ = 0.0;
        gravityfactor_ = EWOMS_GET_PARAM(TypeTag, Scalar, Gravityfactor);
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dimWorld-1]  = gravityfactor_ * (-9.81);

    }

        /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * \param context Reference to the object which represents the
     *                current execution context.
     * \param spaceIdx The local index of spatial entity defined by the context
     * \param timeIdx The index used by the time discretization.
     */

    template <class Context>
    const DimVector& gravity(const Context& context OPM_UNUSED,
                             unsigned spaceIdx OPM_UNUSED,
                             unsigned timeIdx OPM_UNUSED) const
    { 
        return gravity(); 
    }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This method is used for problems where the gravitational
     * acceleration does not depend on the spatial position. The
     * default behaviour is that if the <tt>EnableGravity</tt>
     * property is true, \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$ holds,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     */
    const DimVector& gravity() const
    {
        return gravity_; 
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();
        // initialize fixed parameters; temperature, permeability, porosity
        initPetrophysics();

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
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Gravityfactor,
                             "The gravityfactor [-] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Inflowrate,
                             "The inflow rate [?] on the left boundary of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Initialpressure,
                             "The initial pressure [Pa s] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EpisodeLength,
                             "Time interval [s] for episode length");
    }

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
       // this->model().globalPhaseStorage(storageW, waterPhaseIdx);

        // Write mass balance information for rank 0
      //  if (this->gridView().comm().rank() == 0) {
	  //      std::cout << "Storage: oleic = [" << storageO << "], "
	  //                << "aqueous = [" << storageW << "]" << std::endl << std::flush;
      //  }
    }

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Evaluation, FluidSystem> fs;
        initialFs(fs, context, spaceIdx, timeIdx);
        values.assignNaive(fs);

    }

    // Constant temperature
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    {
        return temperature_;
    }

    // Constant gravityfactor
    template <class Context>
    Scalar gravityfactor(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    {
        return gravityfactor_;
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
    const MaterialLawParams& materialLawParams(const Context& context OPM_UNUSED,
                                               unsigned spaceIdx OPM_UNUSED, unsigned timeIdx OPM_UNUSED) const
    {
        return this->mat_;
    }

    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const Scalar eps = std::numeric_limits<double>::epsilon();
	const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onLeftBoundary_(pos)) {
            Scalar inflowrate = EWOMS_GET_PARAM(TypeTag, Scalar, Inflowrate);
            RateVector massRate(0.0);
            massRate = 0.0;
            massRate[contiCO2EqIdx] = inflowrate; // kg / (m^2 * s)

            values.setMassRate(massRate);
        }
        else if (onRightBoundary_(pos)) {
            Opm::CompositionalFluidState<Evaluation, FluidSystem> fs;
            initialFs(fs, context, spaceIdx, timeIdx);
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else
            values.setNoFlow();  // closed on top and bottom

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
    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[XDIM] < 1e-6; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[XDIM] > this->boundingBoxMax()[XDIM] - 1e-6; }

    bool onLowerBoundary_(const GlobalPosition& pos) const
    { return pos[ZDIM] < 1e-6; }

    bool onUpperBoundary_(const GlobalPosition& pos) const
    { return pos[ZDIM] > this->boundingBoxMax()[ZDIM] - 1e-6; }

    bool aboveMiddle_(const GlobalPosition& pos) const
    { return pos[ZDIM] >= (this->boundingBoxMax()[ZDIM] + this->boundingBoxMin()[ZDIM]) / 2; }

    bool leftMiddle(const GlobalPosition& pos) const
    { return pos[XDIM] > (this->boundingBoxMax()[XDIM] + this->boundingBoxMin()[XDIM]) / 2; }

    DimMatrix K_;
    Scalar porosity_;
    Scalar temperature_;
    Scalar gravityfactor_;
    MaterialLawParams mat_;
    DimVector gravity_;
    
    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class FluidState, class Context>
    void initialFs(FluidState& fs, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        // get capillary pressure
        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // pressure; set simple hydrostatic pressure initially.
        // OBS: If horizontal (NZ = 1), then h = Z_SIZE/2 since boundingBoxMax is cell edge
        Scalar init_pressure = EWOMS_GET_PARAM(TypeTag, Scalar, Initialpressure);
        bool enable_gravity = EWOMS_GET_PARAM(TypeTag, bool, EnableGravity);
        Scalar p_init;
        if (enable_gravity == true) {
            Scalar densityW = Brine::liquidDensity(temperature_, Scalar(init_pressure));
            const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
            Scalar h = this->boundingBoxMax()[ZDIM] - pos[ZDIM];
            p_init = (init_pressure*1e5) + densityW * h * (-this->gravity_[1]);
        }
        else
            p_init = init_pressure*1e5;
        fs.setPressure(oilPhaseIdx, p_init);
        fs.setPressure(gasPhaseIdx, p_init);

        // composition
        fs.setMoleFraction(oilPhaseIdx, Comp0Idx, MFCOMP0);
        fs.setMoleFraction(oilPhaseIdx, Comp1Idx, MFCOMP1); 

        fs.setMoleFraction(gasPhaseIdx, Comp0Idx, MFCOMP0);
        fs.setMoleFraction(gasPhaseIdx, Comp1Idx, MFCOMP1);
        
        // saturation, oil-filled
        fs.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

        // temperature
        fs.setTemperature(temperature_);

        // Density
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updatePhase(fs, oilPhaseIdx);
        paramCache.updatePhase(fs, gasPhaseIdx);
        fs.setDensity(oilPhaseIdx, FluidSystem::density(fs, paramCache, oilPhaseIdx));
        fs.setDensity(gasPhaseIdx, FluidSystem::density(fs, paramCache, gasPhaseIdx));

	// fill in viscosity and enthalpy based on the state set above
        // and the fluid system defined in this class. Oleic phase is the reference        
        // typename FluidSystem::template ParameterCache<Scalar> paramCache;
        //         using CFRP = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;
        // CFRP::solve(fs, paramCache,
        //             /*refPhaseIdx=*/oilPhaseIdx,
        //             /*setViscosity=*/true,
        //             /*setEnthalpy=*/false);

        if (enable_gravity == true) {
            // //
            // Run flash to get new density to correct pressure estimate
            // //
            // Set up z
            ComponentVector zInit(0.0);
            Scalar sumMoles = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar tmp = Opm::getValue(fs.molarity(phaseIdx, compIdx) * fs.saturation(phaseIdx));
                    zInit[compIdx] += Opm::max(tmp, 1e-8);
                    sumMoles += tmp;
                }
            }
            zInit /= sumMoles;
            
            // Flash solver setup
            Scalar flashTolerance = EWOMS_GET_PARAM(TypeTag, Scalar, FlashTolerance);
            int flashVerbosity = EWOMS_GET_PARAM(TypeTag, int, FlashVerbosity);
            std::string flashTwoPhaseMethod = EWOMS_GET_PARAM(TypeTag, std::string, FlashTwoPhaseMethod);
            int spatialIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

            // Set K and L initial
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Evaluation Ktmp = fs.wilsonK_(compIdx);
                    fs.setKvalue(compIdx, Ktmp);
            }
            const Evaluation& Ltmp = -1.0;
            fs.setLvalue(Ltmp);

            // Run flash solver
            FlashSolver::solve(fs, zInit, spatialIdx, flashVerbosity, flashTwoPhaseMethod, flashTolerance);

            // Calculate pressure again
            // if (fs.L(0) == 1)
            Evaluation densityL = fs.density(oilPhaseIdx);
            // if (fs.L(0) == 0)
            // Evaluation densityL = fs.density(gasPhaseIdx);
            const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
            Scalar h = this->boundingBoxMax()[ZDIM] - pos[ZDIM];
            p_init = (init_pressure*1e5) + Opm::getValue(densityL) * h * 9.81;
           fs.setPressure(oilPhaseIdx, p_init);
           fs.setPressure(gasPhaseIdx, p_init);
        }

    }
};
} // namespace Opm

#endif
