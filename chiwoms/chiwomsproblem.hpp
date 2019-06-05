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
#include "onephasefluidsystem.hh"
#include "twophasefluidsystem.hh"

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/start.hh>

#include <ewoms/io/structuredgridvanguard.hh>

#include <ewoms/models/ncp/ncpmodel.hh>
#include <ewoms/models/immiscible/immisciblemodel.hh>

#include <ewoms/disc/ecfv/ecfvdiscretization.hh>

#include <ewoms/linear/parallelamgbackend.hh>

#include <opm/common/Exceptions.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Unused.hpp>

//#include <dune/alugrid/grid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>


namespace Ewoms {
template <class TypeTag>
class ChiwomsProblem;
} // namespace Ewoms

BEGIN_PROPERTIES

NEW_TYPE_TAG(ChiwomsProblem);

SET_TYPE_PROP(ChiwomsProblem, Grid, Dune::YaspGrid<2>);
//SET_TYPE_PROP(ChiwomsProblem, Grid,
//              Dune::ALUGrid</*dim=*/2,
//                            /*dimWorld=*/2,
//                            Dune::cube,
//                            Dune::nonconforming>);

// declare the CO2 finger problem specific property tags
NEW_PROP_TAG(TopResvPres);
NEW_PROP_TAG(Temperature);
NEW_PROP_TAG(SimulationName);
NEW_PROP_TAG(EpisodeLength);
NEW_PROP_TAG(WaveLength);

NEW_PROP_TAG(ConnateWater);
NEW_PROP_TAG(ResidualOil);
NEW_PROP_TAG(EntryPressure);
NEW_PROP_TAG(PoreSizeDist);


// Set the problem property
SET_TYPE_PROP(ChiwomsProblem, Problem,
              Ewoms::ChiwomsProblem<TypeTag>);

// Set fluid configuration
SET_PROP(ChiwomsProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Ewoms::TwoPhaseCo2OctaneBrineFluidSystem<Scalar> type;
};

// Set the material Law
SET_PROP(ChiwomsProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum {
        oilPhaseIdx = FluidSystem::oilPhaseIdx,
        waterPhaseIdx = FluidSystem::waterPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<
	    Scalar,
        /*wettingPhaseIdx=*/ FluidSystem::waterPhaseIdx,
        /*nonWettingPhaseIdx=*/ FluidSystem::oilPhaseIdx> Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedBrooksCorey<Traits> EffMaterialLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::EffToAbsLaw<EffMaterialLaw> type;
};

// Write the Newton convergence behavior to disk?
SET_BOOL_PROP(ChiwomsProblem, NewtonWriteConvergence, false);

// Enable gravity
SET_BOOL_PROP(ChiwomsProblem, EnableGravity, true);

// set the defaults for the problem specific properties
SET_SCALAR_PROP(ChiwomsProblem, TopResvPres, 200 * 1.e5);
SET_SCALAR_PROP(ChiwomsProblem, Temperature, 273.15 + TEMPERATURE);
SET_STRING_PROP(ChiwomsProblem, SimulationName, "chiwoms");

// The default for the end time of the simulation
SET_SCALAR_PROP(ChiwomsProblem, EndTime, SIM_TIME * 24. * 60. * 60.);

// convergence control
SET_SCALAR_PROP(ChiwomsProblem, InitialTimeStepSize, 30);
SET_SCALAR_PROP(ChiwomsProblem, LinearSolverTolerance, 1e-3);
SET_SCALAR_PROP(ChiwomsProblem, LinearSolverAbsTolerance, 0.);
SET_SCALAR_PROP(ChiwomsProblem, NewtonTolerance, 1e-3);
SET_SCALAR_PROP(ChiwomsProblem, MaxTimeStepSize, 60 * 60);
SET_INT_PROP(ChiwomsProblem, NewtonMaxIterations, 10);
SET_INT_PROP(ChiwomsProblem, NewtonTargetIterations, 6);

// output
SET_BOOL_PROP(ChiwomsProblem, VtkWriteFilterVelocities, true);
SET_BOOL_PROP(ChiwomsProblem, VtkWritePotentialGradients, true);
SET_BOOL_PROP(ChiwomsProblem, VtkWriteTotalMassFractions, true);
SET_BOOL_PROP(ChiwomsProblem, VtkWriteTotalMoleFractions, true);
SET_BOOL_PROP(ChiwomsProblem, VtkWriteFugacityCoeffs, true);

// write restart for every hour
SET_SCALAR_PROP(ChiwomsProblem, EpisodeLength, 60. * 60.);

// mesh grid
SET_TYPE_PROP(ChiwomsProblem, Vanguard, Ewoms::StructuredGridVanguard<TypeTag>);
SET_SCALAR_PROP(ChiwomsProblem, DomainSizeX, (X_SIZE / 100.));
SET_INT_PROP(ChiwomsProblem, CellsX, NX);
SET_SCALAR_PROP(ChiwomsProblem, DomainSizeY, (Y_SIZE / 100.));
SET_INT_PROP(ChiwomsProblem, CellsY, NY);
SET_SCALAR_PROP(ChiwomsProblem, DomainSizeZ, 1.);
SET_INT_PROP(ChiwomsProblem, CellsZ, 1);

// compositional, with diffusion
SET_BOOL_PROP(ChiwomsProblem, EnableEnergy, false);
SET_BOOL_PROP(ChiwomsProblem, EnableDiffusion, true);

// injection rate parameter
SET_SCALAR_PROP(ChiwomsProblem, WaveLength, WAVE_LENGTH);

// default hydrology: almost (but not quite) Snohvit-like Brooks-Corey
SET_SCALAR_PROP(ChiwomsProblem, ConnateWater, 0.15);
SET_SCALAR_PROP(ChiwomsProblem, ResidualOil, 0.10);
SET_SCALAR_PROP(ChiwomsProblem, EntryPressure, 1.5e4 /* [Pa] */);
SET_SCALAR_PROP(ChiwomsProblem, PoreSizeDist, 2.0);

END_PROPERTIES

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 */
template <class TypeTag>
class ChiwomsProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { numPhases = FluidSystem::numPhases };
    enum {
        oilPhaseIdx = FluidSystem::oilPhaseIdx,
        waterPhaseIdx = FluidSystem::waterPhaseIdx
    };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { OctaneIdx = FluidSystem::OctaneIdx };
    enum { BrineIdx = FluidSystem::BrineIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + CO2Idx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldMatrix<Scalar, NX, NY> FullField;
    typedef Dune::FieldVector<Scalar, NY> FieldColumn;

    const unsigned XDIM = 0;
    const unsigned YDIM = 1;

    // we need to layout the initial state from top to down, so fill
    // a full matrix with that now and access it (possibly) randomly later
    // note that these matrices are logical to fill from the top, so index
    // 0 is the *top-most* element, and the index increase as we move *down*
    // the column. this is opposite of the orientation of the grid in eWoms
    FieldColumn init_pres;  // oleic phase pressure; ref. pres. w/o cap. pres

    // background level of CO2 and Octane that will always be present initially
    const Scalar co2_tracer = 0.;
    const Scalar c8_tracer = 0.;

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

    // find pressure at the bottom of a cell
    Scalar bottomPressure(Scalar pres_top, unsigned phaseIdx) {
	    Scalar grav = -this->gravity()[YDIM];
	    Scalar tiny = std::numeric_limits<double>::epsilon();
	    Scalar comp = 0.; // non-solvent filled, either oil or brine but not CO2

	    // since we have a regular grid, then we can find the height of each
	    // cell outside of the loop. divide by hundred to get to SI unit.
	    Scalar cell_height = (Y_SIZE / 100.) / NY;

	    // find density at the top, using top pressure, use this as the
	    // average of the cell initially.
        Scalar rho_top = (phaseIdx == oilPhaseIdx)
		    ? EOS::oleic_density(temperature_, pres_top, comp)
		    : EOS::aqueous_density(temperature_, pres_top, comp);

	    Scalar rho_avg = rho_top;

	    Scalar pres_btm, rho_btm, rho_old = 0;
	    while(std::fabs(rho_old - rho_avg) > tiny) {
		    // save for comparison next step
		    rho_old = rho_avg;

		    // estimate bottom boundary pressure by the height of this
		    // cell and the density
		    pres_btm = pres_top + cell_height * grav * rho_avg;

		    // get the density at the bottom, using its pressure. the
		    // cell still has the same mole fraction
            rho_btm = (phaseIdx == oilPhaseIdx)
			    ? EOS::oleic_density(temperature_, pres_top, comp)
			    : EOS::aqueous_density(temperature_, pres_top, comp);

		    // calculate new average density for the cell
		    rho_avg = (rho_top + rho_btm) / 2.;
	    }
	    return pres_btm;
    }

    void initPressure() {
        // start at the top of each column, working downwards calculating the
        // pressure based on the weight above
        Scalar topResvPres = EWOMS_GET_PARAM(TypeTag, Scalar, TopResvPres);

        Scalar pres_top = topResvPres;
        for(unsigned j = 0; j < NY; ++j) {
	        // first half is filled with oil, next half is filled with water
	        Scalar pres_btm = bottomPressure(
                pres_top, (j < NY/2) ? oilPhaseIdx : waterPhaseIdx);

	        // at this point we have a stable average density; now let's
	        // turn this into an average pressure for the cell
	        init_pres[j] = (pres_btm + pres_top) / 2.;

	        // prepare for the next iteration; one cell lower
	        pres_top = pres_btm;
        }
    }

    void initRate() {
	    // always the same temperature in isothermal model
	    Scalar temp = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);

	    // get the pressure condition at the site; ignore that the pressure
	    // increase downwards, as we don't want spatially varying rate
	    Scalar pres = EWOMS_GET_PARAM(TypeTag, Scalar, TopResvPres);

	    // look up viscosity at this pressure
	    Scalar co2_visc = EOS::oleic_viscosity(temp, pres, 1.); // x_CO2 = 1., x_C8 = 0.
	    Scalar c8_visc = EOS::oleic_viscosity(temp, pres, 0.); // x_CO2 = 0., x_C8 = 1.

	    // mobility ratio of co2 going into octane
	    Scalar mob_rate = c8_visc / co2_visc;

	    // critical wavelength for fingers to appear. we set this to 1/10 of
	    // the *height* of the domain; distance across the inflow vertical
	    // boundary.
	    Scalar wav_len = EWOMS_GET_PARAM(TypeTag, Scalar, WaveLength);
	    Scalar lambda_c = (this->boundingBoxMax()[YDIM] -
	                       this->boundingBoxMin()[YDIM]) * wav_len;

	    // velocity is u_x - u_c, where u_c is the critical flux, which we
	    // here set to zero
	    Scalar velo = 4. * M_PI * (mob_rate + 1.) / (mob_rate - 1.) *
		    DIFFUSIVITY / lambda_c;

	    // look up density at this pressure
	    Scalar co2_dens = EOS::oleic_density(temp, pres, 1.); // x_CO2 = 1.

	    // influx
	    this->rate = velo * co2_dens; // [m/s * kg/m3 = kg/(s * m^2)

	    // show rate as distance from a well with known injection rate
	    // of 1 Mt/yr (Snohvit is 0.7 Mt/yr)
	    Scalar inj = 1; /* Mt/yr */
	    Scalar scale = 1.e9 /* kg/Mt */ / (365.25 * 24. * 60. * 60.) /* secs/yr */;
	    Scalar height = Y_SIZE / 100.; /* meter high plume tongue */
	    Scalar radius = (inj * scale / rate) / (2. * M_PI * height);
	    if(this->simulator().gridView().comm().rank() == 0) {
		    std::cout << "Rate is equivalent to the front of a "
		              << boost::format("%4.2f") % height << " m thick plume "
		              << "from a " << boost::format("%1.0f") % inj << " Mt/yr well "
		              << "at " << boost::format("%4.0f") % radius << " m distance."
		              << std::endl;
	    }
    }

    void initHydrology() {
        this->mat_.setResidualSaturation(oilPhaseIdx,
                                         0);
        this->mat_.setResidualSaturation(waterPhaseIdx,
                                         0);
        this->mat_.setEntryPressure(0);
	    this->mat_.setLambda(EWOMS_GET_PARAM(TypeTag, Scalar, PoreSizeDist));
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

        // initialize column with hydrostatic pressure
        initPressure();

        // calculate rate up front; we don't vary it with time and place
        initRate();
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
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        const double dy = (Y_SIZE / 100.) / NY;

        // since the initialization matrices are turned upside-down compared
        // to the grid, we must subtract the index from the top to find the
        // appropriate index in the initialization matrix from the coordinate
        // (i, j) is thus only valid for lookup in init_comp and init_pres
        unsigned j = (NY-1) - static_cast<unsigned>(std::floor(
                    (pos[YDIM] - this->boundingBoxMin()[YDIM]) / dy));

        // get capillary pressure
        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // pressure; oleic phase is the reference
        fs.setPressure(oilPhaseIdx, 20e6);
        fs.setPressure(waterPhaseIdx, 20e6);

        // composition
        Scalar sol_oleic_water = 0;
        Scalar sol_aqueous_oil = 0;
        fs.setMoleFraction(oilPhaseIdx, CO2Idx, 0.049);
        fs.setMoleFraction(oilPhaseIdx, OctaneIdx, 0.95);
        fs.setMoleFraction(oilPhaseIdx, BrineIdx, 0.001);

        fs.setMoleFraction(waterPhaseIdx, CO2Idx, 0.0000);
        fs.setMoleFraction(waterPhaseIdx, OctaneIdx, 0.0000);
        fs.setMoleFraction(waterPhaseIdx, BrineIdx, 1.0);

        // temperature
        fs.setTemperature(temperature_);

        // saturation, oil-filled at upper half, water-filled at lower half
        // DOES NOT CONV WITH fugacity != 1 when co2 is injected
        // conv with fugacity = 1 (no mixing) when co2 is injected
        //fs.setSaturation(FluidSystem::oilPhaseIdx, (j < NY/2) ? 0.8 : 0.2);
        //fs.setSaturation(FluidSystem::waterPhaseIdx, (j < NY/2) ? 0.2 : 0.8);

        // saturation, oil-filled at upper half, water-filled at lower half
        if (true || (j < NY/2)) {
            fs.setSaturation(FluidSystem::oilPhaseIdx, 1.);
            fs.setSaturation(FluidSystem::waterPhaseIdx, 0.);
        }
        else {
            fs.setSaturation(FluidSystem::oilPhaseIdx, 0.);
            fs.setSaturation(FluidSystem::waterPhaseIdx, 1.);
        }

        // saturation, oil-filled
        // conv with small interaction when when co2 is injected
        // conv with fugacity = 1 (no mixing) when co2 is injected
        //fs.setSaturation(FluidSystem::oilPhaseIdx, 1.);
        //fs.setSaturation(FluidSystem::waterPhaseIdx, 0);

        // saturation, water-filled
        // DOES NOT CONV WITH fugacity != 1 co2 injected
        // does conv with fugacity != 1 , co2 injected
     //  fs.setSaturation(FluidSystem::oilPhaseIdx, (j < NY/2) ? 0. : 0.);
     //  fs.setSaturation(FluidSystem::waterPhaseIdx, (j < NY/2) ? 1. : 1.);


        // fill in viscosity and enthalpy based on the state set above
        // and the fluid system defined in this class
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        if (fs.saturation(waterPhaseIdx) > .1){
            CFRP::solve(fs, paramCache,
                        /*refPhaseIdx=*/ waterPhaseIdx,
                        /*setViscosity=*/ true,
                        /*setEnthalpy=*/ true);
        } else {
          CFRP::solve(fs, paramCache,
                      /*refPhaseIdx=*/ oilPhaseIdx,
                      /*setViscosity=*/ true,
                      /*setEnthalpy=*/ true);
        }

        std::cout << "hello" << std::endl;
        // const auto& matParams = this->materialLawParams(context, spaceIdx,
        // timeIdx);
        //values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        values.assignNaive(fs);
    }

    /// Constant temperature
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; }

    /// Constant permeability
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context OPM_UNUSED,
                                           unsigned spaceIdx OPM_UNUSED,
                                           unsigned timeIdx OPM_UNUSED) const
    { return K_; }

    /// Constant porosity
    template <class Context>
    Scalar porosity(const Context& context  OPM_UNUSED,
                    unsigned spaceIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED) const
    { return porosity_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    { return this->mat_; }

    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        //values.setNoFlow();
        //return;

        const Scalar eps = std::numeric_limits<double>::epsilon();
	    const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

	    // left side has a fixed inflow rate in the lower half
	    if((pos[XDIM] < this->boundingBoxMin()[XDIM] + eps) &&
	       (pos[YDIM] < this->boundingBoxMin()[YDIM] +
                (this->boundingBoxMax()[YDIM] - this->boundingBoxMin()[YDIM]))) {
		    // create some noise in the injection rate
            const double noise = 0;//this->norm_dist(this->rand_gen);

		    // assign rate to the CO2 component of the inflow
		    RateVector massRate(0.);
            massRate[contiCO2EqIdx] =  -1e-7;
		    values.setMassRate(massRate);
	    }
#if 0
	    else if(pos[XDIM] > this->boundingBoxMax()[XDIM] - eps) {
		    // get mole fraction in incident element
		    Scalar x_co2_in_oil = Toolbox::value(
			    context.intensiveQuantities(spaceIdx, timeIdx).fluidState()
                .moleFraction(oilPhaseIdx, CO2Idx));
		    Scalar x_co2_in_brine = Toolbox::value(
			    context.intensiveQuantities(spaceIdx, timeIdx).fluidState()
                .moleFraction(waterPhaseIdx, CO2Idx));
		    Scalar x_c8 = Toolbox::value(
			    context.intensiveQuantities(spaceIdx, timeIdx).fluidState()
                .moleFraction(oilPhaseIdx, OctaneIdx));
		    Scalar x_h2o = Toolbox::value(
			    context.intensiveQuantities(spaceIdx, timeIdx).fluidState()
                .moleFraction(waterPhaseIdx, BrineIdx));


		    Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
		    // saturation, only one phase for the time being
            fs.setSaturation(FluidSystem::oilPhaseIdx, 1.);
            fs.setSaturation(FluidSystem::waterPhaseIdx, 0.);

		    // copy mole fraction over
            fs.setMoleFraction(oilPhaseIdx, CO2Idx, x_co2_in_oil);
            fs.setMoleFraction(oilPhaseIdx, OctaneIdx, x_c8);
            fs.setMoleFraction(waterPhaseIdx, CO2Idx, x_co2_in_brine);
            fs.setMoleFraction(waterPhaseIdx, BrineIdx, x_h2o);

		    const double dy = (Y_SIZE / 100.) / NY;
		    // since the initialization matrices are turned upside-down compared
		    // to the grid, we must subtract the index from the top to find the
		    // appropriate index in the initialization matrix from the coordinate
		    // (i, j) is thus only valid for lookup in init_comp and init_pres
		    unsigned j = (NY-1) - static_cast<unsigned>(std::floor(
                    (pos[YDIM] - this->boundingBoxMin()[YDIM]) / dy));

		    // get capillary pressure
		    Scalar pC[numPhases];
		    const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
		    MaterialLaw::capillaryPressures(pC, matParams, fs);

		    // hydrostatic pressure in outlet
            fs.setPressure(oilPhaseIdx, this->init_pres[j]);
            fs.setPressure(waterPhaseIdx, this->init_pres[j] + (pC[waterPhaseIdx] - pC[oilPhaseIdx]));
		    fs.setTemperature(this->temperature_);
		    fs.checkDefined();

		    // update dependent quantities, such as viscosities
		    typename FluidSystem::template ParameterCache<Scalar> paramCache;
		    typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
		    CFRP::solve(fs, paramCache,
                        /*refPhaseIdx=*/oilPhaseIdx,
		                /*setViscosity=*/true,
		                /*setEnthalpy=*/true);

		    values.setFreeFlow(context, spaceIdx, timeIdx, fs);
	    }
#endif
	    else {
		    // closed on top and bottom
		    values.setNoFlow();
	    }
    }

    /// No source terms
    template <class Context>
    void source(RateVector& rate, const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    { rate = Scalar(0.0); }

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

};
} // namespace Ewoms

#endif
