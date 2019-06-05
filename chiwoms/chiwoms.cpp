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
#include "config.h"

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/start.hh>

#include <ewoms/io/structuredgridvanguard.hh>

#include <ewoms/models/ncp/ncpmodel.hh>
#include <ewoms/models/immiscible/immisciblemodel.hh>

#include <ewoms/disc/ecfv/ecfvdiscretization.hh>

#include <ewoms/linear/parallelamgbackend.hh>

#include <opm/common/Exceptions.hpp>
#include <opm/material/IdealGas.hpp>

#include <opm/material/components/Component.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/Brine.hpp>
#include <opm/material/eos/PengRobinsonMixture.hpp>
#include <opm/material/eos/PengRobinsonParamsMixture.hpp>
#include "ChiParameterCache.hpp"


#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>

#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>
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

#include <opm/material/binarycoefficients/H2O_CO2.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>

//#include <dune/alugrid/grid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include "padua.h"
extern struct poly2d_t co2_c8_dens;
extern struct poly2d_t co2_c8_visc;
extern struct poly2d_t co2_c8_enth;
extern struct poly2d_t co2_h2o_dens;
extern struct poly2d_t co2_h2o_visc;
extern struct poly2d_t co2_h2o_enth;
#include "ewoms_compat.h"
#include "chiwoms.h"

namespace Ewoms {
/*!
 * \ingroup Components
 *
 * \brief A simple representation of linear octane
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Octane : public Opm::Component<Scalar, Octane<Scalar> >
{
public:
	/// Chemical name
	static const char* name() { return "C8"; }

	/// Molar mass in \f$\mathrm{[kg/mol]}\f$
	static Scalar molarMass() { return 0.11423; }

	/// Critical temperature in \f$\mathrm[K]}\f$
    static Scalar criticalTemperature() { return 568.7; }

	/// Critical pressure in \f$\mathrm[Pa]}\f$
    static Scalar criticalPressure() { return 2.49e6; }

	/// Acentric factor
	static Scalar acentricFactor() { return 0.398; }
};

template <class Scalar>
class ChiwomsCO2 : public Opm::SimpleCO2<Scalar>
{
public:
    /// Acentric factor
    static Scalar acentricFactor() { return 0.225; }
};

template <class Scalar>
class ChiwomsBrine : public Opm::H2O<Scalar>
{
public:
    /// Acentric factor
    static Scalar acentricFactor() { return 0.344; }
};

struct EOS
{
	template<typename LhsEval>
	static LhsEval eval(const struct poly2d_t & poly,
	                    LhsEval T, LhsEval p, LhsEval x) {
		LhsEval val;
		// if the linear solver overshoots and attempts to pick a solution
		// outside of the domain, just pretend that there is no change from the
		// boundary of the domain; hopefully the solver will regain its senses
		// and bring us back into the domain if we don't trip it up with
		// outlandish values.
		LhsEval fenced_p = Opm::min(Opm::max(p, poly.lo_x), poly.up_x);
		LhsEval fenced_x = Opm::min(Opm::max(x, poly.lo_y), poly.up_y);
        padua_eval(&poly, 1, &fenced_p, &fenced_x, &val);
		return val;
	}

	template<typename LhsEval>
	static LhsEval oleic_density(LhsEval T, LhsEval p, LhsEval x) {
		assert(T == (TEMPERATURE + 273.15));
		return eval(co2_c8_dens, T, p, x);
	}

	template<typename LhsEval>
	static LhsEval aqueous_density(LhsEval T, LhsEval p, LhsEval x) {
		assert(T == (TEMPERATURE + 273.15));
		return eval(co2_h2o_dens, T, p, x);
	}

	template<typename LhsEval>
	static LhsEval oleic_viscosity(LhsEval T, LhsEval p, LhsEval x) {
		assert(T == (TEMPERATURE + 273.15));
		return eval(co2_c8_visc, T, p, x);
	}

	template<typename LhsEval>
	static LhsEval aqueous_viscosity(LhsEval T, LhsEval p, LhsEval x) {
		assert(T == (TEMPERATURE + 273.15));
		return eval(co2_h2o_visc, T, p, x);
	}

	template<typename LhsEval>
	static LhsEval oleic_enthalpy(LhsEval T, LhsEval p, LhsEval x) {
		assert(T == (TEMPERATURE + 273.15));
		return eval(co2_c8_enth, T, p, x);
	}

	template<typename LhsEval>
	static LhsEval aqueous_enthalpy(LhsEval T, LhsEval p, LhsEval x) {
		assert(T == (TEMPERATURE + 273.15));
		return eval(co2_h2o_enth, T, p, x);
	}
};

} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup Fluidsystems
 *
 * \brief A one-phase fluid system with co2 and octane as components.
 */
template <class Scalar>
class OnePhaseCo2OctaneFluidSystem
    : public Opm::BaseFluidSystem<Scalar, OnePhaseCo2OctaneFluidSystem<Scalar> >
{
    typedef OnePhaseCo2OctaneFluidSystem<Scalar> ThisType;
    typedef Opm::BaseFluidSystem<Scalar, ThisType> Base;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    using ParameterCache = Opm::NullParameterCache<Evaluation>;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 1;

    //! Index of the liquid phase
    static const int onlyPhase = 0;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {"o"};  // oleic phase

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        // both CO2 and Octane have associative effects, separately and together
        return false;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;

    //! The component index of the oil; octane
    static const int OctaneIdx = 0;

    //! The component index of the solvent; co2
    static const int CO2Idx = 1;

    //! The component for pure oil
    typedef Ewoms::Octane<Scalar> Octane;

    //! The component for pure solvent
    typedef Opm::CO2<Scalar, Scalar> CO2;

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
        static const char* name[] = {
            Octane::name(),
            CO2::name()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        return (compIdx == OctaneIdx)
            ? Octane::molarMass()
            : (compIdx == CO2Idx)
            ? CO2::molarMass()
            : throw std::invalid_argument("Molar mass component index");
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalTemperature(unsigned compIdx)
    {
        return (compIdx == OctaneIdx)
            ? Octane::criticalTemperature()
            : (compIdx == CO2Idx)
            ? CO2::criticalTemperature()
            : throw std::invalid_argument("Critical temperature component index");
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalPressure(unsigned compIdx)
    {
        return (compIdx == OctaneIdx)
            ? Octane::criticalPressure()
            : (compIdx == CO2Idx)
            ? CO2::criticalPressure()
            : throw std::invalid_argument("Critical pressure component index");
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar acentricFactor(unsigned compIdx)
    {
        return (compIdx == OctaneIdx)
            ? Octane::acentricFactor()
            : (compIdx == CO2Idx)
            ? CO2::acentricFactor()
            : throw std::invalid_argument("Acentric factor component index");
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::density
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));

        return EOS::oleic_density(T, p, x);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));
        assert(T == (TEMPERATURE + 273.15));

        return EOS::oleic_viscosity(T, p, x);
    }

    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval enthalpy(const FluidState& fluidState,
                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                            unsigned phaseIdx)
    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));

        return EOS::oleic_enthalpy(T, p, x);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        // thermodynamic function represented by density
        const LhsEval rho = (compIdx == OctaneIdx ?
                             EOS::oleic_density(T, p, LhsEval(0.)) :  // x_CO2 == 0.
                             EOS::oleic_density(T, p, LhsEval(1.)));  // x_CO2 == 1.

        // molar mass
        const Scalar M = molarMass(compIdx);

        // fugacity coefficient calculated
        const LhsEval phi = rho * Opm::IdealGas<Scalar>::R * T / (M * p);

        return phi;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& /*fluidState*/,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned /*phaseIdx*/,
                                        unsigned /*compIdx*/)
    {
        return DIFFUSIVITY;
    }
};

/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with brine and octane as the main components
 * in each their phase, and CO2 as solvent in both.
 */
template <class Scalar>
class TwoPhaseCo2OctaneBrineFluidSystem
	: public Opm::BaseFluidSystem<Scalar, TwoPhaseCo2OctaneBrineFluidSystem<Scalar> >
{
	typedef TwoPhaseCo2OctaneBrineFluidSystem<Scalar> ThisType;
	typedef Opm::BaseFluidSystem<Scalar, ThisType> Base;
    typedef typename Opm::PengRobinson<Scalar> PengRobinson;
    typedef typename Opm::PengRobinsonMixture<Scalar, ThisType> PengRobinsonMixture;

public:
	//! \copydoc BaseFluidSystem::ParameterCache
    //template <class Evaluation>
    //using ParameterCache = Opm::NullParameterCache<Evaluation>;

    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    using ParameterCache = Opm::ChiParameterCache<Evaluation, ThisType>;

    /****************************************
	 * Fluid phase related static parameters
	 ****************************************/

	//! \copydoc BaseFluidSystem::numPhases
	static const int numPhases = 2;

	//! Index of the liquid phase
    static const int oilPhaseIdx = 0;
    static const int waterPhaseIdx = 1;

	//! \copydoc BaseFluidSystem::phaseName
	static const char* phaseName(unsigned phaseIdx)
	{
		static const char* name[] = {"o",  // oleic phase
		                             "w"}; // aquous phase

		assert(0 <= phaseIdx && phaseIdx < numPhases);
		return name[phaseIdx];
	}

	//! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned phaseIdx)
	{
        if (phaseIdx == oilPhaseIdx)
            return false;

        // CO2 have associative effects, both with octane and with brine
        return true;
	}


	/****************************************
	 * Component related static parameters
	 ****************************************/

	//! \copydoc BaseFluidSystem::numComponents
	static const int numComponents = 3;  // octane, co2 and brine

	//! The component index of the oil; octane
	static const int OctaneIdx = 0;

	//! The component index of the solvent; co2
	static const int CO2Idx = 1;

	//! The component index of the other main component; h2o + salts
	static const int BrineIdx = 2;

	//! The component for pure oil
	typedef Ewoms::Octane<Scalar> Octane;

	//! The component for pure solvent; since we have our own thermodynamic
	//! functions, we use the simple definition for the rest.
    typedef Ewoms::ChiwomsCO2<Scalar> CO2;

	//! The component for brine
	typedef Opm::H2O<Scalar> PureH2O;
    typedef Ewoms::ChiwomsBrine<Scalar> Brine;

    static void init(Scalar minT = 273.15,
                     Scalar maxT = 373.15,
                     Scalar minP = 1e4,
                     Scalar maxP = 100e6)
    {
        Opm::PengRobinsonParamsMixture<Scalar, ThisType, oilPhaseIdx, /*useSpe5=*/true> prParams;

        // find envelopes of the 'a' and 'b' parameters for the range
        // minT <= T <= maxT and minP <= p <= maxP. For
        // this we take advantage of the fact that 'a' and 'b' for
        // mixtures is just a convex combination of the attractive and
        // repulsive parameters of the pure components

        Scalar minA = 1e30, maxA = -1e30;
        Scalar minB = 1e30, maxB = -1e30;

        prParams.updatePure(minT, minP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, minP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(minT, maxP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, maxP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };
      //  PengRobinson::init(/*aMin=*/minA, /*aMax=*/maxA, /*na=*/100,
      //                     /*bMin=*/minB, /*bMax=*/maxB, /*nb=*/200);
    }

	//! \copydoc BaseFluidSystem::componentName
	static const char* componentName(unsigned compIdx)
	{
		static const char* name[] = {
			Octane::name(),
			CO2::name(),
			Brine::name(),
		};

		assert(0 <= compIdx && compIdx < numComponents);
		return name[compIdx];
	}

	//! \copydoc BaseFluidSystem::molarMass
	static Scalar molarMass(unsigned compIdx)
	{
		return (compIdx == OctaneIdx)
			? Octane::molarMass()
			: (compIdx == CO2Idx)
			? CO2::molarMass()
			: (compIdx == BrineIdx)
			? Brine::molarMass()
			: throw std::invalid_argument("Molar mass component index");
	}

	/*!
	 * \brief Critical temperature of a component [K].
	 *
	 * \copydetails Doxygen::compIdxParam
	 */
	static Scalar criticalTemperature(unsigned compIdx)
	{
		return (compIdx == OctaneIdx)
			? Octane::criticalTemperature()
			: (compIdx == CO2Idx)
			? CO2::criticalTemperature()
			: (compIdx == BrineIdx)
			? Brine::criticalTemperature()
			: throw std::invalid_argument("Critical temperature component index");
	}

	/*!
	 * \brief Critical pressure of a component [Pa].
	 *
	 * \copydetails Doxygen::compIdxParam
	 */
	static Scalar criticalPressure(unsigned compIdx)
	{
		return (compIdx == OctaneIdx)
			? Octane::criticalPressure()
			: (compIdx == CO2Idx)
			? CO2::criticalPressure()
			: (compIdx == BrineIdx)
			? Brine::criticalPressure()
			: throw std::invalid_argument("Critical pressure component index");
	}

	/*!
	 * \brief The acentric factor of a component [].
	 *
	 * \copydetails Doxygen::compIdxParam
	 */
	static Scalar acentricFactor(unsigned compIdx)
	{
		return (compIdx == OctaneIdx)
			? Octane::acentricFactor()
			: (compIdx == CO2Idx)
			? CO2::acentricFactor()
			: (compIdx == BrineIdx)
			? Brine::acentricFactor()
			: throw std::invalid_argument("Molar mass component index");
	}

	/****************************************
	 * thermodynamic relations
	 ****************************************/

	/*!
	 * \copydoc BaseFluidSystem::density
	 */
	template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
	static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& paramCache,
	                       unsigned phaseIdx)
	{
		assert(0 <= phaseIdx && phaseIdx < numPhases);

        //return fluidState.averageMolarMass(phaseIdx)/paramCache.molarVolume(phaseIdx);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));
        assert(T == (TEMPERATURE + 273.15));

        if(phaseIdx == oilPhaseIdx) {
            //return 650.;
            return EOS::oleic_density(T, p, x);
        }
        else {
            //return 1000.;
            return EOS::aqueous_density(T, p, x);
        }
	}

	//! \copydoc BaseFluidSystem::viscosity
	template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
	static LhsEval viscosity(const FluidState& fluidState,
	                         const ParameterCache<ParamCacheEval>& /*paramCache*/,
	                         unsigned phaseIdx)
	{
		assert(0 <= phaseIdx && phaseIdx < numPhases);

		const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
		const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
		const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));
		assert(T == (TEMPERATURE + 273.15));

        if(phaseIdx == oilPhaseIdx) {
			return EOS::oleic_viscosity(T, p, x);
		}
		else {
			return EOS::aqueous_viscosity(T, p, x);
		}
	}

	//! \copydoc BaseFluidSystem::enthalpy
	template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
	static LhsEval enthalpy(const FluidState& fluidState,
	                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
	                        unsigned phaseIdx)
	{
		const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
		const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
		const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));

        if(phaseIdx == oilPhaseIdx) {
			return EOS::oleic_enthalpy(T, p, x);
		}
		else {
			return EOS::aqueous_enthalpy(T, p, x);
		}
	}

	// according to <https://srdata.nist.gov/solubility/sol_detail.aspx?sysID=38_103>
	// the solubility of octane in aqueous phase is 2.0g/100g sln. since octane
	// (C8H18) has a molecular weight of 114.23 g/mol and water (H2O) of 18.01528 g/mol
	// we have a total of 2.0g/114.23 g/mol ~= 0.0175 mol of octane and
	// (100g-2.0g)/18.01528 g/mol ~= 5.44 mol of water, for a total of 5.45 mol,
	// of which the mole fraction of octane is 0.0175/5.45 ~= 3.2e-3
    //constexpr static Scalar sol_aqueous_oil = 3.208e-3;  // solution of octane in brine

	// the solubility of water in the oleic (octane) phase is according the same
	// reference above, 7.3g/100g sln, giving 7.3g/18.01528 g/mol ~= 0.41 mol of
	// water and (100g-7.3g)/114.23 g/mol ~= 0.81 mol of octane, in a 100 g solution,
	// yielding a mole fraction of 0.41/(0.41 + 0.81) = 0.33 for water.
    //constexpr static Scalar sol_oleic_water = sol_aqueous_oil; // 3.330e-1;

	// partition coefficients when both oleic and aqueous phases are saturated with
	// maximum dissoluted amount of the other component. these coefficients should
	// give fugacities with maximum volatility of the component in its non-native phase
    //constexpr static Scalar k_aqueous_oil = (1 - sol_oleic_water) / sol_aqueous_oil;  // ~200
    //constexpr static Scalar k_oleic_water = (1 - sol_aqueous_oil) / sol_oleic_water;  // ~3


	//! \copydoc BaseFluidSystem::fugacityCoefficient
	template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
	static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& paramCache,
	                                   unsigned phaseIdx,
	                                   unsigned compIdx)
	{
		assert(0 <= phaseIdx && phaseIdx < numPhases);
		assert(0 <= compIdx && compIdx < numComponents);        

        if (phaseIdx == oilPhaseIdx) {
#if 0
#warning HACK
            if (compIdx == BrineIdx) {
                return 1e7/fluidState.pressure(oilPhaseIdx);

            }else if (compIdx == OctaneIdx)
                return 40e3/fluidState.pressure(oilPhaseIdx);
            else
               return 500e3/fluidState.pressure(oilPhaseIdx);
#else
            if (compIdx == BrineIdx) {
                return 1e7/fluidState.pressure(oilPhaseIdx);

            }else if (compIdx == CO2Idx)
                return 500e3/fluidState.pressure(oilPhaseIdx);
            else {
            return PengRobinsonMixture::computeFugacityCoefficient(fluidState,
                                                                   paramCache,
                                                                   phaseIdx,
                                                                   compIdx);
            }
#endif
        } else {
            assert(phaseIdx == waterPhaseIdx);
            //return PengRobinsonMixture::computeFugacityCoefficient(fluidState,
            //                                                       paramCache,
            //                                                       phaseIdx,
            //                                                      compIdx);
            //if (compIdx == OctaneIdx){
            //    return 0.00005;
            //}
            //return 1;
            //return henryCoeffWater_(compIdx, fluidState.temperature(waterPhaseIdx))
            //  / fluidState.pressure(waterPhaseIdx);
            if (compIdx == OctaneIdx) {
                return 10e6/fluidState.pressure(waterPhaseIdx);

            }else if (compIdx == BrineIdx)
                return 80e3/fluidState.pressure(waterPhaseIdx);
            else
               return 1e6/fluidState.pressure(waterPhaseIdx);
        }

        throw std::invalid_argument("crap!");
	}

	//! \copydoc BaseFluidSystem::diffusionCoefficient
	template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
	static LhsEval diffusionCoefficient(const FluidState& /*fluidState*/,
	                                    const ParameterCache<ParamCacheEval>& /*paramCache*/,
	                                    unsigned /*phaseIdx*/,
	                                    unsigned /*compIdx*/)
	{
		return DIFFUSIVITY;
	}

    /*!
     * \brief Returns the interaction coefficient for two components.
     *
     * The values are given by the SPE5 paper.
     */
    static Scalar interactionCoefficient(unsigned comp1Idx, unsigned comp2Idx)
    {
        unsigned i = std::min(comp1Idx, comp2Idx);
        unsigned j = std::max(comp1Idx, comp2Idx);
        #warning interactionCoefficients from Ivar
        if (i == OctaneIdx && j == CO2Idx)
            return 0.1089;
        else if (i == OctaneIdx && j == BrineIdx)
            return 1.1290;
        else if (i == CO2Idx && j == BrineIdx)
            return -0.0736;
        return 0;
    }

    template <class LhsEval>
    static LhsEval henryCoeffWater_(unsigned compIdx, const LhsEval& temperature)
    {

        // use henry's law for the solutes and the vapor pressure for
        // the solvent.
        switch (compIdx) {
        case BrineIdx: return Brine::vaporPressure(temperature);

            // the values of the Henry constant for the solutes have
            // are faked so far...
#warning set correct Henry constants
        case CO2Idx: return 5.57601e+09;
        case OctaneIdx: return 1e10;
        default: throw std::logic_error("Unknown component index "+std::to_string(compIdx));
        }
    }
};

} // namespace Opm

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
        if (j < NY/2) {
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
	        (this->boundingBoxMax()[YDIM] - this->boundingBoxMin()[YDIM])/2.)) {
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

BEGIN_PROPERTIES

NEW_TYPE_TAG(ChiwomsNcpEcfvProblem, INHERITS_FROM(NcpModel, ChiwomsProblem));
SET_TAG_PROP(ChiwomsNcpEcfvProblem, SpatialDiscretizationSplice, EcfvDiscretization);
SET_TAG_PROP(ChiwomsNcpEcfvProblem, LocalLinearizerSplice, AutoDiffLocalLinearizer);

END_PROPERTIES

int main(int argc, char **argv)
{
    typedef TTAG(ChiwomsNcpEcfvProblem) EcfvProblemTypeTag;
    return Ewoms::start<EcfvProblemTypeTag>(argc, argv);
}
