#ifndef THREEPHASEFLUIDSYSTEM_HH
#define THREEPHASEFLUIDSYSTEM_HH

#include "components.hh"
#include "chiwoms.h"

#include <iostream>
#include <cassert>
#include <stdexcept>  // invalid_argument
#include <sstream>
#include <iostream>
#include <string>
#include <random>    // mt19937, normal_distribution
#include <limits>    // epsilon
#include <boost/format.hpp>  // boost::format

#include <opm/common/Exceptions.hpp>
#include <opm/material/IdealGas.hpp>

#include <opm/material/eos/PengRobinsonMixture.hpp>
#include <opm/material/eos/PengRobinsonParamsMixture.hpp>
#include "ChiParameterCache.hpp"

#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>


namespace Opm {
/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with brine and octane as the main components
 * in each their phase, and CO2 as solvent in both.
 */
template <class Scalar>
class ThreePhaseCo2OctaneBrineFluidSystem
        : public Opm::BaseFluidSystem<Scalar, ThreePhaseCo2OctaneBrineFluidSystem<Scalar> >
{
    using ThisType = ThreePhaseCo2OctaneBrineFluidSystem<Scalar>;
    using Base = Opm::BaseFluidSystem<Scalar, ThisType>;
    using PengRobinson =  typename Opm::PengRobinson<Scalar>;
    using PengRobinsonMixture =  typename Opm::PengRobinsonMixture<Scalar, ThisType>;

public:

    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    using ParameterCache = Opm::ChiParameterCache<Evaluation, ThisType>;

    /****************************************
         * Fluid phase related static parameters
         ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 3;

    //! Index of the liquid phase
    static const int oilPhaseIdx = 0;
    static const int gasPhaseIdx = 1;
    static const int waterPhaseIdx = 2;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {"o",  // oil phase
                                     "g",  // gas phase
                                     "w"}; // water phase

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
    using Octane = Opm::Octane<Scalar>;

    //! The component for pure solvent; since we have our own thermodynamic
    //! functions, we use the simple definition for the rest.
    using CO2 = Opm::ChiwomsCO2<Scalar>;

    //! The component for brine
    using PureH2O = Opm::H2O<Scalar>;
    using Brine = Opm::ChiwomsBrine<Scalar>;

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

        // if(phaseIdx == oilPhaseIdx) {
        if (phaseIdx == oilPhaseIdx || phaseIdx == gasPhaseIdx){
            // paramCache.updatePhase(fluidState, phaseIdx);
            auto dens = fluidState.averageMolarMass(phaseIdx)/paramCache.molarVolume(phaseIdx);
            return dens;
        }
            // return 650.;
            // return EOS::oleic_density(T, p, x);
        // } else if(phaseIdx == gasPhaseIdx) {
            // using IdealGas = Opm::IdealGas<Scalar>;
            // return IdealGas::density(LhsEval(molarMass(CO2Idx)), T, p);
        // }
        else {
            return 1000.;
            // return EOS::aqueous_density(T, p, x);
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

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& paramCache,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        if (phaseIdx == waterPhaseIdx || compIdx == BrineIdx)
            return 1.0;
        else {
            Scalar phi = Opm::getValue(PengRobinsonMixture::computeFugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx));
            return phi;
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
     * The values are given by the IVAR.
     */
    static Scalar interactionCoefficient(unsigned comp1Idx, unsigned comp2Idx)
    {
        // return 0;
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

};

};//namespace opm

#endif // THREEPHASEFLUIDSYSTEM_HH
