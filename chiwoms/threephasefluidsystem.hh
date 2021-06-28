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


namespace Opm {
/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with brine and octane as the main components
 * in each their phase, and CO2 as solvent in both.
 */
template <class Scalar>
class ThreePhaseThreeComponentFluidSystem
        : public Opm::BaseFluidSystem<Scalar, ThreePhaseThreeComponentFluidSystem<Scalar> >
{
    using ThisType = ThreePhaseThreeComponentFluidSystem<Scalar>;
    using Base = Opm::BaseFluidSystem<Scalar, ThisType>;
    using PengRobinson =  typename Opm::PengRobinson<Scalar>;
    using PengRobinsonMixture =  typename Opm::PengRobinsonMixture<Scalar, ThisType>;

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
    static const int numPhases = 3;

    //! Index of the liquid phase
    static const int oilPhaseIdx = 0;
    static const int gasPhaseIdx = 1;
    static const int waterPhaseIdx = 2;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {"o",  // oleic phase
                                     "g",  // gas phase
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
    static const int Comp0Idx = 0;

    //! The component index of the solvent; co2
    static const int Comp1Idx = 1;

    //! The component index of the other main component; h2o + salts
    static const int Comp2Idx = 2;

    //! The component for pure oil
    using Comp0 = Opm::Octane<Scalar>;

    //! The component for pure solvent; since we have our own thermodynamic
    //! functions, we use the simple definition for the rest.
    using Comp1 = Opm::ChiwomsCO2<Scalar>;

    //! The component for brine
    using Comp2 = Opm::ChiwomsBrine<Scalar>;

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
            Comp0::name(),
            Comp1::name(),
            Comp2::name(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        return (compIdx == Comp0Idx)
                ? Comp0::molarMass()
                : (compIdx == Comp1Idx)
                  ? Comp1::molarMass()
                  : (compIdx == Comp2Idx)
                    ? Comp2::molarMass()
                    : throw std::invalid_argument("Molar mass component index");
    }

    /*!
         * \brief Critical temperature of a component [K].
         *
         * \copydetails Doxygen::compIdxParam
         */
    static Scalar criticalTemperature(unsigned compIdx)
    {
        return (compIdx == Comp0Idx)
                ? Comp0::criticalTemperature()
                : (compIdx == Comp1Idx)
                  ? Comp1::criticalTemperature()
                  : (compIdx == Comp2Idx)
                    ? Comp2::criticalTemperature()
                    : throw std::invalid_argument("Critical temperature component index");
    }

    /*!
         * \brief Critical pressure of a component [Pa].
         *
         * \copydetails Doxygen::compIdxParam
         */
    static Scalar criticalPressure(unsigned compIdx)
    {
        return (compIdx == Comp0Idx)
                ? Comp0::criticalPressure()
                : (compIdx == Comp1Idx)
                  ? Comp1::criticalPressure()
                  : (compIdx == Comp2Idx)
                    ? Comp2::criticalPressure()
                    : throw std::invalid_argument("Critical pressure component index");
    }

    /*!
         * \brief The acentric factor of a component [].
         *
         * \copydetails Doxygen::compIdxParam
         */
    static Scalar acentricFactor(unsigned compIdx)
    {
        return (compIdx == Comp0Idx)
                ? Comp0::acentricFactor()
                : (compIdx == Comp1Idx)
                  ? Comp1::acentricFactor()
                  : (compIdx == Comp2Idx)
                    ? Comp2::acentricFactor()
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
        const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, Comp1Idx));
        assert(T == (TEMPERATURE + 273.15));

        if(phaseIdx == oilPhaseIdx) {
            return 650.;
            return EOS::oleic_density(T, p, x);
        } else if(phaseIdx == gasPhaseIdx) {
            using IdealGas = Opm::IdealGas<Scalar>;
            return IdealGas::density(LhsEval(molarMass(Comp1Idx)), T, p);
        }
        else {
            return 1000.;
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
        const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, Comp1Idx));
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
        const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, Comp1Idx));

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

        if (phaseIdx == waterPhaseIdx || compIdx == Comp2Idx)
            return 1.0;
        else {
            Scalar phi = Opm::getValue(PengRobinsonMixture::computeFugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx));
            return phi;
        }
//         if (phaseIdx == oilPhaseIdx) {
// #if 0
// #warning HACK
//             if (compIdx == BrineIdx) {
//                 return 1e7/fluidState.pressure(oilPhaseIdx);

//             }else if (compIdx == OctaneIdx)
//                 return 40e3/fluidState.pressure(oilPhaseIdx);
//             else
//                 return 500e3/fluidState.pressure(oilPhaseIdx);
// #else
//             if (compIdx == BrineIdx) {
//                 return 1e7/fluidState.pressure(oilPhaseIdx);

//             }else if (compIdx == CO2Idx)
//                 return 500e3/fluidState.pressure(oilPhaseIdx);
//             else {
//                 return PengRobinsonMixture::computeFugacityCoefficient(fluidState,
//                                                                        paramCache,
//                                                                        phaseIdx,
//                                                                        compIdx);
//             }
// #endif
//         }
//         else if (phaseIdx == gasPhaseIdx) {
//             return 1.0;
//         }
//         else {
//             assert(phaseIdx == waterPhaseIdx);
//             if (compIdx == OctaneIdx) {
//                 return 10e6/fluidState.pressure(waterPhaseIdx);

//             }else if (compIdx == BrineIdx)
//                 return 80e3/fluidState.pressure(waterPhaseIdx);
//             else
//                 return 1e6/fluidState.pressure(waterPhaseIdx);
//         }

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
        // return 0;
        unsigned i = std::min(comp1Idx, comp2Idx);
        unsigned j = std::max(comp1Idx, comp2Idx);
#warning interactionCoefficients from Ivar
        if (i == Comp0Idx && j == Comp1Idx)
            return 0.1089;
        else if (i == Comp1Idx && j == Comp2Idx)
            return 1.1290;
        else if (i == Comp1Idx && j == Comp2Idx)
            return -0.0736;
        return 0;
    }

};

};//namespace opm

#endif // THREEPHASEFLUIDSYSTEM_HH
