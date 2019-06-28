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
 * \copydoc Opm::ChiFlash
 */
#ifndef OPM_CHI_FLASH_HPP
#define OPM_CHI_FLASH_HPP

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/Constants.hpp>

#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

#include <limits>
#include <iostream>

namespace Opm {

/*!
 * \brief Determines the phase compositions, pressures and saturations
 *        given the total mass of all components for the chiwoms problem.
 *
 */
template <class Scalar, class FluidSystem>
class ChiFlash
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { BrineIdx = FluidSystem::BrineIdx }; //rename for generic ?
    enum { OctaneIdx = FluidSystem::OctaneIdx }; //rename for generic ?
    enum { CO2Idx = FluidSystem::CO2Idx }; //rename for generic ?
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx};
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx};
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx};
    enum { numMiscibleComponents = 2}; //octane, co2
    enum { numMisciblePhases = 2}; //oil, gas
    enum {
        numEq =
        numMisciblePhases+
        numMisciblePhases*numMiscibleComponents
    };//pressure, saturation, composition

    enum {
        p0PvIdx = 0, // pressure first phase primary variable index
        S0PvIdx = 1, // saturation first phase primary variable index
        x00PvIdx = S0PvIdx + 1, // molefraction first phase first component primary variable index
        //numMiscibleComponennets*numMisciblePhases-1 molefractions/primvar follow
    };

public:
//    /*!
//     * \brief Guess initial values for all quantities.
//     */
//    template <class FluidState, class Evaluation = typename FluidState::Scalar>
//    static void guessInitial(FluidState& fluidState,
//                             const Dune::FieldVector<Evaluation, numComponents>& globalComposition)
//    {
//        // water saturation, no wate for now
//        Evaluation waterSaturation = 0.0;
//        fluidState.setSaturation(waterPhaseIdx, waterSaturation);

//        // oil and gas saturation
//        fluidState.setSaturation(oilPhaseIdx, (1.0-waterSaturation)/2.0);
//        fluidState.setSaturation(gasPhaseIdx, (1.0-waterSaturation)/2.0);

//        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
//            if (phaseIdx == waterPhaseIdx)
//                continue;
//            // composition
//            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx){
//                if (compIdx == BrineIdx){
//                    fluidState.setMoleFraction(phaseIdx,
//                                               compIdx,
//                                               0.0);
//                } else {
//                    fluidState.setMoleFraction(phaseIdx,
//                                               compIdx,
//                                               globalComposition[compIdx]);
//                }
//            }

//            // pressure. use atmospheric pressure as initial guess
//            fluidState.setPressure(phaseIdx, 20.0e6);
//        }

//        // set composition in water to only water
//        fluidState.setMoleFraction(waterPhaseIdx, BrineIdx, 1.0);
//        fluidState.setMoleFraction(waterPhaseIdx, OctaneIdx, 0.0);
//        fluidState.setMoleFraction(waterPhaseIdx, CO2Idx, 0.0);

//        // set the fugacity coefficients of all components in all phases
//        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
//        paramCache.updateAll(fluidState);
//        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
//            if (phaseIdx == waterPhaseIdx)
//                continue;
//            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
//                if (compIdx == BrineIdx)
//                    continue;
//                const typename FluidState::Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
//                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
//            }
//        }
//    }

    /*!
     * \brief Calculates the fluid state from the global mole fractions of the components and the phase pressures
     *
     */
    template <class FluidState>
    static void solve(FluidState& fluidState,
                      const Dune::FieldVector<typename FluidState::Scalar, numComponents>& globalComposition,
                      Scalar tolerance = -1.0)
    {
        typedef typename FluidState::Scalar InputEval;
        typedef Dune::FieldMatrix<InputEval, numEq, numEq> Matrix;
        typedef Dune::FieldVector<InputEval, numEq> Vector;
        typedef Opm::DenseAd::Evaluation</*Scalar=*/InputEval, /*numDerivs=*/numEq> FlashEval;
        typedef Dune::FieldVector<FlashEval, numEq> FlashDefectVector;
        typedef Opm::CompositionalFluidState<FlashEval, FluidSystem, /*energy=*/false> FlashFluidState;
        typedef Dune::FieldVector<typename FluidState::Scalar, numComponents> ComponentVector;

#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,7)
        Dune::FMatrixPrecision<InputEval>::set_singular_limit(1e-35);
#endif

        if (tolerance <= 0)
            tolerance = std::min<Scalar>(1e-3, 1e8*std::numeric_limits<Scalar>::epsilon());

        //initial guess for the K value, L value in RachfordRice
        ComponentVector K;
        Scalar L;
        Scalar L_min = 1e10;
        Scalar L_max = -1e10;
        Scalar totalMoles = 0;

        for (int compIdx = 0; compIdx<numComponents; ++compIdx) {
            if (compIdx == BrineIdx)
                continue;
            K[compIdx] = wilsonK_(fluidState, compIdx);

            L_min = Opm::min(L_min, 1/(1-K[compIdx]));
            L_max = Opm::max(L_max, 1/(1-K[compIdx]));
        }

        L = (L_min + L_max)/2;
        L = 0.5;
        std::cout << "global composition:" << globalComposition << std::endl;
        //Rachford Rice equation
        std::cout << "lambda:" << L << std::endl;
        L = solveRachfordRice_g_(K, L, globalComposition);
        std::cout << "lambda ny:" << L << std::endl;
        //Phase stability test
        bool isStable;
        ComponentVector x;
        ComponentVector y;
        phaseStabilityTest_(isStable, x, y, fluidState, globalComposition);

        //update the composition using newton
        newtonCompositionUpdate_(K, L, fluidState, globalComposition);

        // compressibility
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(fluidState, oilPhaseIdx);
        paramCache.updatePhase(fluidState, gasPhaseIdx);

        //compressibility
        const Scalar R = Opm::Constants<Scalar>::R;
        Scalar Z_L = (paramCache.molarVolume(oilPhaseIdx) * fluidState.pressure(oilPhaseIdx) )/
                (R * fluidState.temperature(oilPhaseIdx));
        Scalar Z_V = (paramCache.molarVolume(gasPhaseIdx) * fluidState.pressure(gasPhaseIdx) )/
                (R * fluidState.temperature(gasPhaseIdx));

        Scalar Sw = 0.0; //todo: include water from conservation eq
        Scalar So = L*Z_L/(L*Z_L+(1-L)*Z_V);
        Scalar Sg = 1-So-Sw;

        fluidState.setSaturation(waterPhaseIdx, Sw);
        fluidState.setSaturation(oilPhaseIdx, So);
        fluidState.setSaturation(gasPhaseIdx, Sg);
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * This is a convenience method which assumes that the capillary pressure is
     * zero...
     */
    template <class FluidState, class ComponentVector>
    static void solve(FluidState& fluidState,
                      const ComponentVector& globalMolarities,
                      Scalar tolerance = 0.0)
    {
        typedef NullMaterialTraits<Scalar, numPhases> MaterialTraits;
        typedef NullMaterial<MaterialTraits> MaterialLaw;
        typedef typename MaterialLaw::Params MaterialLawParams;

        MaterialLawParams matParams;
        solve<MaterialLaw>(fluidState, matParams, globalMolarities, tolerance);
    }


protected:

//    template <class FluidState>
//    static void printFluidState_(const FluidState& fluidState)
//    {
//        typedef typename FluidState::Scalar FsScalar;

//        std::cout << "saturations: ";
//        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
//            std::cout << fluidState.saturation(phaseIdx) << " ";
//        std::cout << "\n";

//        std::cout << "pressures: ";
//        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
//            std::cout << fluidState.pressure(phaseIdx) << " ";
//        std::cout << "\n";

//        std::cout << "densities: ";
//        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
//            std::cout << fluidState.density(phaseIdx) << " ";
//        std::cout << "\n";

//        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
//            std::cout << "composition " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
//            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
//                std::cout << fluidState.moleFraction(phaseIdx, compIdx) << " ";
//            }
//            std::cout << "\n";
//        }

//        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
//            std::cout << "fugacities " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
//            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
//                std::cout << fluidState.fugacity(phaseIdx, compIdx) << " ";
//            }
//            std::cout << "\n";
//        }

//        std::cout << "global component molarities: ";
//        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
//            FsScalar sum = 0;
//            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
//                sum += fluidState.saturation(phaseIdx)*fluidState.molarity(phaseIdx, compIdx);
//            }
//            std::cout << sum << " ";
//        }
//        std::cout << "\n";
//    }

    template <class FlashFluidState>
    static typename FlashFluidState::Scalar wilsonK_(const FlashFluidState& fluidState, int compIdx)
    {
        typedef typename FlashFluidState::Scalar FlashEval;
        const auto& acf = FluidSystem::acentricFactor(compIdx);
        const auto& T_crit = FluidSystem::criticalTemperature(compIdx);
        const auto& T = fluidState.temperature(0);
        const auto& p_crit = FluidSystem::criticalPressure(compIdx);
        const auto& p = fluidState.pressure(0); //for now assume no capillary pressure

        const auto& tmp = Opm::exp(5.3727 * (1+acf) * (1-T_crit/T)) * (p_crit/p);
        return tmp;
    }

    template <class Vector>
    static typename Vector::field_type rachfordRice_g_(const Vector& K, const Scalar L, const Vector& globalComposition)
    {
        typename Vector::field_type g=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            g += (globalComposition[compIdx]*(K[compIdx]-1))/(1+L*(K[compIdx]-1));
        }
        return g;
    }

    template <class Vector>
    static typename Vector::field_type rachfordRice_dg_dL_(const Vector& K, const Scalar L, const Vector& globalComposition)
    {
        typename Vector::field_type dg=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            dg += -(globalComposition[compIdx]*(K[compIdx]-1)*(K[compIdx]-1))/((1+L*(K[compIdx]-1))*(1+L*(K[compIdx]-1)));
        }
        return dg;
    }

    template <class Vector>
    static typename Vector::field_type solveRachfordRice_g_(const Vector& K, Scalar L, const Vector& globalComposition)
    {
        for (int iteration=0; iteration<200; ++iteration){
            //L=Lold+g/dg;
            Scalar g = rachfordRice_g_(K, L, globalComposition);
            Scalar dg_dL = rachfordRice_dg_dL_(K, L, globalComposition);
            L -= g/dg_dL;
            //check for convergence
            Scalar delta = g/dg_dL;
            if ( Opm::abs(delta) < 1e-10 )
                return L;
        }
        throw std::runtime_error("Rachford Rice did not converge within maximum number of iterations" );

    }

    template <class FlashFluidState, class ComponentVector>
    static void phaseStabilityTest_(bool& isStable, ComponentVector& x, ComponentVector& y, const FlashFluidState& fluidState, const ComponentVector& globalComposition)
    {
        bool isTrivialL, isTrivialV;
        ComponentVector K_l, K_v;
        Scalar S_l, S_v;

        checkStability_(fluidState, isTrivialV, K_v, y, S_v, globalComposition, /*isGas=*/true);
        bool V_unstable = (S_v < (1.0 + 1e-5)) || isTrivialV;

        checkStability_(fluidState, isTrivialL, K_l, x, S_l, globalComposition, /*isGas=*/false);
        bool L_stable = (S_l < (1.0 + 1e-5)) || isTrivialL;

        isStable = L_stable && V_unstable; //L-stable means succes in making liquid, V-unstable means no success in making vapour
        if (isStable) {
            // single phase, i.e. phase composition is equivalent to the global composition
            x = globalComposition;
            y = globalComposition;
        }
    }

    template <class FlashFluidState, class ComponentVector>
    static void checkStability_(const FlashFluidState& fluidState, bool& isTrivial, ComponentVector& K, ComponentVector& xy_loc, Scalar& S_loc, const ComponentVector& globalComposition, bool isGas)
    {
        typedef typename FlashFluidState::Scalar FlashEval;
        //this is "Michelsens stability test"
        //make two fake phases "inside" one phase and check for positive volume
        FlashFluidState fluidState_fake = fluidState;
        FlashFluidState fluidState_global = fluidState;

        for (int compIdx = 0; compIdx < numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            K[compIdx] = wilsonK_(fluidState, compIdx);
        }
        for (int i = 0; i < 19000; ++i) {
            S_loc = 0.0;
            if (isGas) {
                xy_loc;
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    xy_loc[compIdx] = K[compIdx] * globalComposition[compIdx];
                    S_loc += xy_loc[compIdx];
                }
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    xy_loc[compIdx] /= S_loc;
                    fluidState_fake.setMoleFraction(gasPhaseIdx, compIdx, xy_loc[compIdx]);
                }
            }
            else {
                ComponentVector xy_loc;
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    if (compIdx == BrineIdx)
                        continue;
                    xy_loc[compIdx] = globalComposition[compIdx]/K[compIdx];
                    S_loc += xy_loc[compIdx];
                }
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    if (compIdx == BrineIdx)
                        continue;
                    xy_loc[compIdx] /= S_loc;
                    fluidState_fake.setMoleFraction(oilPhaseIdx, compIdx, xy_loc[compIdx]);
                }
            }

            int phaseIdx = (isGas?gasPhaseIdx:oilPhaseIdx);
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                fluidState_global.setMoleFraction(phaseIdx, compIdx, globalComposition[compIdx]);
            }

            typename FluidSystem::template ParameterCache<FlashEval> paramCache_fake;
            paramCache_fake.updatePhase(fluidState_fake, phaseIdx);

            typename FluidSystem::template ParameterCache<FlashEval> paramCache_global;
            paramCache_global.updatePhase(fluidState_global, phaseIdx);

            //fugacity for fake phases each component
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                Scalar phiFake = FluidSystem::fugacityCoefficient(fluidState_fake, paramCache_fake, phaseIdx, compIdx);
                Scalar phiGlobal = FluidSystem::fugacityCoefficient(fluidState_global, paramCache_global, phaseIdx, compIdx);

                fluidState_fake.setFugacityCoefficient(phaseIdx, compIdx, phiFake);
                fluidState_global.setFugacityCoefficient(phaseIdx, compIdx, phiGlobal);
            }

            ComponentVector R;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                if (isGas)
                    R[compIdx] = fluidState_global.fugacity(oilPhaseIdx, compIdx)/fluidState_fake.fugacity(oilPhaseIdx, compIdx)/S_loc;
                else
                    R[compIdx] = fluidState_fake.fugacity(gasPhaseIdx, compIdx)/fluidState_global.fugacity(gasPhaseIdx, compIdx)*S_loc;
            }

            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                K[compIdx] *= R[compIdx];
            }
            Scalar R_norm = 0.0;
            Scalar K_norm = 0.0;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                auto a = R[compIdx] - 1.0;
                auto b = Opm::log(K[compIdx]);
                R_norm += a*a;
                K_norm += b*b;
            }
            isTrivial = (K_norm < 1e-5);
            if (isTrivial || R_norm < 1e-10)
                return;
            //todo: make sure that no molefraction is smaller than 1e-8 ?
            //todo: take care of water!
        }
        throw std::runtime_error("stability test did not converge");
    }//end checkStability

    template <class FlashFluidState, class ComponentVector>
    static void newtonCompositionUpdate_(ComponentVector& K, Scalar& L, FlashFluidState& fluidState, const ComponentVector& globalComposition)
    {
        if(L == 0 || L == 1) {
            //single-phase gas or single-phase oil
            for(int compIdx; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                fluidState.setMoleFraction(oilPhaseIdx, compIdx, globalComposition[compIdx]);
                fluidState.setMoleFraction(gasPhaseIdx, compIdx, globalComposition[compIdx]);
            }
            return;
        }

        //two-phase case
        ComponentVector x;
        ComponentVector y;
        //normalization
        Scalar sumx=0;
        Scalar sumy=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            x[compIdx] = globalComposition[compIdx]/(L + (1-L)*K[compIdx]);
            sumx += x[compIdx];
            y[compIdx] = (K[compIdx]*globalComposition[compIdx])/(L + (1-L)*K[compIdx]);
            sumy += y[compIdx];
        }
        x /= sumx;
        y /= sumy;

        for(int compIdx; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            fluidState.setMoleFraction(oilPhaseIdx, compIdx, x[compIdx]);
            fluidState.setMoleFraction(gasPhaseIdx, compIdx, y[compIdx]);
        }

        //newton
        typedef Dune::FieldVector<Scalar, numMiscibleComponents*2+1> NewtonVector;
        typedef Dune::FieldMatrix<Scalar, numMiscibleComponents*2+1, numMiscibleComponents*2+1> NewtonMatrix;
        NewtonVector newtonX;
        NewtonVector newtonB;
        NewtonMatrix newtonA;
        NewtonVector newtonDelta;

        newtonX[0] = fluidState.moleFraction(oilPhaseIdx, OctaneIdx);
        newtonX[1] = fluidState.moleFraction(oilPhaseIdx, CO2Idx);
        newtonX[2] = fluidState.moleFraction(gasPhaseIdx, OctaneIdx);
        newtonX[3] = fluidState.moleFraction(gasPhaseIdx, CO2Idx);
        newtonX[4] = L;

        for (int i = 0; i< 100; ++i){
            evalDefect_(newtonB, newtonX, fluidState, globalComposition);
            evalJacobian_(newtonA, newtonX, fluidState, globalComposition);
            newtonA.solve(newtonDelta, newtonB);
            newtonX -= newtonDelta;
            if(std::abs(newtonDelta.one_norm())<1e-6)
                break;
        }
        throw std::runtime_error("Newton composition update did not converge within maxIterations");
    }

    template <class FluidState, class DefectVector, class ComponentVector>
    static void evalDefect_(DefectVector& b,
                            DefectVector& x,
                            const FluidState& fluidStateIn,
                            const ComponentVector& globalComposition)
    {
        FluidState fluidState(fluidStateIn);
        //primary variables (numMisciblePhases*numMiscibleComponents +1 )
        fluidState.setMoleFraction(oilPhaseIdx, OctaneIdx, x[0]);
        fluidState.setMoleFraction(oilPhaseIdx, CO2Idx, x[1]);
        fluidState.setMoleFraction(gasPhaseIdx, OctaneIdx, x[2]);
        fluidState.setMoleFraction(gasPhaseIdx, CO2Idx, x[3]);
        Scalar L = x[4];
        // todo: remove this hardcode for water
        fluidState.setMoleFraction(oilPhaseIdx, BrineIdx, 0);
        fluidState.setMoleFraction(gasPhaseIdx, BrineIdx, 0);

        //compute fugacities
        typedef typename FluidSystem::template ParameterCache<typename FluidState::Scalar> ParamCache;
        ParamCache paramCache;
        for (int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
            if (phaseIdx==waterPhaseIdx)
                continue;
            paramCache.updatePhase(fluidState, phaseIdx);
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }

        //todo: make this AD (x, y and L)
        // numMisciblePhases*numMiscible components +1 primary nknowns:L, x(comp1), x(comp2), y(comp1), y(comp2)
        // numMisciblePhases**numMiscible components +1 equations:
        // eq1: z-Lx-(1-L)y=0 (comp1)
        b[0] = globalComposition[OctaneIdx] - L*fluidState.moleFraction(oilPhaseIdx, OctaneIdx)
                - (1-L)*fluidState.moleFraction(gasPhaseIdx, OctaneIdx);
        // eq1: z-Lx-(1-L)y=0 (comp2)
        b[1] = globalComposition[CO2Idx] - L*fluidState.moleFraction(oilPhaseIdx, CO2Idx)
                - (1-L)*fluidState.moleFraction(gasPhaseIdx, CO2Idx);
        // eq3: f_Gas-f_Oil (comp1)
        b[2] = fluidState.fugacity(oilPhaseIdx, OctaneIdx) - fluidState.fugacity(gasPhaseIdx, OctaneIdx);
        // eq4: f_Gas-f_Oil (comp1)
        b[3] = fluidState.fugacity(oilPhaseIdx, CO2Idx) - fluidState.fugacity(gasPhaseIdx, CO2Idx);
        // eq4: sum(x)=sum(y)
        b[4] = fluidState.moleFraction(oilPhaseIdx, OctaneIdx) - fluidState.moleFraction(gasPhaseIdx, OctaneIdx);
        b[4] += fluidState.moleFraction(oilPhaseIdx, CO2Idx) - fluidState.moleFraction(gasPhaseIdx, CO2Idx);
    }//end valDefect

    template <class FluidState, class DefectVector, class DefectMatrix, class ComponentVector>
    static void evalJacobian_(DefectMatrix& A,
                              const DefectVector& xIn,
                              const FluidState& fluidStateIn,
                              const ComponentVector& globalComposition)
    {
        DefectVector x(xIn);
        DefectVector b0;
        evalDefect_(b0, x, fluidStateIn, globalComposition);
        Scalar epsilon = 1e-8;
        //make the jacobian of Ax=b
        for(int i=0; i<b0.size(); ++i){
            x[i] += epsilon;
            DefectVector bEps;
            evalDefect_(bEps, x, fluidStateIn, globalComposition);
            x[i] -= epsilon;
            //derivative of all eqs wrt primary variable i
            DefectVector derivI(bEps);
            derivI -= b0;
            derivI /= epsilon;

            for(int j=0; j<b0.size(); ++j){
                A[j][i] = derivI[j];
            }
        }
    }//end evalJacobian
};//end ChiFlash

} // namespace Opm

#endif
