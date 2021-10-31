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
#include <opm/material/eos/PengRobinsonMixture.hpp>

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

        // PSEUDOCODE for what the code does:
        // for two-phase --> go directly to flash (TODO)
        // TODO: if two-phase, go directly to flash -> IMPLEMENT ??
        // "cells that are two-phase can be flashed directly"
        // EDIT: cells that are two-phase (skip below)
        // twoPhase = getTwoPhaseFlag(state)
        // if single-phase
            // check if it stays single phase (P, T, globalComp, K)
            // --> x0, y0, updatedSingle
            // K0(updatedSingle) = y0(updatedSingle)./x0(updatedSingle);
            // solveRachfordRice (L0, Ko, z);
            // --> new L0
            // calculate Si_L, Si_V, A_L, A_V, B_L, B_B, Bi
            // calculate ZO_L, ZO_V
        // end if single-phase

            // for "active" cells
            // newtonCompositionUpdate(P, T, z, K, L);
            // --> x, y, K, Z_L, Z_V, L


        // for single-phase : DO STABILITYTEST --> initial K0 (TODO)
        // solve Rachford Rice (using K0, L0, z) for stabilityTest showing that singlePhase --> twoPhase --> L0 (TODO)
        // newtonCompositionUpdate using p, T, z, K and L

        // Initial guess for the K value using Wilson's formula
        std::cout << " ======== WILSON ========" << std::endl;
        ComponentVector K;
        for (int compIdx = 0; compIdx<numComponents; ++compIdx) {
            if (compIdx == BrineIdx)
                continue;
            K[compIdx] = wilsonK_(fluidState, compIdx);
        }

        std::cout << "global composition : " << globalComposition << std::endl;

        std::cout << " +++++++++++++++++++++++++++++" << std::endl;
        std::cout << " ++++++++ START FLASH ++++++++" << std::endl;
        std::cout << " +++++++++++++++++++++++++++++" << std::endl;
        
        // Some declarations
        bool isStable = false;
        Scalar L;
        // Flash loop start
        for (int i = 0; i< 100; ++i) {
            
            std::cout << "Iteration : " << i << std::endl;

            // In the first iteration we do stability test to check if cell is single-phase. If so, we do not need to 
            // do Newton update
            if (i == 0) {
                // Phase stability test
                std::cout << " ======== STABILITY ========" << std::endl;
                phaseStabilityTest_(isStable, K, fluidState, globalComposition);

                //Rachford Rice equation
                std::cout << " ======== RACHFORD-RICE ========" << std::endl;
                L = solveRachfordRice_g_(K, globalComposition);
                std::cout << "L :" << L << std::endl;
            } // end if i == 0
            
            // Update the composition using Newton's method if cell is two-phase
            if (isStable == false) {
                std::cout << " ======== NEWTON ========" << std::endl;
                newtonCompositionUpdate_(K, L, fluidState, globalComposition);
            }  // end if newtonComposition update
            else break; // break out of flash loop
        } // end flash loop
        std::cout << " +++++++++++++++++++++++++++" << std::endl;
        std::cout << " ++++++++ END FLASH ++++++++" << std::endl;
        std::cout << " +++++++++++++++++++++++++++" << std::endl;

        // Update phases
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(fluidState, oilPhaseIdx);
        paramCache.updatePhase(fluidState, gasPhaseIdx);

        // Calculate compressibility factor
        const Scalar R = Opm::Constants<Scalar>::R;
        Scalar Z_L = (paramCache.molarVolume(oilPhaseIdx) * fluidState.pressure(oilPhaseIdx) )/
                (R * fluidState.temperature(oilPhaseIdx));
        Scalar Z_V = (paramCache.molarVolume(gasPhaseIdx) * fluidState.pressure(gasPhaseIdx) )/
                (R * fluidState.temperature(gasPhaseIdx));

        // Update saturation
        Scalar Sw = 0.0; //todo: include water from conservation eq
        Scalar So = L*Z_L/(L*Z_L+(1-L)*Z_V);
        Scalar Sg = 1-So-Sw;

        fluidState.setSaturation(waterPhaseIdx, Sw);
        fluidState.setSaturation(oilPhaseIdx, So);
        fluidState.setSaturation(gasPhaseIdx, Sg);

        // Update densities
        fluidState.setDensity(oilPhaseIdx, FluidSystem::density(fluidState, paramCache, oilPhaseIdx));
        fluidState.setDensity(gasPhaseIdx, FluidSystem::density(fluidState, paramCache, gasPhaseIdx));
        fluidState.setDensity(waterPhaseIdx, FluidSystem::density(fluidState, paramCache, waterPhaseIdx));

        // Update water mole fraction
        fluidState.setMoleFraction(waterPhaseIdx, BrineIdx, 1.0);
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
            g += (globalComposition[compIdx]*(K[compIdx]-1))/(K[compIdx]-L*(K[compIdx]-1));
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
            dg += (globalComposition[compIdx]*(K[compIdx]-1)*(K[compIdx]-1))/((K[compIdx]-L*(K[compIdx]-1))*(K[compIdx]-L*(K[compIdx]-1)));
        }
        return dg;
    }

    template <class Vector>
    static typename Vector::field_type solveRachfordRice_g_(const Vector& K, const Vector& globalComposition)
    {
        // Find min and max K. Have to do a laborious for loop to avoid water component (where K=0)
        // TODO: Replace loop with Dune::min_value() and Dune::max_value() when water component is properly handled
        Scalar Kmin = K[0];
        Scalar Kmax = K[0];
        for (int compIdx=1; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            if (K[compIdx] < Kmin)
                Kmin = K[compIdx];
            else if (K[compIdx] >= Kmax)
                Kmax = K[compIdx];
        }
        // Lower and upper bound for solution
        Scalar Lmin = (Kmin / (Kmin - 1));
        Scalar Lmax = Kmax / (Kmax - 1);

        // Check if Lmin < Lmax, and switch if not
        if (Lmin > Lmax)
        {
            Scalar Ltmp = Lmin;
            Lmin = Lmax;
            Lmax = Ltmp;
        }

        // Initial guess
        Scalar L = (Lmin + Lmax)/2;

        // Newton-Rahpson loop
        for (int iteration=0; iteration<200; ++iteration){
            // Calculate function and derivative values
            Scalar g = rachfordRice_g_(K, L, globalComposition);
            Scalar dg_dL = rachfordRice_dg_dL_(K, L, globalComposition);

            // Lnew = Lold - g/dg;
            Scalar delta = g/dg_dL;
            L -= delta;

            // Check if L is within the bounds, and if not, we apply bisection method
            if (L < Lmin || L > Lmax)
                {
                    L = bisection_g_(K, Lmin, Lmax, globalComposition);
                    return L;
                }

            // Check for convergence
            if ( Opm::abs(delta) < 1e-10 ) {
                L = Opm::min(Opm::max(L, 0), 1);
                return L;
            }
        }
        // Throw error if Rachford-Rice fails
        throw std::runtime_error("Rachford-Rice did not converge within maximum number of iterations" );
    }

    template <class Vector>
    static typename Vector::field_type bisection_g_(const Vector& K, Scalar Lmin, Scalar Lmax, const Vector& globalComposition)
    {
        // Check if g(Lmin) and g(Lmax) have opposite sign
        Scalar gLmin = rachfordRice_g_(K, Lmin, globalComposition);
        Scalar gLmax = rachfordRice_g_(K, Lmax, globalComposition);
        if (Dune::sign(gLmin) == Dune::sign(gLmax))
        {
            throw std::runtime_error("Lmin and Lmax are incorrect for bisection");
        }
            
        // Bisection loop
        for (int iteration=0; iteration<200; ++iteration){
            // New midpoint
            Scalar L = (Lmin + Lmax) / 2;
            Scalar gMid = rachfordRice_g_(K, L, globalComposition);

            // Check if midpoint fulfills g=0 or L - Lmin is sufficiently small
            if (Opm::abs(gMid) < 1e-10 || Opm::abs((Lmax - Lmin) / 2) < 1e-10)
                return L;
            
            // Else we repeat with midpoint being either Lmin og Lmax (depending on the signs)
            else if (Dune::sign(gMid) != Dune::sign(gLmin))
                Lmax = L; 
            else
                Lmin = L;
        }
        throw std::runtime_error("Rachford-Rice with bisection failed!");
    }

    template <class FlashFluidState, class ComponentVector>
    static void phaseStabilityTest_(bool& isStable, ComponentVector& K, FlashFluidState& fluidState, const ComponentVector& globalComposition)
    {
        bool isTrivialL, isTrivialV;
        ComponentVector K_l, K_v, x, y;
        Scalar S_l, S_v;

        // Use equilibrium constants made with Wilson's formula (before flash loop)
        K_l = K;
        K_v = K;

        // Check for vapour instable phase
        checkStability_(fluidState, isTrivialV, K_v, y, S_v, globalComposition, /*isGas=*/true);
        bool V_unstable = (S_v < (1.0 + 1e-5)) || isTrivialV;

        // Check for liquids stable phase
        checkStability_(fluidState, isTrivialL, K_l, x, S_l, globalComposition, /*isGas=*/false);
        bool L_stable = (S_l < (1.0 + 1e-5)) || isTrivialL;

        // L-stable means success in making liquid, V-unstable means no success in making vapour
        isStable = L_stable && V_unstable; 
        if (isStable) {
            // Single phase, i.e. phase composition is equivalent to the global composition
            // Update fluidstate with mole fration
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                fluidState.setMoleFraction(gasPhaseIdx, compIdx, globalComposition[compIdx]);
                fluidState.setMoleFraction(oilPhaseIdx, compIdx, globalComposition[compIdx]);
            }
        }
        // If it is not stable, we use the mole fractions from the Michelsen test to calculate a new K
        else {
            for (int compIdx = 0; compIdx<numComponents; ++compIdx) {
                if (compIdx == BrineIdx)
                    continue;
                K[compIdx] = y[compIdx] / x[compIdx];
            }
        }
    }

    template <class FlashFluidState, class ComponentVector>
    static void checkStability_(const FlashFluidState& fluidState, bool& isTrivial, ComponentVector& K, ComponentVector& xy_loc, Scalar& S_loc, const ComponentVector& globalComposition, bool isGas)
    {
        typedef typename FlashFluidState::Scalar FlashEval;
        typedef ThreePhaseCo2OctaneBrineFluidSystem<Scalar> ThisType;
        typedef typename Opm::PengRobinsonMixture<Scalar, ThisType> PengRobinsonMixture;

        // Declarations 
        FlashFluidState fluidState_fake = fluidState;
        FlashFluidState fluidState_global = fluidState;

        // Michelsens stability test.
        // Make two fake phases "inside" one phase and check for positive volume
        for (int i = 0; i < 19000; ++i) {
            S_loc = 0.0;
            if (isGas) {
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
                Scalar phiFake = PengRobinsonMixture::computeFugacityCoefficient(fluidState_fake, paramCache_fake, phaseIdx, compIdx);
                Scalar phiGlobal = PengRobinsonMixture::computeFugacityCoefficient(fluidState_global, paramCache_global, phaseIdx, compIdx);

                fluidState_fake.setFugacityCoefficient(phaseIdx, compIdx, phiFake);
                fluidState_global.setFugacityCoefficient(phaseIdx, compIdx, phiGlobal);
            }

           
            ComponentVector R;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                if (isGas){
                    Scalar fug_fake = fluidState_fake.fugacity(gasPhaseIdx, compIdx);
                    Scalar fug_global = fluidState_global.fugacity(gasPhaseIdx, compIdx);
                    Scalar fug_ratio = fug_global / fug_fake;
                    R[compIdx] = fug_ratio / S_loc;
                }
                else{
                    Scalar fug_fake = fluidState_fake.fugacity(oilPhaseIdx, compIdx);
                    Scalar fug_global = fluidState_global.fugacity(oilPhaseIdx, compIdx);
                    Scalar fug_ratio = fug_fake / fug_global;
                    R[compIdx] = fug_ratio * S_loc;
                }
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
        if (L == 0 || L == 1) {
            //single-phase gas or single-phase oil
            for(int compIdx; compIdx<numComponents; ++compIdx){
                if (compIdx == BrineIdx)
                    continue;
                fluidState.setMoleFraction(oilPhaseIdx, compIdx, globalComposition[compIdx]);
                fluidState.setMoleFraction(gasPhaseIdx, compIdx, globalComposition[compIdx]);
            }
            return;
        }

        // Calculate x and y, and normalize
        ComponentVector x;
        ComponentVector y;
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

        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            fluidState.setMoleFraction(oilPhaseIdx, compIdx, x[compIdx]);
            fluidState.setMoleFraction(gasPhaseIdx, compIdx, y[compIdx]);
        }

        // Newton
        typedef Dune::FieldVector<Scalar, numMiscibleComponents*numMisciblePhases+1> NewtonVector;
        typedef Dune::FieldMatrix<Scalar, numMiscibleComponents*numMisciblePhases+1, numMiscibleComponents*numMisciblePhases+1> NewtonMatrix;
        NewtonVector newtonX;
        NewtonVector newtonB;
        NewtonMatrix newtonA;
        NewtonVector newtonDelta;

        // Assign primary variables (mole fractions + L)
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            newtonX[compIdx] = fluidState.moleFraction(oilPhaseIdx, compIdx);
            newtonX[compIdx + numMiscibleComponents] = fluidState.moleFraction(gasPhaseIdx, compIdx);
        }
        newtonX[numMisciblePhases*numMiscibleComponents] = L;

        // Main Newton loop
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
        // Set primary variables in FluidState for misc. calculations
        // Primary variables = numMisciblePhases*numMiscibleComponents + 1
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            fluidState.setMoleFraction(oilPhaseIdx, compIdx, x[compIdx]);
            fluidState.setMoleFraction(gasPhaseIdx, compIdx, x[compIdx + numMiscibleComponents]);
        }
        Scalar L = x[numMiscibleComponents*numMisciblePhases];
        fluidState.setMoleFraction(oilPhaseIdx, BrineIdx, 0); /* OBS */
        fluidState.setMoleFraction(gasPhaseIdx, BrineIdx, 0); /* OBS */

        // Compute fugacities
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

        // Compute residuals for Newton update:
        // Primary variables are: x, y and L
        // TODO: Make this AD
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            if (compIdx == BrineIdx)
                continue;
            // z - Lx - (1-L)y = 0
            b[compIdx] = globalComposition[compIdx] - L*fluidState.moleFraction(oilPhaseIdx, compIdx) - (1-L)* fluidState.moleFraction(gasPhaseIdx,compIdx);
            // f_vapor - f_liquid = 0
            b[compIdx + numMiscibleComponents] = fluidState.fugacity(oilPhaseIdx, compIdx) - fluidState.fugacity(gasPhaseIdx, compIdx);
            // sum(x) - sum(y) = 0
            b[numMisciblePhases*numMiscibleComponents] += fluidState.moleFraction(oilPhaseIdx, compIdx) - fluidState.moleFraction(gasPhaseIdx,compIdx);
        }
    }//end valDefect

    template <class FluidState, class DefectVector, class DefectMatrix, class ComponentVector>
    static void evalJacobian_(DefectMatrix& A,
                              const DefectVector& xIn,
                              const FluidState& fluidStateIn,
                              const ComponentVector& globalComposition)
    {
        // TODO: Use AD instead
        // Calculate response of current state x
        DefectVector x(xIn);
        DefectVector b0;
        evalDefect_(b0, x, fluidStateIn, globalComposition);

        // Make the jacobian A in Newton system Ax=b
        Scalar epsilon = 1e-8;
        for(int i=0; i<b0.size(); ++i){
            // Permutate x and calculate response
            x[i] += epsilon;
            DefectVector bEps;
            evalDefect_(bEps, x, fluidStateIn, globalComposition);
            x[i] -= epsilon;

            // Forward difference of all eqs wrt primary variable i
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
