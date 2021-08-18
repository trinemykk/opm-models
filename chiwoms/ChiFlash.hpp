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
#include <iomanip>

namespace Opm {

/*!
 * \brief Determines the phase compositions, pressures and saturations
 *        given the total mass of all components for the chiwoms problem.
 *
 */
template <class Scalar, class Evaluation, class FluidSystem>
class ChiFlash
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { Comp2Idx = FluidSystem::Comp2Idx }; //rename for generic ?
    enum { Comp0Idx = FluidSystem::Comp0Idx }; //rename for generic ?
    enum { Comp1Idx = FluidSystem::Comp1Idx }; //rename for generic ?
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
    /*!
     * \brief Calculates the fluid state from the global mole fractions of the components and the phase pressures
     *
     */
    template <class FluidState>
    static void solve(FluidState& fluidState,
                      const Dune::FieldVector<typename FluidState::Scalar, numComponents>& globalComposition,
                      int verbosity,
                      std::string twoPhaseMethod,
                      Scalar tolerance)
    {

        using InputEval = typename FluidState::Scalar;
        using Matrix = Dune::FieldMatrix<InputEval, numEq, numEq>;
        using Vetor = Dune::FieldVector<InputEval, numEq>;
        using FlashEval = Opm::DenseAd::Evaluation</*Scalar=*/InputEval, /*numDerivs=*/numEq>;
        using FlashDefectVector = Dune::FieldVector<FlashEval, numEq>;
        using FlashFluidState = Opm::CompositionalFluidState<FlashEval, FluidSystem, /*energy=*/false>;
        using ComponentVector = Dune::FieldVector<typename FluidState::Scalar, numComponents>;

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
        ComponentVector K;
        for (int compIdx = 0; compIdx<numComponents; ++compIdx) {
            K[compIdx] = wilsonK_(fluidState, compIdx);
        }

        // Some declarations
        bool isStable = false;
        InputEval L;
        // First we do stability test to check if cell is single-phase. If so, we do not need to do Newton update
        // Phase stability test
        phaseStabilityTest_(isStable, K, fluidState, globalComposition);
        
        // Rachford Rice equation
        L = solveRachfordRice_g_(K, globalComposition);
        
        // Update the composition using Newton's method if cell is two-phase
        if (isStable == false) {
            if (twoPhaseMethod == "newton"){
                if (verbosity == 1) {
                    std::cout << "********" << std::endl;
                    std::cout << "Cell is two-phase; solve using Newton!" << std::endl;
                    std::cout << "Initial variables: K = [" << K[Comp0Idx] << " " << K[Comp1Idx] << "], z = [" << globalComposition[Comp0Idx] << " " << globalComposition[Comp1Idx] <<"], and L = " << L << std::endl;
                }
                newtonCompositionUpdate_(K, L, fluidState, globalComposition, verbosity);
                if (verbosity == 1) {
                    std::cout << "********" << std::endl;
                }
            }
            else if (twoPhaseMethod == "ssi"){
                if (verbosity == 1) {
                    std::cout << "********" << std::endl;
                    std::cout << "Cell is two-phase; solve using Succcessive Substitution!" << std::endl;
                    std::cout << "Initial variables: K = [" << K[Comp0Idx] << " " << K[Comp1Idx] << "], z = [" << globalComposition[Comp0Idx] << " " << globalComposition[Comp1Idx] <<"], and L = " << L << std::endl;
                }
                successiveSubstitutionComposition_(K, L, fluidState, globalComposition, /*standAlone=*/true, verbosity);
                if (verbosity == 1) {
                    std::cout << "********" << std::endl;
                }
            }
        }
        else{
            // Cell is one-phase. Make sure L is either 1 or 0.
            if (L > 0.5)
                L = 1.0;
                if (verbosity >= 2) {std::cout << "Cell is pure liquid (L == 1.0)." << std::endl;}
            else if (L <= 0.5)
                L = 0.0;
                if (verbosity >= 2) {std::cout << "Cell is pure vapor (L == 0.0)." << std::endl;}
        }

        // Update phases
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        paramCache.updatePhase(fluidState, oilPhaseIdx);
        paramCache.updatePhase(fluidState, gasPhaseIdx);

        // Calculate compressibility factor
        const Scalar R = Opm::Constants<Scalar>::R;
        Evaluation Z_L = (paramCache.molarVolume(oilPhaseIdx) * fluidState.pressure(oilPhaseIdx) )/
                (R * fluidState.temperature(oilPhaseIdx));
        Evaluation Z_V = (paramCache.molarVolume(gasPhaseIdx) * fluidState.pressure(gasPhaseIdx) )/
                (R * fluidState.temperature(gasPhaseIdx));

        // Update saturation
        Evaluation So = (L*Z_L/(L*Z_L+(1-L)*Z_V));
        //Evaluation Sg = 1-So;
        Evaluation Sg = (1-L)*Z_V/(L*Z_L+(1-L)*Z_V);
        fluidState.setSaturation(oilPhaseIdx, So);
        fluidState.setSaturation(gasPhaseIdx, Sg);

        // Print info
        if (verbosity >= 2) {
            std::cout << "So = " << So <<std::endl;
            std::cout << "Sg = " << Sg <<std::endl;
            std::cout << "Sw = " << Sw <<std::endl;
        }

        // Update densities
        fluidState.setDensity(oilPhaseIdx, FluidSystem::density(fluidState, paramCache, oilPhaseIdx));
        fluidState.setDensity(gasPhaseIdx, FluidSystem::density(fluidState, paramCache, gasPhaseIdx));
    }//end solve

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
        using MaterialTraits = NullMaterialTraits<Scalar, numPhases>;
        using MaterialLaw = NullMaterial<MaterialTraits>;
        using MaterialLawParams = typename MaterialLaw::Params;

        MaterialLawParams matParams;
        solve<MaterialLaw>(fluidState, matParams, globalMolarities, tolerance);
    }


protected:
    template <class FlashFluidState>
    static typename FlashFluidState::Scalar wilsonK_(const FlashFluidState& fluidState, int compIdx)
    {
        using FlashEval = typename FlashFluidState::Scalar;
        const auto& acf = FluidSystem::acentricFactor(compIdx);
        const auto& T_crit = FluidSystem::criticalTemperature(compIdx);
        const auto& T = fluidState.temperature(0);
        const auto& p_crit = FluidSystem::criticalPressure(compIdx);
        const auto& p = fluidState.pressure(0); //for now assume no capillary pressure

        const auto& tmp = Opm::exp(5.3727 * (1+acf) * (1-T_crit/T)) * (p_crit/p);
        return tmp;
    }

    template <class Vector>
    static typename Vector::field_type rachfordRice_g_(const Vector& K, const Evaluation L, const Vector& globalComposition)
    {
        typename Vector::field_type g=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            g += (globalComposition[compIdx]*(K[compIdx]-1))/(K[compIdx]-L*(K[compIdx]-1));
        }
        return g;
    }

    template <class Vector>
    static typename Vector::field_type rachfordRice_dg_dL_(const Vector& K, const Evaluation L, const Vector& globalComposition)
    {
        typename Vector::field_type dg=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            dg += (globalComposition[compIdx]*(K[compIdx]-1)*(K[compIdx]-1))/((K[compIdx]-L*(K[compIdx]-1))*(K[compIdx]-L*(K[compIdx]-1)));
        }
        return dg;
    }

    template <class Vector>
    static typename Vector::field_type solveRachfordRice_g_(const Vector& K, const Vector& globalComposition)
    {
        // Find min and max K. Have to do a laborious for loop to avoid water component (where K=0)
        // TODO: Replace loop with Dune::min_value() and Dune::max_value() when water component is properly handled
        Evaluation Kmin = K[0];
        Evaluation Kmax = K[0];
        for (int compIdx=1; compIdx<numComponents; ++compIdx){
            if (K[compIdx] < Kmin)
                Kmin = K[compIdx];
            else if (K[compIdx] >= Kmax)
                Kmax = K[compIdx];
        }
        // Lower and upper bound for solution
        Evaluation Lmin = (Kmin / (Kmin - 1));
        Evaluation Lmax = Kmax / (Kmax - 1);

        // Check if Lmin < Lmax, and switch if not
        if (Lmin > Lmax)
        {
            Evaluation Ltmp = Lmin;
            Lmin = Lmax;
            Lmax = Ltmp;
        }

        // Initial guess
        Evaluation L = (Lmin + Lmax)/2;

        // Newton-Rahpson loop
        for (int iteration=0; iteration<200; ++iteration){
            // Calculate function and derivative values
            Evaluation g = rachfordRice_g_(K, L, globalComposition);
            Evaluation dg_dL = rachfordRice_dg_dL_(K, L, globalComposition);

            // Lnew = Lold - g/dg;
            Evaluation delta = g/dg_dL;
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
        throw std::runtime_error(" Rachford-Rice did not converge within maximum number of iterations" );
    }

    template <class Vector>
    static typename Vector::field_type bisection_g_(const Vector& K, Evaluation Lmin, Evaluation Lmax, const Vector& globalComposition)
    {
        // Check if g(Lmin) and g(Lmax) have opposite sign
        Evaluation gLmin = rachfordRice_g_(K, Lmin, globalComposition);
        Evaluation gLmax = rachfordRice_g_(K, Lmax, globalComposition);
        if (Dune::sign(gLmin) == Dune::sign(gLmax))
        {
            throw std::runtime_error("Lmin and Lmax are incorrect for bisection");
        }
            
        // Bisection loop
        for (int iteration=0; iteration<200; ++iteration){
            // New midpoint
            Evaluation L = (Lmin + Lmax) / 2;
            Evaluation gMid = rachfordRice_g_(K, L, globalComposition);

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
        Evaluation S_l, S_v;

        // Use equilibrium constants made with Wilson's formula
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
                fluidState.setMoleFraction(gasPhaseIdx, compIdx, globalComposition[compIdx]);
                fluidState.setMoleFraction(oilPhaseIdx, compIdx, globalComposition[compIdx]);
            }
        }
        // If it is not stable, we use the mole fractions from the Michelsen test to calculate a new K
        else {
            for (int compIdx = 0; compIdx<numComponents; ++compIdx) {
                K[compIdx] = y[compIdx] / x[compIdx];
            }
            // std::cout << "K = " << K << std::endl;
        }
    }

    template <class FlashFluidState, class ComponentVector>
    static void checkStability_(const FlashFluidState& fluidState, bool& isTrivial, ComponentVector& K, ComponentVector& xy_loc, Evaluation& S_loc, const ComponentVector& globalComposition, bool isGas)
    {
        using FlashEval = typename FlashFluidState::Scalar;
        using ThisType = TwoPhaseThreeComponentFluidSystem<Scalar>;
        using PengRobinsonMixture = typename Opm::PengRobinsonMixture<Scalar, ThisType>;

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
                    xy_loc[compIdx] = globalComposition[compIdx]/K[compIdx];
                    S_loc += xy_loc[compIdx];
                }
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    xy_loc[compIdx] /= S_loc;
                    fluidState_fake.setMoleFraction(oilPhaseIdx, compIdx, xy_loc[compIdx]);
                }
            }

            int phaseIdx = (isGas?gasPhaseIdx:oilPhaseIdx);
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                fluidState_global.setMoleFraction(phaseIdx, compIdx, globalComposition[compIdx]);
            }

            typename FluidSystem::template ParameterCache<FlashEval> paramCache_fake;
            paramCache_fake.updatePhase(fluidState_fake, phaseIdx);

            typename FluidSystem::template ParameterCache<FlashEval> paramCache_global;
            paramCache_global.updatePhase(fluidState_global, phaseIdx);

            //fugacity for fake phases each component
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                auto phiFake = PengRobinsonMixture::computeFugacityCoefficient(fluidState_fake, paramCache_fake, phaseIdx, compIdx);
                auto phiGlobal = PengRobinsonMixture::computeFugacityCoefficient(fluidState_global, paramCache_global, phaseIdx, compIdx);

                fluidState_fake.setFugacityCoefficient(phaseIdx, compIdx, phiFake);
                fluidState_global.setFugacityCoefficient(phaseIdx, compIdx, phiGlobal);
            }

           
            ComponentVector R;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (isGas){
                    auto fug_fake = fluidState_fake.fugacity(gasPhaseIdx, compIdx);
                    auto fug_global = fluidState_global.fugacity(gasPhaseIdx, compIdx);
                    auto fug_ratio = fug_global / fug_fake;
                    R[compIdx] = fug_ratio / S_loc;
                }
                else{
                    auto fug_fake = fluidState_fake.fugacity(oilPhaseIdx, compIdx);
                    auto fug_global = fluidState_global.fugacity(oilPhaseIdx, compIdx);
                    auto fug_ratio = fug_fake / fug_global;
                    R[compIdx] = fug_ratio * S_loc;
                }
            }

            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                K[compIdx] *= R[compIdx];
            }
            Scalar R_norm = 0.0;
            Scalar K_norm = 0.0;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                auto a = Opm::getValue(R[compIdx]) - 1.0;
                auto b = Opm::log(Opm::getValue(K[compIdx]));
                R_norm += a*a;
                K_norm += b*b;
            }
            isTrivial = (K_norm < 1e-5);
            if (isTrivial || R_norm < 1e-10)
                return;
            //todo: make sure that no molefraction is smaller than 1e-8 ?
            //todo: take care of water!
        }
        throw std::runtime_error(" Stability test did not converge");
    }//end checkStability

    template <class FlashFluidState, class ComponentVector>
    static void computeLiquidVapor_(FlashFluidState& fluidState, Evaluation& L, ComponentVector& K, const ComponentVector& globalComposition)
    {
        // Calculate x and y, and normalize
        ComponentVector x;
        ComponentVector y;
        Evaluation sumx=0;
        Evaluation sumy=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            x[compIdx] = globalComposition[compIdx]/(L + (1-L)*K[compIdx]);
            sumx += x[compIdx];
            y[compIdx] = (K[compIdx]*globalComposition[compIdx])/(L + (1-L)*K[compIdx]);
            sumy += y[compIdx];
        }
        x /= sumx;
        y /= sumy;

        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            fluidState.setMoleFraction(oilPhaseIdx, compIdx, x[compIdx]);
            fluidState.setMoleFraction(gasPhaseIdx, compIdx, y[compIdx]);
        }
    }

    template <class FlashFluidState, class ComponentVector>
    static void newtonCompositionUpdate_(ComponentVector& K, Evaluation& L, FlashFluidState& fluidState, const ComponentVector& globalComposition,
                                         int verbosity)
    {
        // Newton declarations
        using NewtonVector = Dune::FieldVector<Evaluation, numMisciblePhases*numMiscibleComponents+1>;
        using NewtonMatrix = Dune::FieldMatrix<Evaluation, numMisciblePhases*numMiscibleComponents+1, numMisciblePhases*numMiscibleComponents+1>;
        NewtonVector newtonX;
        NewtonVector newtonB;
        NewtonMatrix newtonA;
        NewtonVector newtonDelta;

        // Header for information output
        if (verbosity == 1) {
            std::cout << std::setw(10) << "Iteration" << std::setw(16) << "Norm2(step)" << std::setw(16) << "Norm2(Residual)" << std::endl;
        }

        // Compute x and y from K, L and Z
        computeLiquidVapor_(fluidState, L, K, globalComposition);

        // Assign primary variables (x, y and L)
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            newtonX[compIdx] = Opm::getValue(fluidState.moleFraction(oilPhaseIdx, compIdx));
            newtonX[compIdx + numMiscibleComponents] = Opm::getValue(fluidState.moleFraction(gasPhaseIdx, compIdx));
        }
        newtonX[numMiscibleComponents*numMiscibleComponents] = Opm::getValue(L);

        // 
        // Main Newton loop
        // 
        bool convFug = false; /* for convergence check */
        for (int i = 0; i< 100; ++i){
            // Evaluate residuals (newtonB)
            evalDefect_(newtonB, newtonX, fluidState, globalComposition);

            // Print iteration info
            if (verbosity == 1) {
                if (i == 0) {
                    std::cout << std::setw(10) << i << std::setw(16) << "N/A" << std::setw(16)  << newtonB.two_norm() << std::endl;
                }
                else {
                    std::cout << std::setw(10) << i << std::setw(16) << newtonDelta.two_norm() << std::setw(16) << newtonB.two_norm() << std::endl;
                }
            }

            // Check fugacity equilibrium for convergence
            convFug = checkFugacityEquil_(newtonB);

            // If convergence have been met, we abort; else we update step and loop once more
            if (convFug == true) {
                // Extract x, y and L together with calculation of K
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    fluidState.setMoleFraction(oilPhaseIdx, compIdx, newtonX[compIdx]);
                    fluidState.setMoleFraction(gasPhaseIdx, compIdx, newtonX[compIdx + numMiscibleComponents]);
                    K[compIdx] = newtonX[compIdx + numMiscibleComponents] / newtonX[compIdx];
                }
                L = newtonX[numMiscibleComponents*numMiscibleComponents];

                // Print info
                if (verbosity == 1) {
                    std::cout << "Solution converged to the following result :" << std::endl;
                    std::cout << std::setw(10) << "Variable" << std::setw(16) << "Comp. 0" << std::setw(16) << "Comp. 1" << std::endl;
                    std::cout << std::setw(10) << "x" <<  std::setw(16) << fluidState.moleFraction(oilPhaseIdx, Comp1Idx) << std::setw(16) << fluidState.moleFraction(oilPhaseIdx, Comp0Idx) << std::endl;
                    std::cout << std::setw(10) << "y" << std::setw(16) << fluidState.moleFraction(gasPhaseIdx, Comp1Idx) << std::setw(16) << fluidState.moleFraction(gasPhaseIdx, Comp0Idx) << std::endl; 
                    std::cout << std::setw(10) << "K" << std::setw(16) << K[Comp0Idx] << std::setw(16) << K[Comp1Idx] << std::endl;
                    std::cout << "L = " << L << std::endl;
                }
                return;
            }
            else {
                // Calculate Jacobian (newtonA)
                evalJacobian_(newtonA, newtonX, fluidState, globalComposition);
              
                // Solve system J * x = -r, which in our case is newtonA*newtonX = newtonB, to get next step (newtonDelta) 
                newtonA.solve(newtonDelta, newtonB);

                // Update current solution (newtonX) with simple relaxation method (percentage of step applied)
                updateCurrentSol_(newtonX, newtonDelta);
            }
        }
        throw std::runtime_error(" Newton composition update did not converge within maxIterations");
    }

    template <class DefectVector>
    static void updateCurrentSol_(DefectVector& x, DefectVector& d)
    {
        // Find smallest percentage update
        Scalar w = 1.0;
        for (int i=0; i<x.size(); ++i){
            Scalar w_tmp = Opm::getValue(Opm::min(Opm::max(x[i] + d[i], 0.0), 1.0) - x[i]) / Opm::getValue(d[i]);
            w = Opm::min(w, w_tmp);
        }

        // Loop over the solution vector and apply the smallest percentage update
        for (int i=0; i<x.size(); ++i){
            x[i] += w*d[i];
        }
    }

    template <class DefectVector>
    static bool checkFugacityEquil_(DefectVector& b)
    {
        // Init. fugacity vector
        DefectVector fugVec;

        // Loop over b and find the fugacity equilibrium
        // OBS: If the equations in b changes in evalDefect_ then you have to change here as well!
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            fugVec[compIdx] = 0.0;
            fugVec[compIdx+numMiscibleComponents] = b[compIdx+numMiscibleComponents];
        }
        fugVec[numMiscibleComponents*numMisciblePhases] = 0.0;
        
        // Check if norm(fugVec) is less than tolerance
        bool conv = fugVec.two_norm() < 1e-6;
        return conv;
    }

    template <class FluidState, class DefectVector, class ComponentVector>
    static void evalDefect_(DefectVector& b,
                            DefectVector& x,
                            const FluidState& fluidStateIn,
                            const ComponentVector& globalComposition)
    {
        // Put x and y in a FluidState instance for fugacity calculations
        FluidState fluidState(fluidStateIn);
        ComponentVector K;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            fluidState.setMoleFraction(oilPhaseIdx, compIdx, x[compIdx]);
            fluidState.setMoleFraction(gasPhaseIdx, compIdx, x[compIdx + numMiscibleComponents]);
        }
        

        // Compute fugacities
        using ParamCache = typename FluidSystem::template ParameterCache<typename FluidState::Scalar>;
        ParamCache paramCache;
        for (int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
            paramCache.updatePhase(fluidState, phaseIdx);
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                Evaluation phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }

        // Compute residuals for Newton update. Primary variables are: x, y, and L
        // TODO: Make this AD
        // Extract L
        Evaluation L = x[numMiscibleComponents*numMisciblePhases];
        
        // Residuals
        // OBS: the residuals are negative in the newton system!
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            // z - L*x - (1-L) * y
            b[compIdx] = -globalComposition[compIdx] + L*x[compIdx] + (1-L)*x[compIdx + numMiscibleComponents];
            
            // (f_liquid/f_vapor) - 1 = 0
            b[compIdx + numMiscibleComponents] = -(fluidState.fugacity(oilPhaseIdx, compIdx) / fluidState.fugacity(gasPhaseIdx, compIdx)) + 1.0;
        
            // sum(x) - sum(y) = 0
            b[numMiscibleComponents*numMisciblePhases] += -x[compIdx] + x[compIdx + numMiscibleComponents];
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
        DefectVector x;
        DefectVector b0;
        for(int j=0; j<xIn.size(); ++j){
            x[j] = xIn[j];
        }
        
        evalDefect_(b0, x, fluidStateIn, globalComposition);

        // Make the jacobian A in Newton system Ax=b
        Scalar epsilon = 1e-10;
        for(int i=0; i<b0.size(); ++i){
            // Permutate x and calculate response
            x[i] += epsilon;
            DefectVector bEps;
            evalDefect_(bEps, x, fluidStateIn, globalComposition);
            x[i] -= epsilon;

            // Forward difference of all eqs wrt primary variable i
            DefectVector derivI;
            for(int j=0; j<b0.size(); ++j){
                derivI[j] = bEps[j];
                derivI[j] -= b0[j];
                derivI[j] /= epsilon;
                A[j][i] = derivI[j];
            }
        }
    }//end evalJacobian
    
    template <class FlashFluidState, class ComponentVector>
    static void successiveSubstitutionComposition_(ComponentVector& K, Evaluation& L, FlashFluidState& fluidState, const ComponentVector& globalComposition, 
                                                   bool standAlone, int verbosity)
    {
        // Determine max. iterations based on if it will be used as a standalone flash or as a pre-process to Newton (or other) method.
        int maxIterations;
        if (standAlone == true)
            maxIterations = 100;
        else
            maxIterations = 10;

        // 
        // Successive substitution loop
        // 
        for (int i=0; i < maxIterations; ++i){
            // Compute (normalized) liquid and vapor mole fractions
            computeLiquidVapor_(fluidState, L, K, globalComposition);

            // Calculate fugacity coefficient
            using ParamCache = typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>;
            ParamCache paramCache;
            for (int phaseIdx=0; phaseIdx<numPhases; ++phaseIdx){
                paramCache.updatePhase(fluidState, phaseIdx);
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    Evaluation phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
            
            // Calculate fugacity ratio
            ComponentVector fugRatio;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                fugRatio[compIdx] = (fluidState.fugacity(oilPhaseIdx, compIdx)/fluidState.fugacity(gasPhaseIdx, compIdx)) - 1.0;
            }

            // Check convergence
            if (fugRatio.two_norm() < 1e-6){
                // Print info
                if (verbosity == 1) {
                    std::cout << "Solution converged to the following result :" << std::endl;
                    std::cout << std::setw(10) << "Variable" << std::setw(16) << "Comp. 0" << std::setw(16) << "Comp. 1" << std::endl;
                    std::cout << std::setw(10) << "x" <<  std::setw(16) << fluidState.moleFraction(oilPhaseIdx, Comp1Idx) << std::setw(16) << fluidState.moleFraction(oilPhaseIdx, Comp0Idx) << std::endl;
                    std::cout << std::setw(10) << "y" << std::setw(16) << fluidState.moleFraction(gasPhaseIdx, Comp1Idx) << std::setw(16) << fluidState.moleFraction(gasPhaseIdx, Comp0Idx) << std::endl; 
                    std::cout << std::setw(10) << "K" << std::setw(16) << K[Comp0Idx] << std::setw(16) << K[Comp1Idx] << std::endl;
                    std::cout << "L = " << L << std::endl;
                }
                return;
            }
            
            //  If convergence is not met, K is updated in a successive substitution manner
            else{
                // Update K
                for (int compIdx=0; compIdx<numComponents; ++compIdx){
                    K[compIdx] *= (fugRatio[compIdx] + 1.0);
                }

                // Solve Rachford-Rice to get L from updated K
                L = solveRachfordRice_g_(K, globalComposition);
            }

        }
        throw std::runtime_error("Successive substitution composition update did not converge within maxIterations");
    }
    
};//end ChiFlash

} // namespace Opm

#endif
