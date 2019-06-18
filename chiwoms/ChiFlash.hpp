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
    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static void guessInitial(FluidState& fluidState,
                             const Dune::FieldVector<Evaluation, numComponents>& globalMolarities)
    {
        // water saturation
        Evaluation brineMass = globalMolarities[BrineIdx] * FluidSystem::molarMass(BrineIdx);
        Evaluation waterSaturation = Opm::min(1.0, brineMass/1000); //mass/density

        fluidState.setSaturation(waterPhaseIdx, waterSaturation);
        // oil and gas saturation
        fluidState.setSaturation(oilPhaseIdx, (1.0-waterSaturation)/2.0);
        fluidState.setSaturation(gasPhaseIdx, (1.0-waterSaturation)/2.0);


        // the sum of all molarities
        Evaluation sumMoles = 0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (compIdx == BrineIdx)
                    continue;
            sumMoles += globalMolarities[compIdx];
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (phaseIdx == waterPhaseIdx)
                continue;
            // composition
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx){
                if (compIdx == BrineIdx){
                    fluidState.setMoleFraction(phaseIdx,
                                               compIdx,
                                               0.0);
                } else {
                    fluidState.setMoleFraction(phaseIdx,
                                               compIdx,
                                               globalMolarities[compIdx]/sumMoles);
                }
            }

            // pressure. use atmospheric pressure as initial guess
            fluidState.setPressure(phaseIdx, 20.0e6);
        }

        // set composition in water to only water
        fluidState.setMoleFraction(waterPhaseIdx,
                                   BrineIdx,
                                   1.0);
        fluidState.setMoleFraction(waterPhaseIdx,
                                   OctaneIdx,
                                   0.0);
        fluidState.setMoleFraction(waterPhaseIdx,
                                   CO2Idx,
                                   0.0);

        // set the fugacity coefficients of all components in all phases
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updateAll(fluidState);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (phaseIdx == waterPhaseIdx)
                continue;
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
                if (compIdx == BrineIdx)
                    continue;
                const typename FluidState::Scalar phi =
                    FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }
    }

    /*!
     * \brief Calculates the fluid state from the total mass of the components
     *
     */
    template <class MaterialLaw, class FluidState>
    static void solve(FluidState& fluidState,
                      const typename MaterialLaw::Params& matParams,
                      typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                      const Dune::FieldVector<typename FluidState::Scalar, numComponents>& globalMolarities,
                      Scalar tolerance = -1.0)
    {
        typedef typename FluidState::Scalar InputEval;

        typedef Dune::FieldMatrix<InputEval, numEq, numEq> Matrix;
        typedef Dune::FieldVector<InputEval, numEq> Vector;
        typedef Opm::DenseAd::Evaluation</*Scalar=*/InputEval,
                                         /*numDerivs=*/numEq> FlashEval;

        typedef Dune::FieldVector<FlashEval, numEq> FlashDefectVector;
        typedef Opm::CompositionalFluidState<FlashEval, FluidSystem, /*energy=*/false> FlashFluidState;
        typedef Dune::FieldVector<typename FluidState::Scalar, numComponents> ComponentVector;

#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,7)
        Dune::FMatrixPrecision<InputEval>::set_singular_limit(1e-35);
#endif

        if (tolerance <= 0)
            tolerance = std::min<Scalar>(1e-3,
                                         1e8*std::numeric_limits<Scalar>::epsilon());

        typename FluidSystem::template ParameterCache<FlashEval> flashParamCache;
        flashParamCache.assignPersistentData(paramCache);



        //initial guess for the K value, L value and global composition z in RachfordRice
        ComponentVector K;
        ComponentVector z;
        Scalar L;
        Scalar L_min = 1e10;
        Scalar L_max = -1e10;
        Scalar totalMoles = 0;

        for (int compIdx = 0; compIdx<numComponents; ++compIdx) {
            K[compIdx] = wilsonK_(fluidState, compIdx);
            L_min = Opm::min(L_min, 1/(1-K[compIdx]));
            L_max = Opm::max(L_max, 1/(1-K[compIdx]));

            totalMoles += globalMolarities[compIdx];

        }
        L = (L_min + L_max)/2;
        z = globalMolarities;
        z /= totalMoles;

        //Rachford Rice equation
        L = solveRachfordRice_g_(K, L, z);

        //Phase stability test
        bool isStable;
        ComponentVector x;
        ComponentVector y;
        phaseStabilityTest_(isStable, x, y, fluidState, z);

        /////////////////////////
        // Newton method
        /////////////////////////

        // Jacobian matrix
        Matrix J;
        // solution, i.e. phase composition
        Vector deltaX;
        // right hand side
        Vector b;

        Valgrind::SetUndefined(J);
        Valgrind::SetUndefined(deltaX);
        Valgrind::SetUndefined(b);

        FlashFluidState flashFluidState;
        assignFlashFluidState_<MaterialLaw>(fluidState, flashFluidState, matParams, flashParamCache);

        // copy the global molarities to a vector of evaluations. Remember that the
        // global molarities are constants. (but we need to copy them to a vector of
        // FlashEvals anyway in order to avoid getting into hell's kitchen.)
        Dune::FieldVector<FlashEval, numMiscibleComponents> flashGlobalMolarities;
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx){
            if (compIdx == BrineIdx)
                continue;
            flashGlobalMolarities[compIdx] = globalMolarities[compIdx];
        }

        FlashDefectVector defect;
        const unsigned nMax = 50; // <- maximum number of newton iterations
        for (unsigned nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate the defect of the flash equations and their derivatives
            evalDefect_(defect, flashFluidState, flashGlobalMolarities);
            Valgrind::CheckDefined(defect);

            // create field matrices and vectors out of the evaluation vector to solve
            // the linear system of equations.
            for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
                for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx)
                    J[eqIdx][pvIdx] = defect[eqIdx].derivative(pvIdx);

                b[eqIdx] = defect[eqIdx].value();
            }
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

            // Solve J*x = b
            deltaX = 0.0;
            try { J.solve(deltaX, b); }
            catch (const Dune::FMatrixError& e) {
                throw Opm::NumericalIssue(e.what());
            }
            Valgrind::CheckDefined(deltaX);

            // update the fluid quantities.
            Scalar relError = update_<MaterialLaw>(globalMolarities, flashFluidState, matParams, flashParamCache, deltaX);

            if (relError < tolerance) {
                assignOutputFluidState_(flashFluidState, fluidState);
                return;
            }
        }

        std::ostringstream oss;
        oss << "ChiFlash solver failed:"
            << " {c_alpha^kappa} = {" << globalMolarities << "}, "
            << " T = " << fluidState.temperature(/*phaseIdx=*/0);
        throw NumericalIssue(oss.str());
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
    template <class FluidState>
    static void printFluidState_(const FluidState& fluidState)
    {
        typedef typename FluidState::Scalar FsScalar;

        std::cout << "saturations: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.saturation(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "pressures: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.pressure(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "densities: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.density(phaseIdx) << " ";
        std::cout << "\n";

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "composition " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fluidState.moleFraction(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "fugacities " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fluidState.fugacity(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        std::cout << "global component molarities: ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            FsScalar sum = 0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                sum += fluidState.saturation(phaseIdx)*fluidState.molarity(phaseIdx, compIdx);
            }
            std::cout << sum << " ";
        }
        std::cout << "\n";
    }

    template <class MaterialLaw, class InputFluidState, class FlashFluidState>
    static void assignFlashFluidState_(const InputFluidState& inputFluidState,
                                       FlashFluidState& flashFluidState,
                                       const typename MaterialLaw::Params& matParams,
                                       typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& flashParamCache)
    {
        typedef typename FlashFluidState::Scalar FlashEval;

        // copy the temperature: even though the model which uses the flash solver might
        // be non-isothermal, the flash solver does not consider energy. (it could be
        // modified to do so relatively easily, but it would come at increased
        // computational cost and normally temperature instead of "total internal energy
        // of the fluids" is specified.)
        flashFluidState.setTemperature(inputFluidState.temperature(/*phaseIdx=*/0));

        // copy the saturations: the first N-1 phases are primary variables, the last one
        // is one minus the sum of the former.
        FlashEval Slast = 1.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            FlashEval S = inputFluidState.saturation(phaseIdx);
            S.setDerivative(S0PvIdx + phaseIdx, 1.0);

            Slast -= S;

            flashFluidState.setSaturation(phaseIdx, S);
        }
        flashFluidState.setSaturation(numPhases - 1, Slast);

        // copy the pressures: the first pressure is the first primary variable, the
        // remaining ones are given as p_beta = p_alpha + p_calpha,beta
        FlashEval p0 = inputFluidState.pressure(0);
        p0.setDerivative(p0PvIdx, 1.0);

        std::array<FlashEval, numPhases> pc;
        MaterialLaw::capillaryPressures(pc, matParams, flashFluidState);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            flashFluidState.setPressure(phaseIdx, p0 + (pc[phaseIdx] - pc[0]));

        // copy the mole fractions: all of them are primary variables
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                FlashEval x = inputFluidState.moleFraction(phaseIdx, compIdx);
                x.setDerivative(x00PvIdx + phaseIdx*numComponents + compIdx, 1.0);
                flashFluidState.setMoleFraction(phaseIdx, compIdx, x);
            }
        }

        flashParamCache.updateAll(flashFluidState);

        // compute the density of each phase and the fugacity coefficient of each
        // component in each phase.
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const FlashEval& rho = FluidSystem::density(flashFluidState, flashParamCache, phaseIdx);
            flashFluidState.setDensity(phaseIdx, rho);

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                const FlashEval& fugCoeff = FluidSystem::fugacityCoefficient(flashFluidState, flashParamCache, phaseIdx, compIdx);
                flashFluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }
    }

    template <class FlashFluidState, class OutputFluidState>
    static void assignOutputFluidState_(const FlashFluidState& flashFluidState,
                                        OutputFluidState& outputFluidState)
    {
        outputFluidState.setTemperature(flashFluidState.temperature(/*phaseIdx=*/0).value());

        // copy the saturations, pressures and densities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const auto& S = flashFluidState.saturation(phaseIdx).value();
            outputFluidState.setSaturation(phaseIdx, S);

            const auto& p = flashFluidState.pressure(phaseIdx).value();
            outputFluidState.setPressure(phaseIdx, p);

            const auto& rho = flashFluidState.density(phaseIdx).value();
            outputFluidState.setDensity(phaseIdx, rho);
        }

        // copy the mole fractions and fugacity coefficients
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                const auto& moleFrac =
                    flashFluidState.moleFraction(phaseIdx, compIdx).value();
                outputFluidState.setMoleFraction(phaseIdx, compIdx, moleFrac);

                const auto& fugCoeff =
                    flashFluidState.fugacityCoefficient(phaseIdx, compIdx).value();
                outputFluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }
    }

    template <class FlashFluidState>
    static typename FlashFluidState::Scalar wilsonK_(const FlashFluidState& fluidState, int compIdx)
    {
        typedef typename FlashFluidState::Scalar FlashEval;
        const auto& acf = FluidSystem::acentricFactor(compIdx);
        const auto& T_crit = FluidSystem::criticalTemperature(compIdx);
        const auto& T = fluidState.temperature(0);
        const auto& p_crit = FluidSystem::criticalPressure(compIdx);
        const auto& p = fluidState.pressure(0); //for now assume no capillary pressure

        return Opm::exp(5.3727 * (1+acf) * (1-T_crit/T) * p_crit/p);

    }

    template <class Vector>
    static typename Vector::field_type rachfordRice_g_(const Vector& K, const Scalar L, const Vector& z)
    {
        typename Vector::field_type g=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            g += (z[compIdx]*(K[compIdx]-1))/(1+L*(K[compIdx]-1));
        }
        return g;
    }

    template <class Vector>
    static typename Vector::field_type rachfordRice_dg_dL_(const Vector& K, const Scalar L, const Vector& z)
    {
        typename Vector::field_type dg=0;
        for (int compIdx=0; compIdx<numComponents; ++compIdx){
            dg += -(z[compIdx]*(K[compIdx]-1)*(K[compIdx]-1))/((1+L*(K[compIdx]-1))*(1+L*(K[compIdx]-1)));
        }
        return dg;
    }

    template <class Vector>
    static typename Vector::field_type solveRachfordRice_g_(const Vector& K, Scalar L, const Vector& z)
    {
        for (int iteration=0; iteration<100; ++iteration){
            //L=Lold+g/dg;
            Scalar g = rachfordRice_g_(K, L, z);
            Scalar dg_dL = rachfordRice_dg_dL_(K, L, z);
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
        bool V_stable = (S_v < (1.0 + 1e-5)) || isTrivialV;

        checkStability_(fluidState, isTrivialL, K_l, x, S_l, globalComposition, /*isGas=*/false);
        bool L_stable = (S_l < (1.0 + 1e-5)) || isTrivialL;

        isStable = L_stable && V_stable; //todo: understand this. should maybe called {L,V}_unstable??
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

        //make two fake phases inside one phase and check for positive volume
        FlashFluidState fluidState_fake = fluidState;
        FlashFluidState fluidState_global = fluidState;

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            K[compIdx] = wilsonK_(fluidState, compIdx);

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

                Scalar phiFake = FluidSystem::fugacityCoefficient(fluidState_fake, paramCache_fake, phaseIdx, compIdx);
                Scalar phiGlobal = FluidSystem::fugacityCoefficient(fluidState_global, paramCache_global, phaseIdx, compIdx);

                fluidState_fake.setFugacityCoefficient(phaseIdx, compIdx, phiFake);
                fluidState_global.setFugacityCoefficient(phaseIdx, compIdx, phiGlobal);
            }

            ComponentVector R;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                if (isGas)
                    R[compIdx] = fluidState_global.fugacity(oilPhaseIdx, compIdx)/fluidState_fake.fugacity(oilPhaseIdx, compIdx)/S_loc;
                else
                    R[compIdx] = fluidState_fake.fugacity(gasPhaseIdx, compIdx)/fluidState_global.fugacity(gasPhaseIdx, compIdx)*S_loc;
            }

            for (int compIdx=0; compIdx<numComponents; ++compIdx)
                K[compIdx] *= R[compIdx];

            Scalar R_norm = 0.0;
            Scalar K_norm = 0.0;
            for (int compIdx=0; compIdx<numComponents; ++compIdx){
                auto a = R[compIdx] - 1.0;
                auto b = Opm::log(K[compIdx]);

                R_norm += a*a;
                K_norm += b*b;
            }

            isTrivial = (K_norm < 1e-5);
            if (isTrivial || R_norm < 1e-10)
                return;
            //note: make sure that no molefraction is smaller than 1e-8 ?
            //note: take care of water!
        }

        throw std::runtime_error("stability test did not converge");
    }


    template <class FlashFluidState, class FlashDefectVector, class FlashComponentVector>
    static void evalDefect_(FlashDefectVector& b,
                            const FlashFluidState& fluidState,
                            const FlashComponentVector& globalMolarities)
    {
        typedef typename FlashFluidState::Scalar FlashEval;

        unsigned eqIdx = 0;


        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (compIdx == BrineIdx)
                continue;


#warning WTF
//            b[eqIdx] =
//                    compute_W_oilPhase(fluidState, compIdx)
//                    - compute_W_gasPhase(fluidState, compIdx);
            ++eqIdx;
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (phaseIdx == waterPhaseIdx)
                continue;

#warning WTF2
        }
        assert(eqIdx == numMiscibleComponents*numMisciblePhases);

        // the fact saturations must sum up to 1 is included implicitly and also,
        // capillary pressures are treated implicitly!

        // global molarities of the miscible components/phases are given
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (compIdx==BrineIdx)
                continue;
            b[eqIdx] = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (phaseIdx == waterPhaseIdx)
                    continue;
                b[eqIdx] +=
                    fluidState.saturation(phaseIdx)
                    * fluidState.molarity(phaseIdx, compIdx);
            }

            b[eqIdx] -= globalMolarities[compIdx];
            ++eqIdx;
        }

        assert(eqIdx == numEq);
    }

    template <class MaterialLaw, class FlashFluidState, class FlashEvalVector, class EvalVector>
    static Scalar update_(const EvalVector& globalMolarities,
                          FlashFluidState& fluidState,
                          const typename MaterialLaw::Params& matParams,
                          typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& paramCache,
                          const FlashEvalVector& deltaX)
    {
        // note that it is possible that FlashEval::Scalar is an Evaluation itself
        typedef typename FlashFluidState::Scalar FlashEval;
        typedef typename FlashEval::ValueType InnerEval;

#ifndef NDEBUG
        // make sure we don't swallow non-finite update vectors
        assert(deltaX.dimension == numEq);
        for (unsigned i = 0; i < numEq; ++i)
            assert(std::isfinite(Opm::scalarValue(deltaX[i])));
#endif

        Scalar relError = 0;
        for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            FlashEval tmp = getQuantity_(fluidState, pvIdx);
            InnerEval delta = deltaX[pvIdx];

            relError = std::max(relError,
                                std::abs(Opm::scalarValue(delta))
                                * quantityWeight_(fluidState, pvIdx));

            if (isSaturationIdx_(pvIdx)) {
                // dampen to at most 25% change in saturation per iteration
                delta = Opm::min(0.25, Opm::max(-0.25, delta));
            }
            else if (isMoleFracIdx_(pvIdx)) {
                // dampen to at most 20% change in mole fraction per iteration
                delta = Opm::min(0.20, Opm::max(-0.20, delta));
            }
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 50% change in pressure per iteration
                delta = Opm::min(0.5*fluidState.pressure(0).value(),
                                 Opm::max(-0.5*fluidState.pressure(0).value(),
                                          delta));
            }

            tmp -= delta;
            setQuantity_(fluidState, pvIdx, tmp);
        }

        completeFluidState_<MaterialLaw>(globalMolarities, fluidState, paramCache, matParams);

        return relError;
    }

    template <class MaterialLaw, class FlashFluidState, class EvalVector>
    static void completeFluidState_(const EvalVector& globalMolarities,
                                    FlashFluidState& flashFluidState,
                                    typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& paramCache,
                                    const typename MaterialLaw::Params& matParams)
    {
        typedef typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar> ParamCache;

        typedef typename FlashFluidState::Scalar FlashEval;


        // calculate the saturation of the last phase as a function of
        // the other saturations

        // water saturation
        FlashEval brineMass = globalMolarities[BrineIdx] * FluidSystem::molarMass(BrineIdx);
        FlashEval waterSaturation = Opm::min(1.0, brineMass/1000.0); //mass/density
        const FlashEval& oilSaturation = flashFluidState.saturation(oilPhaseIdx);
        const FlashEval& gasSaturation = 1 - waterSaturation - oilSaturation;

        flashFluidState.setSaturation(waterPhaseIdx, waterSaturation);
        flashFluidState.setSaturation(oilPhaseIdx, oilSaturation);
        flashFluidState.setSaturation(gasPhaseIdx, gasSaturation);

        // update the pressures using the material law (saturations
        // and first pressure are already set because it is implicitly
        // solved for.)
        Dune::FieldVector<FlashEval, numPhases> pC;
        MaterialLaw::capillaryPressures(pC, matParams, flashFluidState);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            flashFluidState.setPressure(phaseIdx,
                                        flashFluidState.pressure(oilPhaseIdx)
                                        + (pC[phaseIdx] - pC[oilPhaseIdx]));

        // update the parameter cache
        paramCache.updateAll(flashFluidState, /*except=*/ParamCache::Temperature);

        // update all densities and fugacity coefficients
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            const FlashEval& rho = FluidSystem::density(flashFluidState, paramCache, phaseIdx);
            flashFluidState.setDensity(phaseIdx, rho);

            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
                const FlashEval& phi = FluidSystem::fugacityCoefficient(flashFluidState, paramCache, phaseIdx, compIdx);
                flashFluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }
    }

    static bool isPressureIdx_(unsigned pvIdx)
    { return pvIdx == p0PvIdx; }

    static bool isSaturationIdx_(unsigned pvIdx)
    { return pvIdx == S0PvIdx;}

    static bool isMoleFracIdx_(unsigned pvIdx)
    { return x00PvIdx <= pvIdx && pvIdx < x00PvIdx + numMisciblePhases*numMiscibleComponents; }

    // retrieves a quantity from the fluid state
    template <class FluidState>
    static const typename FluidState::Scalar& getQuantity_(const FluidState& fluidState, unsigned pvIdx)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx == p0PvIdx)
            return fluidState.pressure(oilPhaseIdx);
        // first saturation
        else if (pvIdx == S0PvIdx)
            return fluidState.saturation(oilPhaseIdx);
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
        {
            assert(pvIdx >= x00PvIdx && pvIdx < x00PvIdx+numMisciblePhases*numMiscibleComponents);
            unsigned misciblePhaseIdx = (pvIdx - x00PvIdx)/numMiscibleComponents;
            unsigned miscibleCompIdx = (pvIdx - x00PvIdx)%numMiscibleComponents;
            unsigned phaseIdx;
            switch(misciblePhaseIdx) {
            case 0:
                phaseIdx = oilPhaseIdx;
                break;
            case 1:
                phaseIdx = gasPhaseIdx;
                break;
            default:
                throw std::logic_error("index does not exist, only two miscible phases");
            }
            unsigned compIdx;
            switch(miscibleCompIdx) {
            case 0:
                compIdx = OctaneIdx;
                break;
            case 1:
                compIdx = CO2Idx;
                break;
            default:
                throw std::logic_error("index does not exist, only two miscible components");
            }
            return fluidState.moleFraction(phaseIdx, compIdx);
        }
    }

    // set a quantity in the fluid state
    template <class FluidState>
    static void setQuantity_(FluidState& fluidState,
                             unsigned pvIdx,
                             const typename FluidState::Scalar& value)
    {
        assert(pvIdx < numEq);

        Valgrind::CheckDefined(value);
        // first pressure
        if (pvIdx == p0PvIdx) {
            unsigned phaseIdx = oilPhaseIdx;
            fluidState.setPressure(phaseIdx, value);
        }
        // first (non-water) saturation
        else if (pvIdx == S0PvIdx) {
            unsigned phaseIdx = oilPhaseIdx;
            fluidState.setSaturation(phaseIdx, value);
        }
        // mole fractions
        else {
            assert(pvIdx < numMisciblePhases*numMiscibleComponents + numMisciblePhases);
            unsigned misciblePhaseIdx = (pvIdx - x00PvIdx)/numMiscibleComponents;
            unsigned miscibleCompIdx = (pvIdx - x00PvIdx)%numMiscibleComponents;
            unsigned phaseIdx;
            switch(misciblePhaseIdx) {
            case 0:
                phaseIdx = oilPhaseIdx;
                break;
            case 1:
                phaseIdx = gasPhaseIdx;
                break;
            default:
                throw std::logic_error("index does not exist, only two miscible phases");
            }
            unsigned compIdx;
            switch(miscibleCompIdx) {
            case 0:
                compIdx = OctaneIdx;
                break;
            case 1:
                compIdx = CO2Idx;
                break;
            default:
                throw std::logic_error("index does not exist, only two miscible components");
            }
            fluidState.setMoleFraction(phaseIdx, compIdx, value);
        }
    }

    template <class FluidState>
    static Scalar quantityWeight_(const FluidState& /*fluidState*/, unsigned pvIdx)
    {
        // first pressure
        if (pvIdx == p0PvIdx)
            return 1e-6;
        // first (non-water) saturation
        else if (pvIdx == S0PvIdx)
            return 1.0;
        // mole fractions
        else
            return 1.0;
    }
};

} // namespace Opm

#endif
