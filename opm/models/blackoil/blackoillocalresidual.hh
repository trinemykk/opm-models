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
 *
 * \copydoc Opm::BlackOilLocalResidual
 */
#ifndef EWOMS_BLACK_OIL_LOCAL_RESIDUAL_HH
#define EWOMS_BLACK_OIL_LOCAL_RESIDUAL_HH

#include "blackoilproperties.hh"
#include "blackoilsolventmodules.hh"
#include "blackoilextbomodules.hh"
#include "blackoilpolymermodules.hh"
#include "blackoilenergymodules.hh"
#include "blackoilfoammodules.hh"
#include "blackoilbrinemodules.hh"
#include "blackoildiffusionmodule.hh"
#include "blackoilmicpmodules.hh"
#include <opm/material/fluidstates/BlackOilFluidState.hpp>

namespace Opm {
/*!
 * \ingroup BlackOilModel
 *
 * \brief Calculates the local residual of the black oil model.
 */
template <class TypeTag>
class BlackOilLocalResidual : public GetPropType<TypeTag, Properties::DiscLocalResidual>
{
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };

    static const bool waterEnabled = Indices::waterEnabled;
    static const bool gasEnabled = Indices::gasEnabled;
    static const bool oilEnabled = Indices::oilEnabled;
    static const bool compositionSwitchEnabled = (compositionSwitchIdx >= 0);

    static constexpr bool blackoilConserveSurfaceVolume = getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>();
    static constexpr bool enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>();
    static constexpr bool enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>();

    using Toolbox = MathToolbox<Evaluation>;
    using SolventModule = BlackOilSolventModule<TypeTag>;
    using ExtboModule = BlackOilExtboModule<TypeTag>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using EnergyModule = BlackOilEnergyModule<TypeTag>;
    using FoamModule = BlackOilFoamModule<TypeTag>;
    using BrineModule = BlackOilBrineModule<TypeTag>;
    using DiffusionModule = BlackOilDiffusionModule<TypeTag, enableDiffusion>;
    using MICPModule = BlackOilMICPModule<TypeTag>;

public:
    /*!
     * \copydoc FvBaseLocalResidual::computeStorage
     */
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext& elemCtx,
                        unsigned dofIdx,
                        unsigned timeIdx) const
    {
        // retrieve the intensive quantities for the SCV at the specified point in time
        const IntensiveQuantities& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        storage = 0.0;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                if (Indices::numPhases == 3) { // add trivial equation for the pseudo phase
                    unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                    if (timeIdx == 0)
                        storage[conti0EqIdx + activeCompIdx] = variable<LhsEval>(0.0, conti0EqIdx + activeCompIdx);
                    else
                        storage[conti0EqIdx + activeCompIdx] = 0.0;
                }
                continue;
            }

            unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            LhsEval surfaceVolume =
                Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx))
                * Toolbox::template decay<LhsEval>(fs.invB(phaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());

            storage[conti0EqIdx + activeCompIdx] += surfaceVolume;

            // account for dissolved gas
            if (phaseIdx == oilPhaseIdx && FluidSystem::enableDissolvedGas()) {
                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                storage[conti0EqIdx + activeGasCompIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.fluidState().Rs())
                    * surfaceVolume;
            }

            // account for dissolved gas in water phase
            if (phaseIdx == waterPhaseIdx && FluidSystem::enableDissolvedGasInWater()) {
                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                storage[conti0EqIdx + activeGasCompIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.fluidState().Rsw())
                    * surfaceVolume;
            }

            // account for vaporized oil
            if (phaseIdx == gasPhaseIdx && FluidSystem::enableVaporizedOil()) {
                unsigned activeOilCompIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
                storage[conti0EqIdx + activeOilCompIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.fluidState().Rv())
                    * surfaceVolume;
            }

            // account for vaporized water
            if (phaseIdx == gasPhaseIdx && FluidSystem::enableVaporizedWater()) {
                unsigned activeWaterCompIdx = Indices::canonicalToActiveComponentIndex(waterCompIdx);
                storage[conti0EqIdx + activeWaterCompIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.fluidState().Rvw())
                    * surfaceVolume;
            }
        }

        adaptMassConservationQuantities_(storage, intQuants.pvtRegionIndex());

        // deal with solvents (if present)
        SolventModule::addStorage(storage, intQuants);

        // deal with zFracton (if present)
        ExtboModule::addStorage(storage, intQuants);

        // deal with polymer (if present)
        PolymerModule::addStorage(storage, intQuants);

        // deal with energy (if present)
        EnergyModule::addStorage(storage, intQuants);

        // deal with foam (if present)
        FoamModule::addStorage(storage, intQuants);

        // deal with salt (if present)
        BrineModule::addStorage(storage, intQuants);

        // deal with micp (if present)
        MICPModule::addStorage(storage, intQuants);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector& flux,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx) const
    {
        assert(timeIdx == 0);

        flux = 0.0;

        const ExtensiveQuantities& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
            unsigned pvtRegionIdx = up.pvtRegionIndex();
            if (upIdx == focusDofIdx)
                evalPhaseFluxes_<Evaluation>(flux, phaseIdx, pvtRegionIdx, extQuants, up.fluidState());
            else
                evalPhaseFluxes_<Scalar>(flux, phaseIdx, pvtRegionIdx, extQuants, up.fluidState());
        }
		
		const auto& problem = elemCtx.problem();
		const auto& oilVaporizationControl = problem.simulator().vanguard().schedule()[problem.episodeIndex()].oilvap();

		if(oilVaporizationControl.drsdtConvective()) {
			const auto& stencil = elemCtx.stencil(timeIdx);
			const auto& scvf = stencil.interiorFace(scvfIdx);
			Scalar zIn = problem.dofCenterDepth(elemCtx, scvf.interiorIndex(), timeIdx);
			Scalar zEx = problem.dofCenterDepth(elemCtx, scvf.exteriorIndex(), timeIdx);
			// the distances from the DOF's depths. (i.e., the additional depth of the
			// exterior DOF)
			Scalar distZ = zIn - zEx;
			Scalar g = problem.gravity()[2];
			const auto& intQuantsIn = elemCtx.intensiveQuantities(scvf.interiorIndex(), timeIdx);
			const auto& intQuantsEx = elemCtx.intensiveQuantities(scvf.exteriorIndex(), timeIdx);
			const auto& rho_sat_in = FluidSystem::saturatedDensity(intQuantsIn.fluidState(), oilPhaseIdx, intQuantsIn.pvtRegionIndex());
			const auto& rho_sat_ex = Opm::getValue(FluidSystem::saturatedDensity(intQuantsEx.fluidState(), oilPhaseIdx, intQuantsEx.pvtRegionIndex()));
			const auto& rho_in = intQuantsIn.fluidState().density(oilPhaseIdx);
			const auto& rho_ex = Opm::getValue(intQuantsEx.fluidState().density(oilPhaseIdx));
			const auto delta_rho = (rho_sat_ex + rho_sat_in - rho_in -rho_ex)/2;		
			const auto& pressureDiff = extQuants.pressureDifference(oilPhaseIdx);
            const auto pressure_difference_convective_mixing =  delta_rho * g * distZ;
			if (Opm::abs(pressure_difference_convective_mixing) > 1e-12){ // 1e-9
		
				// find new upstream direction
				unsigned upIdx = scvf.interiorIndex();
				unsigned downIdx = scvf.exteriorIndex();
				if (pressure_difference_convective_mixing > 0) {
					upIdx = scvf.exteriorIndex();
					downIdx = scvf.interiorIndex();
				}
                Scalar trans = problem.transmissibility(elemCtx, scvf.interiorIndex(), scvf.exteriorIndex());
                Scalar faceArea = scvf.area();
		
				const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
				const IntensiveQuantities& down = elemCtx.intensiveQuantities(downIdx, timeIdx);
				const auto& Rs =  up.fluidState().Rs();	
				const Evaluation SoMax = 0.0;
				const auto& RsSat = FluidSystem::saturatedDissolutionFactor(up.fluidState(),
                                                            oilPhaseIdx,
                                                            up.pvtRegionIndex(),
                                                            SoMax);
				
                const Evaluation& transMult = up.rockCompTransMultiplier();												
				Evaluation Sm = Opm::min(0.999, Opm::max(0.001, Rs/RsSat));
				const Scalar Xhi = oilVaporizationControl.getMaxDRSDT(intQuantsIn.pvtRegionIndex());
				Scalar Smo = 0.34;
				Evaluation sg = up.fluidState().saturation(FluidSystem::gasPhaseIdx);
				Evaluation S = (Rs - RsSat * sg) / (RsSat * ( 1.0 - sg));
				 
				if ( (S > Smo || down.fluidState().Rs() > 0) ) {
					
				    const auto& invB = up.fluidState().invB(oilPhaseIdx);
                    const auto& visc = up.fluidState().viscosity(oilPhaseIdx);

				    // what will be the flux when muliplied with trans_mob
				    const auto convectiveFlux = -trans*transMult*Xhi*invB*g*distZ*delta_rho*Rs/(visc*faceArea); 
				    unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
				    // Since the upwind direction may have changed from what was used to compute the mobility. 
				    // We keep the derivative of the trans_mob
				if (upIdx == focusDofIdx)
					flux[conti0EqIdx + activeGasCompIdx] += convectiveFlux;
				else 
					flux[conti0EqIdx + activeGasCompIdx] += Opm::getValue(convectiveFlux);
				
				}
        
			}
		}	
        // deal with solvents (if present)
        SolventModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with zFracton (if present)
        ExtboModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with polymer (if present)
        PolymerModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with energy (if present)
        EnergyModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with foam (if present)
        FoamModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with salt (if present)
        BrineModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with micp (if present)
        MICPModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        DiffusionModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeSource
     */
    void computeSource(RateVector& source,
                       const ElementContext& elemCtx,
                       unsigned dofIdx,
                       unsigned timeIdx) const
    {
        // retrieve the source term intrinsic to the problem
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);

        // deal with MICP (if present)
        MICPModule::addSource(source, elemCtx, dofIdx, timeIdx);

        // scale the source term of the energy equation
        if (enableEnergy)
            source[Indices::contiEnergyEqIdx] *= getPropValue<TypeTag, Properties::BlackOilEnergyScalingFactor>();
    }

    /*!
     * \brief Helper function to calculate the flux of mass in terms of conservation
     *        quantities via specific fluid phase over a face.
     */
    template <class UpEval, class FluidState>
    static void evalPhaseFluxes_(RateVector& flux,
                                 unsigned phaseIdx,
                                 unsigned pvtRegionIdx,
                                 const ExtensiveQuantities& extQuants,
                                 const FluidState& upFs)
    {
        const auto& invB = getInvB_<FluidSystem, FluidState, UpEval>(upFs, phaseIdx, pvtRegionIdx);
        const auto& surfaceVolumeFlux = invB*extQuants.volumeFlux(phaseIdx);
        unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));

        if (blackoilConserveSurfaceVolume)
            flux[conti0EqIdx + activeCompIdx] += surfaceVolumeFlux;
        else
            flux[conti0EqIdx + activeCompIdx] += surfaceVolumeFlux*FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);

        if (phaseIdx == oilPhaseIdx) {
            // dissolved gas (in the oil phase).
            if (FluidSystem::enableDissolvedGas()) {
                const auto& Rs = BlackOil::getRs_<FluidSystem, FluidState, UpEval>(upFs, pvtRegionIdx);

                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                if (blackoilConserveSurfaceVolume)
                    flux[conti0EqIdx + activeGasCompIdx] += Rs*surfaceVolumeFlux;
                else
                    flux[conti0EqIdx + activeGasCompIdx] += Rs*surfaceVolumeFlux*FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
        } else if (phaseIdx == waterPhaseIdx) {
            // dissolved gas (in the water phase).
            if (FluidSystem::enableDissolvedGasInWater()) {
                const auto& Rsw = BlackOil::getRsw_<FluidSystem, FluidState, UpEval>(upFs, pvtRegionIdx);

                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                if (blackoilConserveSurfaceVolume)
                    flux[conti0EqIdx + activeGasCompIdx] += Rsw*surfaceVolumeFlux;
                else
                    flux[conti0EqIdx + activeGasCompIdx] += Rsw*surfaceVolumeFlux*FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
        }
        else if (phaseIdx == gasPhaseIdx) {
            // vaporized oil (in the gas phase).
            if (FluidSystem::enableVaporizedOil()) {
                const auto& Rv = BlackOil::getRv_<FluidSystem, FluidState, UpEval>(upFs, pvtRegionIdx);

                unsigned activeOilCompIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
                if (blackoilConserveSurfaceVolume)
                    flux[conti0EqIdx + activeOilCompIdx] += Rv*surfaceVolumeFlux;
                else
                    flux[conti0EqIdx + activeOilCompIdx] += Rv*surfaceVolumeFlux*FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            }
             // vaporized water (in the gas phase).
            if (FluidSystem::enableVaporizedWater()) {
                const auto& Rvw = BlackOil::getRvw_<FluidSystem, FluidState, UpEval>(upFs, pvtRegionIdx);

                unsigned activeWaterCompIdx = Indices::canonicalToActiveComponentIndex(waterCompIdx);
                if (blackoilConserveSurfaceVolume)
                    flux[conti0EqIdx + activeWaterCompIdx] += Rvw*surfaceVolumeFlux;
                else
                    flux[conti0EqIdx + activeWaterCompIdx] += Rvw*surfaceVolumeFlux*FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            }
        }
    }

    /*!
     * \brief Helper function to convert the mass-related parts of a Dune::FieldVector
     *        that stores conservation quantities in terms of "surface-volume" to the
     *        conservation quantities used by the model.
     *
     * Depending on the value of the BlackoilConserveSurfaceVolume property, the model
     * either conserves mass by means of "surface volume" of the components or mass
     * directly. In the former case, this method is a no-op; in the latter, the values
     * passed are multiplied by their respective pure component's density at surface
     * conditions.
     */
    template <class Scalar>
    static void adaptMassConservationQuantities_(Dune::FieldVector<Scalar, numEq>& container, unsigned pvtRegionIdx)
    {
        if (blackoilConserveSurfaceVolume)
            return;

        // convert "surface volume" to mass. this is complicated a bit by the fact that
        // not all phases are necessarily enabled. (we here assume that if a fluid phase
        // is disabled, its respective "main" component is not considered as well.)

        if (waterEnabled) {
            unsigned activeWaterCompIdx = Indices::canonicalToActiveComponentIndex(waterCompIdx);
            container[conti0EqIdx + activeWaterCompIdx] *=
                FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
        }

        if (gasEnabled) {
            unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
            container[conti0EqIdx + activeGasCompIdx] *=
                FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
        }

        if (oilEnabled) {
            unsigned activeOilCompIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
            container[conti0EqIdx + activeOilCompIdx] *=
                FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
        }
    }
};

} // namespace Opm

#endif
