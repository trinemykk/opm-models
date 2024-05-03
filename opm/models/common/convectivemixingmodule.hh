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
 * \brief Classes required for convective mixing.
 */
#ifndef EWOMS_CONVECTIVEMIXING_MODULE_HH
#define EWOMS_CONVECTIVEMIXING_MODULE_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

namespace Opm {

/*!
 * \ingroup Diffusion
 * \class Opm::ConvectiveMixingModule
 * \brief P
 */
template <class TypeTag>
class ConvectiveMixingModule;

/*!
 * \copydoc Opm::ConvectiveMixingModule
 */
template <class TypeTag>
class ConvectiveMixingModule<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;

public:
    /*!
     * \brief Register all run-time parameters for the convective mixing module.
     */
    static void registerParameters()
    {}

    /*!
     * \brief Adds the diffusive mass flux flux to the flux vector over a flux
     *        integration point.
      */
    template <class Context>
    static void addConvectiveMixingFlux(RateVector&,
                                 const Context&,
                                 unsigned,
                                 unsigned)
    {}
};
}