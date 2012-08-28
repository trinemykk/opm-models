// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Markus Wolff                                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

/*!
 * \ingroup MPFA
 * \ingroup Properties
 * \file
 *
 * \brief Properties for the MPFA-O method.
 */
#ifndef DUMUX_FVMPFAOPROPERTIES2P_HH
#define DUMUX_FVMPFAOPROPERTIES2P_HH

// dumux environment
#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>
#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>

namespace Dumux
{
namespace Properties
{
NEW_TYPE_TAG(FVMPFAOPressureTwoP, INHERITS_FROM(PressureTwoP, MPFAProperties));
}
}

#include <dumux/decoupled/2p/diffusion/fvmpfa/omethod/fvmpfao2pfaboundvelocity2p.hh>

namespace Dumux
{
namespace Properties
{
SET_TYPE_PROP(FVMPFAOPressureTwoP, PressureModel, Dumux::FVMPFAO2PFABoundVelocity2P<TypeTag>);
}
}// end of Dune namespace
#endif