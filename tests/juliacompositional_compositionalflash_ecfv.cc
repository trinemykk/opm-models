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
 * \brief Box problem with two phases and multiple components
 */
#include "config.h"

#include <opm/models/utils/start.hh>
//#include <opm/models/immiscible/immisciblemodel.hh>
//#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include "problems/oldjuliacompositionalproblem.hh"
//#include "problems/threephasecompositionalproblem.hh"


namespace Opm::Properties {

namespace TTag {
    struct Co2InjectionCompositionalEcfvProblem 
{
    using InheritsFrom = std::tuple<Co2InjectionCompositional, FlashModel>;
};
}

template <class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::Co2InjectionCompositionalEcfvProblem>
{
    using type = TTag::EcfvDiscretization;
};
//TESTAD
template <class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::Co2InjectionCompositionalEcfvProblem>
{
    using type = TTag::AutoDiffLocalLinearizer;
};


} // namespace Opm::Properties

int main(int argc, char **argv)
{
    using EcfvProblemTypeTag = Opm::Properties::TTag::Co2InjectionCompositionalEcfvProblem;
    return Opm::start<EcfvProblemTypeTag>(argc, argv);
}