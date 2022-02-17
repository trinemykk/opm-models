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
 * \brief Test for the immisicible VCVF discretization with only a single phase
 */
#include "config.h"

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include "problems/smeaheiaProblem.hh"

#include <opm/models/immiscible/immisciblemodel.hh>
//#include <opm/models/flash/flashmodel.hh>
#include <opm/models/io/unstructuredgridvanguard.hh>

#include <opm/models/utils/start.hh>
#include <opm/grid/polyhedralgrid.hh>
//#include <opm/simulators/linalg/umfpackbackend.hh>

namespace Opm::Properties {

// Create new type tags
namespace TTag {
//struct SmeaheiaProblem { using InheritsFrom = std::tuple<SmeaheiaBaseProblem, FlashModel>; };
struct SmeaheiaProblem { using InheritsFrom = std::tuple<SmeaheiaBaseProblem, ImmiscibleModel>; };

} // end namespace TTag

// Set the grid file name for all grids

// template<class TypeTag>
// struct WorldDim<TypeTag, TTag::SmeaheiaProblem> { static constexpr auto value = 3; };
// template<class TypeTag>
// struct DomainDim<TypeTag, TTag::SmeaheiaProblem> { static constexpr auto value = 3; };
// template<class TypeTag>
// struct GridDim<TypeTag, TTag::SmeaheiaProblem> { static constexpr auto value = 3; };

template <class TypeTag>
 struct Grid<TypeTag, TTag::SmeaheiaProblem> {
     using type = Dune::PolyhedralGrid<1, 2>;
 };
//sruct Grid<TypeTag, TTag::SmeaheiaProblem> { using type = Dune::YaspGrid<2>; };

// template<class TypeTag>
// struct Grid<TypeTag, TTag::SmeaheiaProblem> {
// private:
//     //enum{ worldDim = getPropValue<TypeTag, Properties::WorldDim>() };
//     //enum { worldDim = getPropValue<TypeTag, GridView::dimensionworld() };    
// public:
//     using type = Dune::PolyhedralGrid<3, 3>;
// };

template <class TypeTag>
struct Vanguard<TypeTag, TTag::SmeaheiaProblem> { using type = Opm::UnstructuredGridVanguard<TypeTag>; };

template<class TypeTag>
struct GridFile<TypeTag, TTag::SmeaheiaProblem> { static constexpr auto value = "./home/AD.NORCERESEARCH.NO/tmyk/OPM/opm-data/smeaheia/smeaheia_3.txt"; };

// Define time stepping
// The default for the end time of the simulation
template <class TypeTag>
struct EndTime<TypeTag, TTag::SmeaheiaProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1 * 365 * 24 * 60 * 60;
};

// Add output file name
//template<class TypeTag>
//struct CsvFilePath<TypeTag, TTag::SmeaheiaProblem> { static constexpr auto value = "../python/resultsSmeaheia/"; };

// The default for the initial time step size of the simulation
template <class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::SmeaheiaProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1000;
};



template <class TypeTag>
struct PressureDependentPermeability<TypeTag,  TTag::SmeaheiaProblem> {
    static constexpr bool value = false; 
};

template <class TypeTag>
struct InjectionRate<TypeTag,  TTag::SmeaheiaProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 54.0;
};

template <class TypeTag>
struct NewtonTolerance<TypeTag,  TTag::SmeaheiaProblem> {
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e-3;
};
template <class TypeTag>
struct SimulationName<TypeTag, TTag::SmeaheiaProblem> {
    static constexpr auto value = "SmeaheiaSingle3";
};
template <class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::SmeaheiaProblem> {
    static constexpr bool value = true;
};
// template<class TypeTag>
// struct LinearSolverBackend<TypeTag, TTag::SmeaheiaProblem>
// { using type = Opm::Linear::UMFPACKBackend<TypeTag>; };
template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::SmeaheiaProblem>
{ using type = Opm::Linear::ParallelBiCGStabSolverBackend<TypeTag>; };

} // end namespace Opm::Properties

int main(int argc, char **argv)
{
    using SmeaheiaBenchmarkTypeTag = Opm::Properties::TTag::SmeaheiaProblem;
    Opm::start<SmeaheiaBenchmarkTypeTag>(argc, argv);
}
