#include "config.h"
#include "compositionalproblem.hpp"

#include <opm/models/utils/start.hh>

namespace Opm::Properties {

namespace TTag {
    struct CompositionalNcpEcfvProblem {using InheritsFrom = std::tuple<CompositionalProblem, FlashModel>;};

}
template <class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::CompositionalNcpEcfvProblem>
{
    using type = TTag::EcfvDiscretization;
};
template <class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::CompositionalNcpEcfvProblem>
{
    using type = TTag::AutoDiffLocalLinearizer;
};


} // namespace Opm::Properties

int main(int argc, char **argv)
{
    // typedef TTAG(CompositionalNcpEcfvProblem) EcfvProblemTypeTag;
    using EcfvProblemTypeTag = Opm::Properties::TTag::CompositionalNcpEcfvProblem;
    return Opm::start<EcfvProblemTypeTag>(argc, argv);
}
