#include "config.h"
#include "chiwomsproblem.hpp"

#include <opm/models/utils/start.hh>

namespace Opm::Properties {

// NEW_TYPE_TAG(ChiwomsNcpEcfvProblem, INHERITS_FROM(FlashModel, ChiwomsProblem));
// SET_TAG_PROP(ChiwomsNcpEcfvProblem, SpatialDiscretizationSplice, EcfvDiscretization);
// SET_TAG_PROP(ChiwomsNcpEcfvProblem, LocalLinearizerSplice, AutoDiffLocalLinearizer);
namespace TTag {
    struct ChiwomsNcpEcfvProblem {using InheritsFrom = std::tuple<ChiwomsProblem, FlashModel>;};

}
template <class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::ChiwomsNcpEcfvProblem>
{
    using type = TTag::EcfvDiscretization;
};
template <class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::ChiwomsNcpEcfvProblem>
{
    using type = TTag::AutoDiffLocalLinearizer;
};


} // namespace Opm::Properties

int main(int argc, char **argv)
{
    // typedef TTAG(ChiwomsNcpEcfvProblem) EcfvProblemTypeTag;
    using EcfvProblemTypeTag = Opm::Properties::TTag::ChiwomsNcpEcfvProblem;
    return Opm::start<EcfvProblemTypeTag>(argc, argv);
}