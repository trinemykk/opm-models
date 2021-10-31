#include "config.h"
#include "chiwomsproblem.hpp"

#include <opm/models/utils/start.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(ChiwomsNcpEcfvProblem, INHERITS_FROM(FlashModel, ChiwomsProblem));
SET_TAG_PROP(ChiwomsNcpEcfvProblem, SpatialDiscretizationSplice, EcfvDiscretization);
SET_TAG_PROP(ChiwomsNcpEcfvProblem, LocalLinearizerSplice, AutoDiffLocalLinearizer);

END_PROPERTIES

int main(int argc, char **argv)
{
    typedef TTAG(ChiwomsNcpEcfvProblem) EcfvProblemTypeTag;
    return Opm::start<EcfvProblemTypeTag>(argc, argv);
}
