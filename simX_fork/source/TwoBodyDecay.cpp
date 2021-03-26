#include "simX/TwoBodyDecay.h"
#include <simX/generator/TwoBodyKinematics.h>

#include "ausa/util/memory"

using namespace std;
using namespace simX;
using namespace simX::angular;


TwoBodyDecay::TwoBodyDecay( Particle& parent, Particle d1, Particle d2, unique_ptr<AngularCorrelation> ac, int lorb, unique_ptr<WeightCalculator> ampCalc )
 :  NBodyDecay(parent, {d1, d2}, nullptr, move(ac), move(ampCalc), lorb)
{
    auto& d = getDaughters();
    setFinalStateGenerator(std::make_unique<TwoBodyKinematics>(parent, *d[0], *d[1]));
}


