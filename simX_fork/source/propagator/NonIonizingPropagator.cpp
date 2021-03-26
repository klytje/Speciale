#include "simX/propagator/NonIonizingPropagator.h"
#include "simX/propagator/Util.h"

#include "simX/Layer.h"
#include "simX/Particle.h"

using namespace simX;
using namespace simX::propagator;

void NonIonizingPropagator::propagate(const Layer& layer, Particle& part, double range, double * nonIonizingEloss) {
    if (range == -1) range = findMaxRange(layer, part);

    auto& pos = part.getPosition();
    auto& dir = part.getDirectionLab();

    part.setPosition(pos+range*dir);

    // All energy loss is assumed to be ionizing
    if (nonIonizingEloss != nullptr) *nonIonizingEloss = 0;
}
