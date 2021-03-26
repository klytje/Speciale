//
// Created by munk on 07-08-15.
//

#include "simX/propagator/NoTrackPropagator.h"
#include "simX/Particle.h"
#include "simX/Layer.h"

void simX::propagator::NoTrackPropagator::propagate(const simX::Layer& layer, simX::Particle& part, double range, double * nonIonizingEloss) {
    part.setFourMomentumLab(0, part.getDirectionLab());
    
    // All energy loss is assumed to be ionizing
    if (nonIonizingEloss != nullptr) *nonIonizingEloss = 0;
}
