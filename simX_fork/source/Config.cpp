//
// Created by munk on 07-08-15.
//

#include "simX/Config.h"
#include "simX/propagator/IonizingPropagator.h"
#include "simX/Layer.h"

#include <ausa/util/Resource.h>
#include <ausa/eloss/SRIMTabulation.h>
#include <ausa/eloss/RangeInterpolator.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>
#include <ausa/util/memory>

using namespace simX;
using namespace simX::propagator;
using namespace AUSA::EnergyLoss;

void Config::applyTargetPropagator(simX::Particle& particle) {
    particle.setPropagator(factoryTarget(particle));
}

void Config::applyDetectionPropagator(Particle& particle) {
    particle.setPropagator(factoryDetection(particle));
}

void Config::applyBeamPropagator(Particle &particle) {
    particle.setPropagator(factoryBeam(particle));
}
