//
// Created by munk on 06-08-15.
//

#include "simX/propagator/IonizingPropagator.h"
#include "simX/propagator/Util.h"

#include "simX/Particle.h"
#include "simX/Layer.h"

#include "simX/Random.h"

using namespace std;
using namespace simX;
using namespace simX::propagator;

IonizingPropagator::IonizingPropagator(CalculatorFactory factory)
        : factory(factory)
{

}

void IonizingPropagator::propagate(const Layer& layer, Particle& part, double range, double * nonIonizingEloss) {
    if (range == -1) range = findMaxRange(layer, part);

    // If we have not worked with the Layer (ie. Material) before create a LossCalc, but cache it for next time.
    auto key = make_pair(&layer, &part);
    if (!map.count(key)) map[key] = factory(layer, part);

    auto& calc = map[key];
    auto loss = calc -> getTotalEnergyLoss(part.getKineticEnergyLab(), range, range);

    double E = part.getKineticEnergyLab() - loss;
    if (E < calc->getCutOff()) E = 0;

    part.setFourMomentumLab(E, part.getDirectionLab());

    part.setPosition(part.getPosition()+range*part.getDirectionLab());

    // All energy loss is assumed to be ionizing
    if (nonIonizingEloss != nullptr) *nonIonizingEloss = 0;
}
