//
// Created by munk on 28-06-15.
//

#include "simX/CompoundFormation.h"
#include <assert.h>


using namespace std;
using namespace simX;

namespace {
    Particle makeProduct(Particle& p1, Particle& p2) {
        using AUSA::EnergyLoss::Ion;
        auto A = p1.getA() + p2.getA();
        auto Z = p1.getZ() + p2.getZ();
        // nominal excitation energy = m1+m2-m(1+2)
        auto Ex0 = p1.getMass()+p2.getMass()-Ion(Z,A).getMass(); 
        return Particle(Ion(Z, A), Ex0, 0, &p1);  // compound has p1 (beam) as parent
    }
}

CompoundFormation::CompoundFormation(Particle& beam, Particle target)
    : beam(beam), target(target), product(makeProduct(beam, target))
{
    daughters.push_back(&product);
}

void CompoundFormation::runProcess() {
    product.setFourMomentumLab( beam.getFourMomentumLab() + target.getFourMomentumLab() );
    // excitation energy is automatically set by call to setFourMomentumLab( TLorentzVector )
    product.setPosition( beam.getPosition() );
}

std::vector<Particle*>& CompoundFormation::getDaughters() {
    return daughters;
}

const std::vector<Particle *>& CompoundFormation::getDaughters() const {
    return daughters;
}

void CompoundFormation::setWeightCalculator( unique_ptr<WeightCalculator> wCalc ) {
    weightCalculator = move(wCalc); 
}

double CompoundFormation::getWeight() const { 
//    assert(weightCalculator!=nullptr);
    double w = 1;
    if (weightCalculator!=nullptr) w = weightCalculator -> getWeight();
    return w;
}

int CompoundFormation::getNumberOfNucleons() const {
    return beam.getA()+target.getA();
}

int CompoundFormation::getNumberOfProtonsInInitialState() const {
    return beam.getZ()+target.getZ();
}



