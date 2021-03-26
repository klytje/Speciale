//
// Created by jesper on 4/4/16.
//

#include "simX/generator/SingleBodyKinematics.h"

using namespace std;
using namespace simX;


SingleBodyKinematics::SingleBodyKinematics( const Particle &particle)
        : particle(particle), momenta(1)
{
}

std::vector<TLorentzVector>& SingleBodyKinematics::getFourMomenta() {
    momenta[0] = particle.getFourMomentumLab();
    return momenta;

}