
#include "simX/generator/TwoBodyKinematics.h"

#include "ausa/util/memory"
#include "ausa/util/stream.h"

#include <iostream>
#include <simX/Logger.h>

using namespace std;
using namespace simX;


TwoBodyKinematics::TwoBodyKinematics( const Particle &parent, const Particle &d1, const Particle &d2 )
 : parent(parent), d1(d1), d2(d2), pDaughters(2)
{
}

std::vector<TLorentzVector>& TwoBodyKinematics::getFourMomenta() {
    static auto log = log::getLogger("TwoBodyKinematics");

    // Total energy
    double E = parent.getMass() + parent.getExcitationEnergy();

    
    log->debug("Total energy {}", E);

    // Masses of daughters
    double X1 = d1.getExcitationEnergy();
    double X2 = d2.getExcitationEnergy();
    double mx1 = d1.getMass() + X1;
    double mx2 = d2.getMass() + X2;

    // Total energy of daughters
    double E1 = (E * E + mx1 * mx1 - mx2 * mx2) / (2* E);
    double E2 = E - E1; // Remaining energy

    log->debug("mx1: {}, mx2: {}, E1: {}, E2: {}, X1: {}, X2: {}", mx1, mx2, E1, E2, X1, X2);

    // Momentum magnitude of daughters
    double p = sqrt(E1*E1 - mx1*mx1);
    log->debug("p: {}", p);

    // Let daughter 1 travel along positive Z-axis 
    // and daughter 2 in opposite direction.
    // (rotation takes place in TwoBodyDecay.cpp)
    pDaughters[0] = TLorentzVector(0.,0.,p,E1);
    pDaughters[1] = TLorentzVector(0.,0.,-p,E2);

    return pDaughters;
}

