
#include "simX/Particle.h"
#include <vector>
#include <stdexcept>
#include <assert.h>

using namespace std;
using namespace simX;
using namespace AUSA::EnergyLoss;

Particle::Particle(const Ion ion, bool tracking, double ex0, double g0, Particle* parent)
        : ion(ion),
          tracking(tracking),                      // this sets ionizing=true if Z>0 and ionizing=false otherwise.
          position(0.,0.,0.),                // particle created at origo
          direction(0,0,0),                // particle created at origo
          fourMomentumLab(0.,0.,0.,ion.getMass()),     // particle created at rest
          nominalExcitationEnergy(ex0),
          nominalWidth(g0),
          parent(parent),
          minWidth(5.)                                 // minimum width for excitation sampling (keV)
{
    excitationEnergy = nominalExcitationEnergy;
    minExcitationEnergy = nominalExcitationEnergy - 5.0*nominalWidth;
    maxExcitationEnergy = nominalExcitationEnergy + 5.0*nominalWidth;
    // So, by default, we sample excitation energies in the range [E0-5*G0,E0+5*G0]
    // for a breit-wigner with G0=FWHM, this sample the distribution down to 1e-2 of
    // its peak value. However, note that the sampling range actually used in the 
    // simulation covers the entire range allowed by energy conservation.
    // See EventSimulator::EventSimulator and ExcitationSampler::determineLimits for 
    // implementation details. 
}

Particle::Particle( const Ion ion, double ex0, double g0, Particle* parent )
    : Particle(ion, true, ex0, g0, parent)
{

}

void Particle::setFourMomentumLab(const TVector3& p) {
    double mx = ion.getMass() + excitationEnergy;
    double E = sqrt( pow(mx,2.) + p.Mag2() );
    fourMomentumLab.SetVect(p);
    fourMomentumLab.SetE(E);

    direction = p.Unit();
}


void Particle::setFourMomentumLab(double ekin, const TVector3& d) {
    if (d.Mag()==0) throw invalid_argument("Particle direction cannot be null-vector.");
    Particle::direction = d.Unit();

    double mx = ion.getMass() + excitationEnergy;
    double E = mx + ekin;
    double pmag = sqrt( pow(E,2.) - pow(mx,2.) );
    TVector3 p = pmag * direction;
    fourMomentumLab.SetVect(p);
    fourMomentumLab.SetE(E);
}


void Particle::setFourMomentumLab(const TLorentzVector& p4) {
    fourMomentumLab = p4;
    excitationEnergy = fourMomentumLab.Mag() - ion.getMass();
    direction = fourMomentumLab.Vect().Unit();
}


void Particle::updateKineticEnergyLab( double ekin ) {
    double mx = ion.getMass() + excitationEnergy;
    double E = mx+ekin;
    double p = sqrt(pow(E,2)-pow(mx,2));
    TVector3 plab = fourMomentumLab.Vect();
    double p0 = plab.Mag();
    if (p0==0) throw invalid_argument("updateKineticEnergy requires particle to have non-zero momentum.");
    plab *= p/p0;
    fourMomentumLab.SetVect(plab);
    fourMomentumLab.SetE(E);

    direction = fourMomentumLab.Vect().Unit();
}


double Particle::getKineticEnergyLab() const {
    return fourMomentumLab.E()-(ion.getMass()+excitationEnergy);
}


const TVector3& Particle::getDirectionLab() const {
    return direction;
}


TVector3 Particle::getDirectionRFOP() const {
    TVector3 dir;
    if (parent == nullptr)
        dir = getDirectionLab();
    else {
        TLorentzVector v = fourMomentumLab;
        TVector3 b = parent -> getVelocityLab();
        v.Boost(b);
        dir = v.Vect().Unit();
    }
    return dir;
}


bool Particle::hasFiniteWidth() const {
    return nominalWidth>minWidth;
}

void Particle::propagate(const Layer& layer, Particle& part, double range, double * nonIonizing) {
    assert(propagator != nullptr);
    propagator -> propagate(layer, part, range, nonIonizing);
}


