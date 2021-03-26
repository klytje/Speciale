#include "simX/weight/GammaPhaseSpace.h"
#include "simX/TwoBodyDecay.h"
#include "simX/Particle.h"

#include "math.h"
#include <assert.h>

using namespace std;
using namespace simX;


GammaPhaseSpace::GammaPhaseSpace(NBodyDecay& proc )
 : lorb(proc.getL())
{
    assert(proc.getDaughters().size() == 2);

    // Set pointers
    daughter1 = proc.getDaughters()[0];
    daughter2 = proc.getDaughters()[1];
    parent = daughter1->getParent();

   // Properties of parent
    isBroad = parent -> hasFiniteWidth();
    EX0     = parent -> getNominalExcitationEnergy();
    G0      = parent -> getNominalWidth();
    
    // Difference in rest mass (not including excitation energies) 
    // between initial and final states
    Q = parent->getMass() - daughter1->getMass() - daughter2->getMass();    
}

GammaPhaseSpace::~GammaPhaseSpace()
{}


double GammaPhaseSpace::getWeight() const {
    // excitation energy of parent
	double EX = parent->getExcitationEnergy();
    // total excitation energy of daughters
    double exsum = daughter1->getExcitationEnergy() + daughter2->getExcitationEnergy();
    // kinetic energy to be shared among daughters
	double E = Q + EX - exsum;		
	// calculate weight
    double a = 2.*lorb+1.;
    double f = 0;
	if (E>0.) f = pow( E/1000., a ); // gamma phase-space factor is proportional to E^2l+1
	return f;
}

