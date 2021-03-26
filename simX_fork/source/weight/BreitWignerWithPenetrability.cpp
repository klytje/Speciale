#include "simX/weight/BreitWignerWithPenetrability.h"
#include "simX/TwoBodyDecay.h"
#include "simX/Particle.h"

#include <assert.h>
#include <CoulombFunctions.h>

using namespace std;
using namespace simX;


BreitWignerWithPenetrability::BreitWignerWithPenetrability( NBodyDecay& proc )
 : buffer(1.)
{
    assert(proc.getDaughters().size() == 2);

    // Set pointers
    daughter1 = proc.getDaughters()[0];
    daughter2 = proc.getDaughters()[1];
    parent = daughter1->getParent();
    
    // Initiate coulomb functions
    double AM1 = daughter1 -> getA();
    double AM2 = daughter2 -> getA();
    double Z1  = daughter1 -> getZ();
    double Z2  = daughter2 -> getZ();
    int lorb   = proc.getL();
    double R0  = 1.4;  // reduced channel radius (fm)
    CF = unique_ptr<CoulombFunctions>( new CoulombFunctions(AM1,AM2,Z1,Z2,lorb,R0) );

    // Properties of parent
    isBroad = parent -> hasFiniteWidth();
    EX0     = parent -> getNominalExcitationEnergy();
    G0      = parent -> getNominalWidth();
    
    // Difference in rest mass (not including excitation energies) 
    // between initial and final states
    Q = parent->getMass() - daughter1->getMass() - daughter2->getMass();    
}

BreitWignerWithPenetrability::~BreitWignerWithPenetrability()
{}


double BreitWignerWithPenetrability::getWeight() const {
    // excitation energy of parent
	double EX = parent->getExcitationEnergy();
    // total excitation energy of daughters
    double exsum = daughter1->getExcitationEnergy() + daughter2->getExcitationEnergy();
    // kinetic energy to be shared among daughters
	double E = Q + EX - exsum;		
	// calculate weight
    double w = 0;
    if (E>buffer) {	
        double pen = CF -> penetrability(E);
        if (!isBroad) w = pen;
        else w = pen / ( pow(EX0-EX,2) + pow(G0/2,2) );
	}
    return w;
}    
    //--- TO DO ---//
    // Consider implementing a more sophisticated expression
    // that includes the shift function and takes R-matrix 
    // level energies and reduced widths as inputs
    //--- ---- ---// 


