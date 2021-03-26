
#include "simX/weight/SimpleBreitWigner.h"
#include "simX/Particle.h"
#include "simX/NBodyDecay.h"

using namespace std;
using namespace simX;


SimpleBreitWigner::SimpleBreitWigner( NBodyDecay& proc )
{
    // Set pointer
    parent = proc.getDaughters()[0]->getParent();
    
    // Properties of parent
    isBroad = parent -> hasFiniteWidth();
    EX0     = parent -> getNominalExcitationEnergy();
    G0      = parent -> getNominalWidth();
}

SimpleBreitWigner::~SimpleBreitWigner()
{}

double SimpleBreitWigner::getWeight() const {
    double EX = parent -> getExcitationEnergy();
    double w = 1;
    if (isBroad) w = 1. / ( pow(EX0-EX,2) + pow(G0/2,2) );
    return w;
}

