
#include "simX/weight/TF1WeightCalculator.h"
#include "simX/Particle.h"
#include "simX/NBodyDecay.h"

using namespace std;
using namespace simX;


TF1WeightCalculator::TF1WeightCalculator( NBodyDecay& proc, std::string f)
{
    // Set pointer
    parent = proc.getDaughters()[0]->getParent();
    func = TF1("", f.c_str(), -1e12, 1e12);
}

TF1WeightCalculator::~TF1WeightCalculator()
{}

double TF1WeightCalculator::getWeight() const {
    double ex = parent -> getExcitationEnergy();
    return func.Eval(ex);
}

