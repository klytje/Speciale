//
// Created by jesper on 9/18/19.
//

#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include "simX/weight/TF1WeightCalculator.h"
#include "simX/NBodyDecay.h"

using namespace std;
using namespace simX;
using namespace AUSA::EnergyLoss;

SUITE(CustomWeightSamplerTest) {


    TEST(Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST(CheckSinus) {

        std::string sinfunc = "sin(x)";

        Particle parent(Ion("Li10"));
        vector<Particle> daughters = {Particle(Ion("Li9"), 0, 0, &parent), Particle(Ion("n"), 0, 0, &parent)};

        NBodyDecay nBodyDecay(parent, daughters);

        parent.setExcitationEnergy(0.5);
        TF1WeightCalculator calc(nBodyDecay, sinfunc);

        CHECK_CLOSE(sin(0.5), calc.getWeight(), 1e-9);
    }

}
