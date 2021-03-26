
#include <unittest++/UnitTest++.h>

#include "simX/angular/GraphDistribution.h"

#include <math.h>
#include <TMath.h>
#include <tuple>
#include <iostream>

using namespace simX;
using namespace simX::angular;

SUITE(GraphDistributionTest) {

    TEST(LoadingDWBAFileFromNickWorks) {
        GraphDistribution g(0, 360, 1, 160, "test/_res/dwba_output.dat");

        double theta, phi;
        g.sampleAngles(theta, phi);
    }

    TEST(DistWillAutomaticallyClampToMinMaxValuesFromFile) {
        GraphDistribution(0, 360, 0, 160, "test/_res/dwba_output.dat");
    }
}
