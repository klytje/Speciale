//
// Created by munk on 29-06-15.
//

#include "simX/angular/IsotropicAngularCorrelation.h"
#include "simX/Random.h"

#include <boost/math/constants/constants.hpp>
#include <cmath>

using namespace std;
using namespace simX::angular;

namespace {
    double toRadians(double deg) {
        return deg * boost::math::constants::pi<double>()/180;
    }
}

IsotropicAngularCorrelation::IsotropicAngularCorrelation(const Limits& phi, const Limits& theta)
    : IsotropicAngularCorrelation(phi.first, phi.second, theta.first, theta.second)
{ }

IsotropicAngularCorrelation::IsotropicAngularCorrelation()
    : IsotropicAngularCorrelation(0, 360, 0, 180)
{

}

IsotropicAngularCorrelation::IsotropicAngularCorrelation(degree phiMin, degree phiMax, degree thetaMin,
                                                               degree thetaMax) :
    phiLow(toRadians(std::min(phiMin, phiMax))), phiDiff(toRadians(std::abs(phiMax - phiMin)))
{
    double a = cos(toRadians(thetaMin));
    double b = cos(toRadians(thetaMax));

    thetaLow = std::min(a, b);
    thetaDiff = std::abs(a - b);
}

void IsotropicAngularCorrelation::sampleAngles(double& theta, double& phi) const {
    phi = phiLow + rnd() * phiDiff;
    theta = std::acos(thetaLow + rnd() * thetaDiff);
}



