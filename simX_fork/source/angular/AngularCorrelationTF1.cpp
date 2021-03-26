
#include "simX/angular/AngularCorrelationTF1.h"
#include <TMath.h>

using namespace std;
using namespace simX;
using namespace simX::angular;

namespace {
    double degToRad(double deg) {
        return deg*TMath::Pi()/180;
    }

    TF1 makeF1(const string& f, double a, double b) {
        a = degToRad(a);
        b = degToRad(b);
        auto min_ = std::min(a,b);
        auto max_ = std::max(a,b);
        return TF1("", f.c_str(), min_, max_);
    }
}

AngularCorrelationTF1::AngularCorrelationTF1(std::string theta, std::string phi, double phiMin,
                                             double phiMax, double thetaMin, double thetaMax)
    : fPhi(makeF1(phi, phiMin, phiMax)), fTheta(makeF1(theta, thetaMin, thetaMax))
{
}

void AngularCorrelationTF1::sampleAngles( double& theta, double& phi) const {
    phi = fPhi.GetRandom();
    theta = fTheta.GetRandom();
}



