//
// Created by munk on 05-08-15.
//

#include "simX/angular/FixedAngularCorrelation.h"
#include <TMath.h>

using namespace simX::angular;

void FixedAngularCorrelation::sampleAngles(double& theta, double& phi) const {
    theta = FixedAngularCorrelation::theta;
    phi   = FixedAngularCorrelation::phi;
}

FixedAngularCorrelation::FixedAngularCorrelation(double theta, double phi) {
    FixedAngularCorrelation::phi   = phi * TMath::Pi()/180.;
    FixedAngularCorrelation::theta = theta * TMath::Pi()/180.;
}
