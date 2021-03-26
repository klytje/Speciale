//
// Created by munk on 19-11-15.
//

#include "simX/angular/RutherfordDistribution.h"
#include "simX/Logger.h"
#include <cmath>
#include <iostream>
#include <Math/Factory.h>
#include <boost/math/constants/constants.hpp>
#include <Math/DistSamplerOptions.h>


namespace {
    double ruther(double* xx, double* par) {
        auto x = *xx;
        return std::sin(x)/std::pow(std::sin(x/2), 4);
    }

    namespace {
        double toRadians(double deg) {
            return deg * boost::math::constants::pi<double>()/180;
        }
    }
}

using namespace simX::angular;

RutherfordDistribution::RutherfordDistribution(RutherfordDistribution::Limits phi,
                                                     RutherfordDistribution::Limits theta) {
    if (theta.first == 0) {
        auto logger = log::getLogger("RutherfordDistribution");
        logger->warn("Rutherford diverges at 0 degrees. Incrementing lower limit to 1.");
        theta.first = 1;
    }


    function = TF1("", ruther);

    ROOT::Math::DistSamplerOptions::SetDefaultSampler("Foam");
    sampler = ROOT::Math::Factory::CreateDistSampler();
    sampler -> SetFunction( function, 2 );
    double min[2], max[2];
    min[1] = toRadians(phi.first);
    min[0] = toRadians(theta.first);

    max[1] = toRadians(phi.second);
    max[0] = toRadians(theta.second);

    sampler->SetRange(min, max);
    sampler->Init();
}

void RutherfordDistribution::sampleAngles(double& theta, double& phi) const {
    double v[2];
    sampler->Sample(v);

    theta = v[0];
    phi = v[1];
}


