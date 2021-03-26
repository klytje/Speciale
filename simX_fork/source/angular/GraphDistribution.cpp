//
// Created by munk on 19-11-15.
//

#include "simX/angular/GraphDistribution.h"
#include "simX/Logger.h"
#include "simX/Random.h"

#include <TGraph.h>
#include <Math/Factory.h>

#include <ausa/util/memory>
#include <boost/math/constants/constants.hpp>

using namespace simX::angular;
using namespace ROOT::Math;
using namespace boost::math::constants;

namespace {
    double toRadians(double deg) {
        return deg * pi<double>()/180;
    }
}

GraphDistribution::GraphDistribution(double phiMin, double phiMax,
                                     double thetaMin, double thetaMax,
                                     std::string file, std::string format)
    : phiMin(toRadians(phiMin)), phiDiff(toRadians(phiMax)-toRadians(phiMin)),
      file(file), format(format)
{
    TGraph g(file.c_str(), format.c_str());

    auto N = g.GetN();
    auto x = g.GetX();
    auto x0 = x[0];
    auto xN = x[N - 1];

    // Check that if the graph is given in radians or degrees
    if (xN > pi<double>() || x0 > pi<double>()) {
        log::getLogger("GraphDistribution")->notice("Graph {} is givin is degrees. Converting to radians.", file);

        for (int i = 0; i < N; i++) {
            x[i] = toRadians(x[i]);
        }

        x0 = x[0];
        xN = x[N - 1];
    }

    // Do not exterpolate beyound the bounds
    auto x0User = toRadians(thetaMin);
    auto xNUser = toRadians(thetaMax);
    if (x0User < x0) {
        log::getLogger("GraphDistribution")->warn("User low limit for theta {} smaller than minimum theta {} from {}", x0User, x0, file);
    } else {
        x0 = x0User;
    }

    if (xNUser > xN) {
        log::getLogger("GraphDistribution")->warn("User high limit for theta {} higher than maximum theta {} from {}", xNUser, xN, file);
    } else {
        xN = xNUser;
    }


    // Build interpolator
    auto inter = std::make_shared<Interpolator>(N);
    inter->SetData(N, g.GetX(), g.GetY());
    func = std::make_unique<AUSA::IGenInterpolator>(inter);


    // Construct
    sampler = std::unique_ptr<DistSampler>(Factory::CreateDistSampler("Foam"));
    sampler->SetFunction(*func, 1);
    sampler->SetRange(&x0, &xN);
    sampler->Init();
}

void GraphDistribution::sampleAngles(double& theta, double& phi) const {
    sampler->Sample(&theta);

    phi = phiMin + rnd()*phiDiff;
}

void GraphDistribution::printData() const {
    TGraph g(file.c_str(), format.c_str());
    g.Print();
}
