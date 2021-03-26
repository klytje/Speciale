//
// Created by munk on 19-11-15.
//

#ifndef SIMX_GRAPHDISTRIBUTION_H
#define SIMX_GRAPHDISTRIBUTION_H

#include "AngularCorrelation.h"

#include <utility>
#include <Math/DistSampler.h>
#include <TF1.h>
#include <memory>
#include <ausa/util/IGenInterpolator.h>

namespace simX {
    namespace angular {
        class GraphDistribution : public AngularCorrelation {
        public:
            GraphDistribution(double phiMin, double phiMax, double thetaMin, double thetaMax, std::string file, std::string format = "%lg %lg");

            virtual void sampleAngles(double& theta, double& phi) const override;

            void printData() const;

        private:
            using Sampler = std::shared_ptr<ROOT::Math::DistSampler>;
            using Func = std::unique_ptr<AUSA::IGenInterpolator>;
            Sampler sampler;
            Func func;

            double phiMin, phiDiff;
            std::string file, format;
        };
    }
}
#endif //SIMX_GRAPHDISTRIBUTION_H
