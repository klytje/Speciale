//
// Created by munk on 19-11-15.
//

#ifndef SIMX_RUTHERFORDDISTRIBUTION_H
#define SIMX_RUTHERFORDDISTRIBUTION_H

#include "AngularCorrelation.h"

#include <utility>
#include <Math/DistSampler.h>
#include <TF1.h>

namespace simX {
    namespace angular {
        class RutherfordDistribution : public AngularCorrelation {
        public:
            using Limits = std::pair<double, double>;
            RutherfordDistribution(Limits phi, Limits theta);

            virtual void sampleAngles(double& theta, double& phi) const;

        private:
            ROOT::Math::DistSampler *sampler;
            TF1 function;
        };
    }
}
#endif //SIMX_RUTHERFORDDISTRIBUTION_H
