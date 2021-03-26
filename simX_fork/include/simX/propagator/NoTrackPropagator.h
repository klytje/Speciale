//
// Created by munk on 07-08-15.
//

#ifndef SIMX_NOTRACKPROPAGATOR_H
#define SIMX_NOTRACKPROPAGATOR_H

#include "ParticlePropagator.h"

namespace simX {

    namespace propagator {
        /**
        * Derived class for propagating particles that are not tracked.
        * This will remove all kinetic energy but not change position.
        */
        class NoTrackPropagator : public ParticlePropagator {
        public:
            NoTrackPropagator() = default;

            virtual void propagate(const Layer& layer, Particle& part, double range = -1, double * nonIonizingEloss = nullptr) override;
        };
    }
}


#endif //SIMX_NOTRACKPROPAGATOR_H
