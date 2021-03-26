
#ifndef NONIONIZINGPROPAGATOR_H
#define	NONIONIZINGPROPAGATOR_H

#include "ParticlePropagator.h"

namespace simX {
    namespace propagator {
        /**
        * Derived class for propagating non-ionizing particles
        * such as gamma rays and neutrons.
        */
        class NonIonizingPropagator : public ParticlePropagator {

        public:
            NonIonizingPropagator() = default;

            virtual void propagate(const Layer& layer, Particle& part, double range = -1, double * nonIonizingEloss = nullptr) override;
        };
    }
}

#endif	/* NONIONIZINGPROPAGATOR_H */

