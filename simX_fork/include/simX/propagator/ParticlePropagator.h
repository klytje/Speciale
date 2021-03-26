#ifndef SIMX_PARTICLEPROPAGATOR_H
#define	SIMX_PARTICLEPROPAGATOR_H

#include "simX/Layer.h"

namespace simX {
    class Particle;

    namespace propagator {
        /**
        * Base class for propagation of particles in materials.
        */
        class ParticlePropagator {
        public:
            ParticlePropagator() = default;

            virtual ~ParticlePropagator() = default;

            /**
            * Propagates particle through layer until it:
            * 1) exists
            * 2) has lost all its energy
            * 3) decays (if it has a finite lifetime)
            * Alternatively, the user can specify a fixed range.
            * @param layer Layer that particle is propagating in.
            * @param part Particle.
            * @param range Range of propagation (optional).
            * @param nonIonizingEloss Non-ionizing energy loss.
             */
            virtual void propagate(const Layer& layer, Particle& part, double range = -1, double * nonIonizingEloss = nullptr) = 0;
        };
    }
}

#endif	/* SIMX_PARTICLEPROPAGATOR_H */

