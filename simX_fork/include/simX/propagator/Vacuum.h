//
// Created by munk on 06-08-15.
//

#ifndef SIMX_VACUUM_H
#define SIMX_VACUUM_H

#include "simX/Layer.h"

#include <vector>
#include <stddef.h>

namespace simX {
    class Particle;

    namespace propagator {
        /**
         * Propagates particle in vacuum until it encounters next layer of material.
         * @param layers The layers that the particle can hit.
         * @param part Particle.
         * @return Layer that particle hits next or -1.
         */
        size_t propagateInVacuum(const std::vector<const Layer*>& layers, Particle& part);

        /**
         * Propagates particle in vacuum until it encounters this layer of material.
         * @param layer The layer it will hit
         * @param part The particle to have its position updated.
         *
         * @pre The Particle will hit this particle.
         */
        void propagateInVacuum(const Layer& layer, Particle& part);
    }
}
#endif //SIMX_VACUUM_H
