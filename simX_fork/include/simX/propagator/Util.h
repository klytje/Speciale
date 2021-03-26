//
// Created by munk on 06-08-15.
//

#ifndef SIMX_PROPAGATOR_UTIL_H
#define SIMX_PROPAGATOR_UTIL_H

#include "simX/Layer.h"

#include <cmath>


namespace simX {
    class Particle;

    namespace propagator {

        /**
         * Determine maximal range in Layer assuming no straggling or eloss.
         */
        double findMaxRange(const Layer& layer, const Particle& p);

        /**
         * Determine energy-straggling standard deviation based on 
         * empirical relation
         * @param layer Layer of material that particle is propagating in
         * @param part Particle
         * @param dist Distance propagated
         * @param energy Particle energy (typically, taken as the average of initial and final energy)
         */        
        double getEnergyStdDev(const Layer& layer, Particle& part, double dist, double energy);

        /*
         * Computes x^(2/3)
         */
        template <class T>
        T cbrt2(T t) {
          return std::pow(std::cbrt(t), 2);
        }
    }
}
#endif //SIMX_PROPAGATOR_UTIL_H
