//
// Created by munk on 06-08-15.
//

#ifndef SIMX_IONIZINGPROPAGATOR_H
#define SIMX_IONIZINGPROPAGATOR_H

#include "ParticlePropagator.h"

#include <memory>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include <ausa/eloss/EnergyLossCalculator.h>

namespace simX {
    namespace propagator {
        /**
        * Derived class for propagating ionizing particles
        */
        class IonizingPropagator : public ParticlePropagator {
        public:
            using LossCalc = std::unique_ptr<AUSA::EnergyLoss::EnergyLossCalculator>;
            using CalculatorFactory = std::function<LossCalc(const Layer&, const Particle&)>;

            IonizingPropagator(CalculatorFactory factory);

            virtual void propagate(const Layer& layer, Particle& part, double range = -1, double * nonIonizingEloss = nullptr) override;

        private:
            CalculatorFactory factory;
            using Key = std::pair<const Layer*, const Particle*>;
            std::unordered_map<Key, LossCalc, boost::hash<Key>> map;
        };
    }
}


#endif //SIMX_IONIZINGPROPAGATOR_H
