//
// Created by jesper on 4/1/16.
//

#ifndef SIMX_ELOSSFACTORY_H
#define SIMX_ELOSSFACTORY_H

#include <simX/Config.h>
#include "ConfigParser.h"
#include "simX/Layer.h"
#include "simX/Particle.h"
#include "simX/propagator/MCStragglingPropagator.h"

namespace simX {
    namespace parser {
        /**
         * Factory that creates a NonIonizingPropagator or a NoTrackPropagator,
         * depending on the particle.
         */
        Config::PropagatorFactory noLossFactory(ConfigParser::Options& o);

        /**
         * Factory that creates an ionizing propagator. If the supplied options contains the
         * tag "tabulation" it will look for either SRIM13 or GEANT. Otherwise it will throw. It
         * defaults to GEANT if the tag is not present.
         */
        Config::PropagatorFactory elossFactory(ConfigParser::Options& o);

        /**
         * Factory that creates a Gaussian straggling propagator. If the supplied options contains the
         * tag "tabulation" it will look for either SRIM13 or GEANT. Otherwise it will throw. It
         * defaults to GEANT if the tag is not present.
         */
        Config::PropagatorFactory gaussianStragglingFactory(ConfigParser::Options& o);

        /**
         * Factory that creates a Monte-Carlo straggling propagator. If the supplied options contains the
         * tag "tabulation" it will look for either SRIM13 or GEANT. Otherwise it will throw. It
         * defaults to GEANT if the tag is not present.
         */
        Config::PropagatorFactory mcStragglingFactory(ConfigParser::Options& o);

        simX::propagator::IonizingPropagator::CalculatorFactory elossCalculator(ConfigParser::Options& o);

        /** 
         * Factory for creating energy-loss calculators (needed by mcStragglingFactory).
         */
        simX::propagator::MCStragglingPropagator::LossCalc lossFactory(const Layer& layer, const Particle& p);

        /** 
         * Factory for creating range calculators (needed by mcStragglingFactory).
         */
        simX::propagator::MCStragglingPropagator::RangeCalc rangeFactory(const Layer& layer, const Particle& p);
    }
}

#endif //SIMX_ELOSSFACTORY_H
