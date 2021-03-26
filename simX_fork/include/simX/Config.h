//
// Created by munk on 07-08-15.
//

#ifndef SIMX_CONFIG_H
#define SIMX_CONFIG_H

#include "Particle.h"
#include "Beam.h"

#include <boost/optional.hpp>

namespace simX {
    class Config {
    public:
        using PropagatorFactory = std::function<Particle::Propagator(Particle&)>;

        Config(PropagatorFactory factoryTarget, PropagatorFactory factoryDetection, PropagatorFactory factoryBeam,
                        boost::optional<std::string> beam,
                        boost::optional<std::string> target, boost::optional<std::string> reaction,
                        boost::optional<std::string> detectionSystem)
                : factoryTarget(factoryTarget), factoryDetection(factoryDetection), factoryBeam(factoryBeam),
                  beam(beam), target(target),
                  reaction(reaction), detectionSystem(detectionSystem) { }

        void applyBeamPropagator(Particle& particle);

        void applyTargetPropagator(Particle& particle);

        void applyDetectionPropagator(Particle& particle);

        bool hasTargetFile() {
            return target.is_initialized();
        }

        std::string& getTargetFile() {
            return target.get();
        }

        bool hasReactionFile() {
            return reaction.is_initialized();
        }

        std::string& getReactionFile() {
            return reaction.get();
        }

        bool hasDetectionSystemFile() {
            return detectionSystem.is_initialized();
        }

        std::string& getDetectionSystemFile() {
            return detectionSystem.get();
        }

        bool hasBeamFile() {
            return beam.is_initialized();
        }

        std::string& getBeamFile() {
            return beam.get();
        }
    private:
        PropagatorFactory factoryTarget, factoryDetection, factoryBeam;
        boost::optional<std::string> beam, target, reaction, detectionSystem;
    };
}

#endif //SIMX_CONFIG_H
