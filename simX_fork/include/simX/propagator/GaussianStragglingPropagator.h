
#ifndef SIMX_PROPAGATOR_ANGULARSTRAGLINGPROPAGATOR_H
#define	SIMX_PROPAGATOR_ANGULARSTRAGLINGPROPAGATOR_H

#include "ParticlePropagator.h"

#include <TF1.h>
#include <TVector3.h>

#include <vector>
#include <memory>

namespace simX {
    class Particle;

    namespace propagator {
        /**
        * Derived class for simulating angular straggling of ions
        */
        class GaussianStragglingPropagator : public ParticlePropagator {
        public:

            enum class SCREENING {
                THOMAS_FERMI = 0,
                LENZ_JENSEN = 1
            };

            GaussianStragglingPropagator(std::shared_ptr<ParticlePropagator> inner, 
                                         bool angularStraggling = true,
                                         bool energyStraggling = false);

            virtual void propagate(const Layer& layer, Particle& part, double range = -1, double * nonIonizingEloss = nullptr) override;
            
            /**
            * Set screening model (only affects angular straggling for tau < 0.1)
            * @param screening Screening model (THOMAS_FERMI or LENZ_JENSEN)
            */
            void setScreeningModel(SCREENING screening);

            /**
            * Returns most recently computed half-scattering angle in degrees
            * (mostly for testing purposes)
            */
            double getHalfScatteringAngle() const;

        private:
            double getHalfScatteringAngle(const Layer& layer, Particle& part, double dist, double energy);
            void doAngularStraggling(const Layer& layer, Particle& part, double dist, double ei, double ef);

//            double getEnergyStdDev(const Layer& layer, Particle& part, double dist, double energy);
            void doEnergyStraggling(const Layer& layer, Particle& part, double dist, double ei, double ef);

            double getAlphaTilde(double tau);

            std::vector<double>& getInterpolationFactors(double tau);
            
            SCREENING screen;
            std::shared_ptr<ParticlePropagator> inner;
            TF1 scattFunc;
            double halfScattAngle;
            double nsig;
            
            bool angularStraggling;
            bool energyStraggling;
            
            std::vector<double> ws;
        };
    }
}

#endif	/* SIMX_PROPAGATOR_ANGULARSTRAGLINGPROPAGATOR_H */

