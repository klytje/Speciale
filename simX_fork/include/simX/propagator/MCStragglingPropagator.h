
#ifndef SIMX_PROPAGATOR_MCSTRAGLINGPROPAGATOR_H
#define	SIMX_PROPAGATOR_MCSTRAGLINGPROPAGATOR_H

#include "simX/propagator/ParticlePropagator.h"
#include "simX/propagator/ScatteringIntegral.h"

#include "ausa/eloss/EnergyLossCalculator.h"
#include "ausa/eloss/RangeCalculator.h"

#include <vector>
#include <memory>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include <TH1.h>
#include <TF1.h>
#include <TVector3.h>

namespace simX {
    class Particle;

    namespace propagator {

        /**
        * Derived class for simulating multiple scattering of ions
        * due to screened Coulomb potential of target atoms.
        */
        class MCStragglingPropagator : public ParticlePropagator {
            public:
                using LossCalc = std::unique_ptr<AUSA::EnergyLoss::EnergyLossCalculator>;
                using LossCalcFactory = std::function<LossCalc(const Layer&, const Particle&)>;

                using RangeCalc = std::unique_ptr<AUSA::EnergyLoss::RangeCalculator>;
                using RangeCalcFactory = std::function<RangeCalc(const Layer&, const Particle&)>;

                MCStragglingPropagator(LossCalcFactory lossFactory, RangeCalcFactory rangeFactory);

                virtual void propagate(const Layer& layer, Particle& part, double range = -1, double * nonIonizingEloss = nullptr) override;
                
                /**
                * Estimate of energy lost to ionizing processes
                */
                double getIonizingEnergyLoss();

                /**
                * Estimate of energy lost to non-ionizing processes
                */
                double getNonIonizingEnergyLoss();
                
                /**
                * Enable/disable straggling for electronic energy-loss.
                * @param status true/false
                */
                void setEnergyStraggling(bool status);
                
                /**
                * Use this method to increase/decrease the total scattering cross section.
                * @param multiplier Multiplicative factor (equal to 1 by default)
                */
                void setCrossSectionMultiplier(double multiplier);

                /**
                * Get radius multiplicative factor
                */
                double getCrossSectionMultiplier();
                
                /**
                * Use this method to modify the low-enery cut-off.
                * MC simulation of scattering process is terminated when ion energy 
                * falls below the cut off.
                * @param cutOff Low-energy cut-off (equal to 10 keV by default)
                */
                void setLowEnergyCutOff(double cutOff); // default is 10 keV

                /**
                * Get low-enery cut-off.
                */
                double getLowEnergyCutOff();

                /**
                * Get number of collisions
                */
                size_t getNumberOfCollisions();

                /**
                * Get mean free path
                */
                double getMeanFreePath();

            private:
                double getEnergyUnitKeV(double Z1, double M1, double Z2, double M2);
                double getScreeningRadius(double Z1, double Z2);

                LossCalcFactory lossFactory;
                RangeCalcFactory rangeFactory;
                using Key = std::pair<const Layer*, const Particle*>;
                std::unordered_map<Key, LossCalc, boost::hash<Key>> lossMap;
                std::unordered_map<Key, RangeCalc, boost::hash<Key>> rangeMap;
                std::vector<Key> recoilKeys;
                
                double lowECutOff;
                bool energyStraggling;
                double csMultiplier;
                double meanFreePath;
                size_t collisionCounter;

                double ionizingEnergyLoss, nonIonizingEnergyLoss;
            
                std::unique_ptr<ScatteringIntegral> Jz;
                
                TH1F * hn;
        };
    }
}

#endif	/* SIMX_PROPAGATOR_ANGULARSTRAGLINGPROPAGATOR_H */

