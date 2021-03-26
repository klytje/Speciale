
#ifndef PARTICLE_h
#define PARTICLE_h 1

#include "simX/propagator/ParticlePropagator.h"

// C++ and ROOT header files
#include <string>
#include <memory>
#include <TVector3.h>
#include <TLorentzVector.h>

// AUSAlib header files
#include <ausa/eloss/Ion.h>



namespace simX {

    /**
    * Particle class.
    * @param ion Particle type.
    * @param nominalExcitationEnergy Nominal excitation energy (keV).
    * @param nominalWidth Nominal width (keV).
    */
    class Particle {

        public:
            Particle( const AUSA::EnergyLoss::Ion ion, double nominalExcitationEnergy=0, double nominalWidth=0, Particle* parent = nullptr );
            Particle( const AUSA::EnergyLoss::Ion ion, bool tracking, double nominalExcitationEnergy=0, double nominalWidth=0, Particle* parent = nullptr );
            ~Particle() = default;


            /**
             * True if the user want this particle to be tracked after creation.
             */
            const bool getTracking() const {
                return tracking;
            }

            /**
            * Get particle name.
            */
            inline std::string getName() const { return ion.getName(); }

            inline int getA() const {return ion.getA();}

            inline int getZ() const {return ion.getZ();}
            
            /**
            * Get rest mass.
            */
            inline double getMass() const { return ion.getMass(); }

            /**
            * Set current excitation energy.
            */
            inline void setExcitationEnergy( double ex ) { excitationEnergy=ex; }

            /**
            * Get current excitation energy.
            */
            inline double getExcitationEnergy() const { return excitationEnergy; }

            /**
             * Set the nominal excitation energy.
             */
            inline void setNominalExcitationEnergy( double ex ) { nominalExcitationEnergy=ex; }

            /**
             * Get the nominal excitation energy.
             */
            inline double getNominalExcitationEnergy() const { return nominalExcitationEnergy;}

            /**
             * Set the minimum possible excitation energy.
             * @param ex Minimum possible excitation energy (keV).
             */
            inline void setMinExcitationEnergy( double ex ) { minExcitationEnergy=ex; }

            /**
             * Set the maximum possible excitation energy.
             * @param ex Maximum possible excitation energy (keV).
             */
            inline void setMaxExcitationEnergy( double ex ) { maxExcitationEnergy=ex; }

            /**
             * Get the minimum possible excitation energy.
             */
            inline double getMinExcitationEnergy() const { return minExcitationEnergy; }

            /**
             * Get maximum possible excitation energy.
             */
            inline double getMaxExcitationEnergy() const { return maxExcitationEnergy; }

            /**
             * Get the nominal width of the state.
             */
            double getNominalWidth() const { return nominalWidth; }
        
            /**
            * Get position.
            */
            inline const TVector3& getPosition() const { return position; }

            /**
            * Set position.
            * @param pos XYZ Position vector.
            */
            inline void setPosition(TVector3 pos) { position=pos; }

            /**
            * Set four-momentum in lab frame.
            * (The kinetic energy, T, is determined from (mx+T)^2 = mx^2+p^2,
            * where mx is the rest mass plus the excitation energy and p is 
            * the momentum.)
            * @param p Momentum vector in keV/c.
            */
            void setFourMomentumLab(const TVector3& p);

            /**
            * Set four-momentum in lab frame.
            * (The momentum, p, is determined from (mx+T)^2 = mx^2+p^2,
            * where mx is the rest mass plus the excitation energy and T=ekin is 
            * the kinetic energy.)
            * @param T Kinetic energy in keV.
            * @param direction Direction of motion.
            */
            void setFourMomentumLab(double ekin, const TVector3& d);

            /**
            * Set four-momentum in lab frame.
            * (The excitation energy is determined from E^2 = mx^2+p^2,
            * where (E,p) is the four-momentum and mx is the rest mass plus
            * the excitation energy.)
            * @param p4 Four-momentum in (keV/c,keV).
            */
            void setFourMomentumLab(const TLorentzVector& p4);

            /**
            * Get four-momentum vector in lab frame.
            */
            inline const TLorentzVector& getFourMomentumLab() const { return fourMomentumLab; }

            /**
            * Get kinetic energy in lab frame.
            */
            double getKineticEnergyLab() const;

            /**
            * Get direction in lab frame.
            * Returns null-vector if particle is at rest in lab.
            */
            const TVector3& getDirectionLab() const;

            /**
            * Get direction in rest frame of parent (RFOP).
            * Returns null-vector if particle is at rest in parent frame.
            * Returns direction in lab if particle has no parent.
            */
            TVector3 getDirectionRFOP() const;

            /**
            * Get momentum vector in lab frame.
            */
            inline TVector3 getMomentumLab() const { return fourMomentumLab.Vect(); }

            /**
            * Get velocity vector in lab frame.
            */
            inline TVector3 getVelocityLab() const { return fourMomentumLab.BoostVector(); }

            /**
             * Get parent.
             * This will be nullptr if this Particle has no parent.
             */
            inline const Particle* getParent() const { return parent; }

            /**
             * Returns true if particle has width larger than 'minWidth'.
             * Returns false otherwise.
             */
            bool hasFiniteWidth() const;
            
            /**
             * Get parent.
             * This will be nullptr if this Particle has no parent.
             */
            inline Particle* getParent() { return parent; }


            using Propagator = std::shared_ptr<propagator::ParticlePropagator>;
            Propagator getPropagator() const {
                return propagator;
            }

            void setPropagator(Propagator propagator) {
                Particle::propagator = propagator;
            }

            /**
             * Convenience method - will delegate to underlying ParticlePropagator.
             * @copydoc simX::propagtor::ParticlePropagator::propagate()
             * @pre propagator != nullptr
             */
            void propagate(const Layer& layer, Particle& part, double range = -1, double * nonIonizing = nullptr);

        private:
            const AUSA::EnergyLoss::Ion ion;
            const bool tracking;
            TVector3 position, direction;
            TLorentzVector fourMomentumLab; 
            double excitationEnergy;
            double nominalExcitationEnergy, nominalWidth;
            Particle* const parent;
            const double minWidth;
            double minExcitationEnergy, maxExcitationEnergy;

            Propagator propagator;

            /**
            * Update lab-frame kinetic energy (the momentum, p, is scaled accordingly).
            * @param s Kinetic energy in keV.
            */
            void updateKineticEnergyLab(double ekin);

    };
}

#endif	/* PARTICLE_H */
