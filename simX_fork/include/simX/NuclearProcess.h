
#ifndef DECAY_H
#define	DECAY_H

#include "simX/Particle.h"
#include "simX/weight/WeightCalculator.h"
#include "simX/generator/FinalStateGenerator.h"
#include <vector>
#include <memory>


namespace simX {
    
    /**
    * Base class for reactions and multi-particle breakups.
    */
    class NuclearProcess {
        public:
            NuclearProcess() = default;
            virtual ~NuclearProcess() = default;

            /**
            * Determines the momenta of the decay products
            * and boosts to the lab frame.
            */
            virtual void runProcess() = 0;

            /**
            * Returns pointers to reaction/decay products.
            */
            virtual std::vector<Particle*>& getDaughters() = 0;

            /**
             * Returns pointers to reaction/decay products.
             */
            virtual const std::vector<Particle*>& getDaughters() const = 0;
            
            /**
             * Return the statistical weight of the process
             */
            virtual double getWeight() const = 0;

            /**
             * Set weight calculator
             */
            virtual void setWeightCalculator(std::unique_ptr<WeightCalculator>) = 0;
            
            /**
             * Returns nucleon number (conserved in all nuclear processes)
             */
            virtual int getNumberOfNucleons() const = 0;
            
            /**
             * Returns (initial) proton number (conserved in non-weak nuclear processes)
             */
            virtual int getNumberOfProtonsInInitialState() const = 0;

    };
}

#endif	/* DECAY_H */

