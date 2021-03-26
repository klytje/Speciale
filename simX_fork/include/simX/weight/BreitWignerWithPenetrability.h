
#ifndef BREITWIGNERWITHPENETRABILITY_H
#define	BREITWIGNERWITHPENETRABILITY_H

#include <memory>
#include "simX/weight/WeightCalculator.h"

class CoulombFunctions;

namespace simX {

    // Forward class declaration
    class NBodyDecay;
    class Particle;
    
    /**
    * Derived class for calculating statistical weights of two-body decays
    * based on (simplified) R-matrix formalism.
    */
    class BreitWignerWithPenetrability : public WeightCalculator {

        public:
            BreitWignerWithPenetrability( NBodyDecay& proc );
            virtual ~BreitWignerWithPenetrability();
    
            virtual double getWeight() const override;

        protected:

        private:
            bool isBroad;
            double EX0, G0, Q;
            Particle *parent, *daughter1, *daughter2;            
            double buffer; // to avoid numerical problems close to singularities
            std::unique_ptr<CoulombFunctions> CF;
    };
}

#endif	/* BREITWIGNERWITHPENETRABILITY_H */

