
#ifndef BREITWIGNERWITHGAMMAPHASESPACE_H
#define	BREITWIGNERWITHGAMMAPHASESPACE_H

#include <memory>
#include "simX/weight/WeightCalculator.h"


namespace simX {
    
    // Forward class declaration
    class NBodyDecay;
    class Particle;

    /**
    * Derived class for calculating statistical weights of gamma decays
    */
    class GammaPhaseSpace : public WeightCalculator {

        public:
            GammaPhaseSpace(NBodyDecay& proc );
            virtual ~GammaPhaseSpace();
    
            virtual double getWeight() const override;

        protected:

        private:
            int lorb;
            bool isBroad;
            double EX0, G0, Q;
            Particle *parent, *daughter1, *daughter2;            
    };
}

#endif	/* BREITWIGNERWITHGAMMAPHASESPACE_H */

