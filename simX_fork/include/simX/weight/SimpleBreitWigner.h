
#ifndef SIMPLEBREITWIGNER_H
#define	SIMPLEBREITWIGNER_H

#include "simX/weight/WeightCalculator.h"


namespace simX {

    // Forward class declaration
    class NBodyDecay;
    class Particle;
        
    /**
    * Derived class for calculating statistical weights of N-body decays
    * using a simple Breit-Wigner profile for the excitaiton-energy of
    * the parent nucleus.
    */
    class SimpleBreitWigner : public WeightCalculator {

        public:
            SimpleBreitWigner( NBodyDecay& proc );
            virtual ~SimpleBreitWigner();
    
            virtual double getWeight() const override;

        protected:

        private:
            bool isBroad;
            double EX0, G0;
            Particle * parent;
    };
}

#endif	/* SIMPLEBREITWIGNER_H */

