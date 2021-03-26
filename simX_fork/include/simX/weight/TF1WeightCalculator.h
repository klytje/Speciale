
#ifndef CUSTOMWEIGHTCALCULATOR_H
#define	CUSTOMWEIGHTCALCULATOR_H

#include <TF1.h>
#include "simX/weight/WeightCalculator.h"


namespace simX {

    // Forward class declaration
    class NBodyDecay;
    class Particle;
        
    /**
    * Derived class for calculating statistical weights using a custom TF1 function.
    */
    class TF1WeightCalculator : public WeightCalculator {

        public:
        TF1WeightCalculator( NBodyDecay& proc , std::string f);
            virtual ~TF1WeightCalculator();
    
            virtual double getWeight() const override;

        protected:

        private:
            Particle * parent;
            TF1 func;
    };
}

#endif	/* CUSTOMWEIGHTCALCULATOR_H */

