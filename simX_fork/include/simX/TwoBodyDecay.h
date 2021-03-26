
#ifndef TWOBODYDECAY_H
#define	TWOBODYDECAY_H

#include "NBodyDecay.h"

namespace simX {
    
    /**
    * Derived class for simulating two-body breakups.
    * @param parent Reference to parent nucleus.
    * @param daughter1 1st decay product.
    * @param daughter2 2nd decay product.
    * @param ac Angular correlation function.
    * @param lorb Orbital angular momentum (used for calculating penetrability and shift function).
    * @param ampCalc Weight calculator.
    */
    class TwoBodyDecay : public NBodyDecay {

        public:
            TwoBodyDecay( Particle& parent, Particle daughter1, Particle daughter2, std::unique_ptr<angular::AngularCorrelation> ac=nullptr, int lorb=0, std::unique_ptr<WeightCalculator> wCalc=nullptr );
            virtual ~TwoBodyDecay() = default;
    };
}

#endif	/* TWOBODYDECAY_H */

