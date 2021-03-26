
#ifndef TWOBODYKINEMATICS_H
#define	TWOBODYKINEMATICS_H

#include "simX/generator/FinalStateGenerator.h"

#include <memory>
#include <simX/Particle.h>

namespace simX {
    
    /**
    * Derived class for generating two-body final states according to
    * two-body kinematics.
    * @param m0 Rest mass of parent.
    * @param ex0 Excitation energy of parent.
    * @param m1 Rest mass of 1st decay product.
    * @param ex1 Excitation energy of 1st decay product.
    * @param m2 Rest mass of 2nd decay product.
    * @param ex2 Excitation energy of 2nd decay product.
    */
    class TwoBodyKinematics : public FinalStateGenerator {

        public:
            TwoBodyKinematics( const Particle &parent, const Particle &d1, const Particle &d2 );
            virtual ~TwoBodyKinematics() = default;

            virtual std::vector<TLorentzVector>& getFourMomenta() override;
            
        protected:

        private:
//            double massParent, massDaughter1, massDaughter2;
//            double excitationParent, excitationDaughter1, excitationDaughter2;
            const Particle &parent, &d1, &d2;

            std::vector<TLorentzVector> pDaughters;
    };
}

#endif	/* TWOBODYKINEMATICS_H */

