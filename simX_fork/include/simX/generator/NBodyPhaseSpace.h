
#ifndef NBODYPHASESPACE_H
#define	NBODYPHASESPACE_H

#include "simX/generator/FinalStateGenerator.h"

// ROOT libraries
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <memory>


namespace simX {
    class Particle;

    /**
    * Derived class for generating N-body final state according to
    * phase-space distribution.
    * @param etot Total available energy (in keV) in CM frame.
    * @param mult Number of daughters.
    * @param masses Rest masses of daughters (including excitation).
    */
    class NBodyPhaseSpace : public FinalStateGenerator {

        public:
            NBodyPhaseSpace(const Particle& p, const std::vector<Particle*>& d);

            virtual ~NBodyPhaseSpace() = default;
    
            virtual std::vector<TLorentzVector>& getFourMomenta() override;


            /**
             * How many samples should be use to determine maximum W.
             */
            static const size_t N_SAMPLES = 1000000;
        protected:

        private:
            bool isInitialized;
            int mult;
            const Particle& parent;
            const std::vector<Particle*>& daughtersPtr;
            TGenPhaseSpace phaseSpaceGenerator;
            double Wmax;
            std::vector<TLorentzVector> pDaughters;

            void init();
    };
}

#endif	/* NBODYPHASESPACE_H */

