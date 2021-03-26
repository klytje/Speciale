
#ifndef FINALSTATEGENERATOR_H
#define	FINALSTATEGENERATOR_H

#include <TLorentzVector.h>
#include <vector>


namespace simX {
    
    /**
    * Base class for generating final state of nuclear process
    */
    class FinalStateGenerator {
        public:
            FinalStateGenerator() = default;
            virtual ~FinalStateGenerator() = default;

            /**
            * Returns four-momenta of daughters
            * @param arg Input arguments
            */
            virtual std::vector<TLorentzVector>& getFourMomenta() = 0;

        protected:

        private:
    };
}

#endif	/* FINALSTATEGENERATOR_H */

