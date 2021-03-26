
#ifndef ANGULARCORRELATIONTF1_H
#define	ANGULARCORRELATIONTF1_H

#include "simX/angular/AngularCorrelation.h"

#include <TF1.h>
#include <string>


namespace simX {
    namespace angular {
        /**
         * Angular distributions and correlations defined using ROOT's TF1 class.
         */
        class AngularCorrelationTF1 : public AngularCorrelation {

        public:
            AngularCorrelationTF1(std::string theta, std::string phi, double phiMin, double phiMax, double thetaMin, double thetaMax);

            virtual void sampleAngles( double& theta, double& phi ) const override;

        private:
            mutable TF1 fPhi, fTheta;
        };
    }
}

#endif	/* ANGULARCORRELATION_H */

