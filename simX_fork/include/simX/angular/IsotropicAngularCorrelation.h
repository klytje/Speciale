//
// Created by munk on 29-06-15.
//

#ifndef SIMX_ISOTROPICANGULARCORRELATION_H
#define SIMX_ISOTROPICANGULARCORRELATION_H

#include "AngularCorrelation.h"

#include <utility>

namespace simX {
    namespace angular {
        class IsotropicAngularCorrelation : public AngularCorrelation {
        public:
            typedef double degree;
            typedef double radian;
            typedef std::pair<degree, degree> Limits;

            IsotropicAngularCorrelation();
            IsotropicAngularCorrelation(const Limits& phi, const Limits& theta);
            IsotropicAngularCorrelation(degree phiMin, degree phiMax, degree thetaMin, degree thetaMax);

            virtual void sampleAngles(double& theta, double& phi) const;

        private:
            radian phiLow, phiDiff;
            radian thetaLow, thetaDiff;
        };
    }
}
#endif //SIMX_ISOTROPICANGULARCORRELATION_H
