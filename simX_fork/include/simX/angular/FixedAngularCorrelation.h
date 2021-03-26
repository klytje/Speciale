//
// Created by munk on 05-08-15.
//

#ifndef SIMX_FIXEDANGULARCORRELATION_H
#define SIMX_FIXEDANGULARCORRELATION_H

#include "AngularCorrelation.h"

namespace simX {
    namespace angular {
        /**
         * AngularCorrelation that will always return the same angles.
         */
        class FixedAngularCorrelation : public AngularCorrelation {
        public:

            /**
             * Constructor that will convert from degrees to radians.
             */
            FixedAngularCorrelation(double theta, double phi);

            virtual void sampleAngles(double& theta, double& phi) const;

            double getPhi() const {
                return phi;
            }

            /**
             * Set phi directly.
             * No conversion is performed!
             */
            void setPhi(double phi) {
                FixedAngularCorrelation::phi = phi;
            }

            double getTheta() const {
                return theta;
            }

            /**
             * Set theta directly.
             * No conversion is performed!
             */
            void setTheta(double theta) {
                FixedAngularCorrelation::theta = theta;
            }

        private:
            double phi, theta;
        };
    }
}


#endif //SIMX_FIXEDANGULARCORRELATION_H
