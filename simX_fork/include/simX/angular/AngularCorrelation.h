
#ifndef ANGULARCORRELATION_H
#define	ANGULARCORRELATION_H


namespace simX {
    namespace angular {
        /**
         * Base class for angular distributions and correlations.
         */
        class AngularCorrelation {

        public:
            AngularCorrelation() = default;
            virtual ~AngularCorrelation() = default;

            virtual void sampleAngles( double& theta, double& phi ) const = 0;
        };
    }
}

#endif	/* ANGULARCORRELATION_H */

