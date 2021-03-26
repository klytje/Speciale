
#ifndef SIMX_SCATTERINGINTEGRAL_H
#define SIMX_SCATTERINGINTEGRAL_H

#include "ParticlePropagator.h"

#include <memory>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include <ausa/eloss/EnergyLossCalculator.h>

#include "Math/Interpolator.h"
#include <TF1.h>


namespace simX {
    namespace propagator {
        class ScatteringIntegral {
            public:
                /**
                * Class for evaluating scattering cross-section integral, J(z), based on
                * W. MOLLER, G. POSPIECH and G. SCHRIEDER, Nucl. Instr. Meth. 130 (1975) 265
                * @param zmin Lower bound for tabulation of scattering integral (must be > 0)
                * @param zmax Upper bound for tabulation of scattering integral (must be > zmin)
                */
                ScatteringIntegral(double zmin = 0, double zmax = 0);
                ~ScatteringIntegral();

                /**
                 * Evaluate scattering integral
                 */
                double eval(double z);

                /**
                 * Evaluate inverse of scattering integral
                 * @param zlow Lower bound for determination of inverse
                 * @param zhigh Higher bound for determination of inverse
                 */
                double invert(double j, double zlow, double zhigh);

            private:
                /** 
                 * Load Lindhard's tabulation of f(z)
                 */
                void loadTabulation();

                /** 
                  * Integrand in Eq. (5) 
                  */
                Double_t finteq5(Double_t *x, Double_t *par);

                /** 
                  * J(z) as given in Eq. (5) 
                  */
                Double_t fjz(Double_t *z, Double_t *par);
                
                TF1 * Jz;

                double zmin, zmax;
                double delta;

                double zlh[15], flh[15];

                std::unique_ptr<ROOT::Math::Interpolator> fzInter; // Interpolation of f(z)
                std::unique_ptr<ROOT::Math::Interpolator> jzInter; // Interpolation of J(z)
        };
    }
}


#endif //SIMX_SCATTERINGINTEGRAL_H
