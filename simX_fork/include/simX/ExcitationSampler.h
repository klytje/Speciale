
#ifndef EXCITATIONSAMPLER_h
#define EXCITATIONSAMPLER_h 1

// simX header files
#include "simX/ProcessChain.h"
#include "simX/Particle.h"

// C++ and ROOT header files
#include <vector>
#include <memory>

// AUSAlib header files
#include "ausa/util/memory"

// ROOT header files
#include "Math/DistSampler.h"
#include "Math/DistSamplerOptions.h"
#include "Math/Factory.h"
#include <TF1.h>


namespace simX {

    /**
    * Class for generating excitation energies (Ex).
    * @param chain Chain of nuclear processes to be simulated.
    */
    class ExcitationSampler {

        public:
            ExcitationSampler( ProcessChain& chain );
            ~ExcitationSampler();

            /**
            * Sample all variable excitation energies
            */
            void sampleExcitationEnergies();
            
            /**
            * Return number of variable excitation energies
            */            
            inline int getDimension() { return DIM; }

        private:
            int DIM;
            ProcessChain& chain;
            std::vector<Particle*> fatDaughters, slimDaughters;
            ROOT::Math::DistSampler *sampler;
            bool verbose, isInitialized;
            
            /**
            * Determine the sampling limits for each variable
            */            
            void determineLimits();
            
            /**
            * Product of weight functions (typically, breit wigners)
            */            
            double bwprod(double *x, double *p);

            /**
             * Initializes the excitation sampler
             */
            void init();

    };
}

#endif	/* EXCITATIONSAMPLER_H */
