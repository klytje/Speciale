
#ifndef EVENTSIMULATOR_h
#define EVENTSIMULATOR_h 1

// simX header files
#include "simX/Beam.h"
#include "simX/Target.h"
#include "ProcessChain.h"
#include "simX/Layer.h"
#include "simX/PhysicsEvent.h"
#include "simX/ExcitationSampler.h"
#include "Logger.h"

// C++ and ROOT header files
#include <vector>
#include <memory>


namespace simX {

    /**
    * Class for simulating physics events.
    * @param beam Beam.
    * @param target Target.
    */
    class EventSimulator {

        public:
            EventSimulator( Beam& beam, Target& target, ProcessChain& chain );
            ~EventSimulator();

            PhysicsEvent run();

        private:
            Beam& beam;
            Particle& beamPart;
            std::vector<const Layer*> layers;
            std::vector<Particle*> finalStates, empty;
            size_t indexActiveLayer;
            ProcessChain& chain;
            std::unique_ptr<ExcitationSampler> exSampler;
            Long64_t event;

            /**
            * Initializes the beam particle, i.e., samples the energy, position, 
            * and direction from the user-specified distributions.
            * @param mode If mode=0, only nominal values are returned.
            */
            void sampleBeam( int mode);

            /**
             * Propagates the beam into the active layer of the target (and
             * through the inactive layers it encounters before reaching the
             * the active layer).
             * @param f Depth in active layer (fraction of full thickness).
             * @pre f is in the range [0, 1]
             */
            bool propagateBeam( double f );

            void orderTargetLayers(const Target& target);
            void propagateParticleOut(Particle& p);

            const Layer* findNextLayer(Particle& p, const Layer* last);




//            void runDecayChain();


    };
}

#endif	/* EVENTSIMULATOR_H */
