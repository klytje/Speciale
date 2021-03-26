
#ifndef PHYSEVENT_h
#define PHYSEVENT_h 1

#include "simX/Particle.h"

namespace simX {

    class PhysicsEvent {

        public:
            PhysicsEvent(const std::vector<Particle*> &particles) : particles(particles) {}
            ~PhysicsEvent() {}

            size_t size() const {return particles.size(); }

            // Iterators for the isotopes.
            typedef std::vector<Particle*>::const_iterator iterator;
            typedef std::vector<Particle*>::const_iterator const_iterator;
            friend iterator begin(PhysicsEvent& s);
            friend iterator end(PhysicsEvent& s);
            friend const_iterator begin(const PhysicsEvent& s);
            friend const_iterator end(const PhysicsEvent& s);

        private:
            const std::vector<Particle*>& particles;
    };
}

#endif	/* PHYSEVENT_H */
