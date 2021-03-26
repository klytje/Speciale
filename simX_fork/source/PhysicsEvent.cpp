//
// Created by munk on 10-08-15.
//

#include "simX/PhysicsEvent.h"

using namespace simX;

/*
	Iterators
*/
PhysicsEvent::iterator simX::begin(PhysicsEvent& s) {
	return begin(s.particles);
}

PhysicsEvent::iterator simX::end(PhysicsEvent& s) {
	return end(s.particles);
}

PhysicsEvent::const_iterator simX::begin(PhysicsEvent const& s) {
	return begin(s.particles);
}

PhysicsEvent::const_iterator simX::end(PhysicsEvent const& s) {
	return end(s.particles);
}