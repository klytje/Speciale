
#include <unittest++/UnitTest++.h>

#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include <simX/propagator/NonIonizingPropagator.h>
#include <ausa/util/memory>
#include <ausa/eloss/Material.h>
#include <ausa/geometry/Box.h>

#include "simX/Particle.h"
#include "simX/Layer.h"

using namespace simX;
using namespace simX::propagator;
using namespace AUSA::EnergyLoss;

using AUSA::Geometry::Box;

SUITE(NonIonizingTest) {

    class SetupFixture {
    public:
        SetupFixture()
                : silicon("Silicon", 2.3290, 14, AUSA::Constants::atomicMass(14)),
                  p(Ion(0,0)),
                  b0(Material::predefined("Silicon"),std::make_unique<Box>(100.,100.,100.,TVector3(0.,0,0.),TVector3(0,0,-1),TVector3(0,1,0)), false)
        {
            p.setPosition({0,0,-50});
            p.setFourMomentumLab(1E3, {0,0,1});
        }

        Material silicon;
        Particle p;
        Layer b0;
        NonIonizingPropagator prop;
    };

    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST_FIXTURE(SetupFixture, HeadOnParticleWillTransverseLayer) {
        prop.propagate(b0, p, 100);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(50, pos.Z(), 1E-5);
    }

    TEST_FIXTURE(SetupFixture, IfNoRangeSpecifiedItWillMoveEntireLength) {
        prop.propagate(b0, p);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(50, pos.Z(), 1E-5);
    }

    TEST_FIXTURE(SetupFixture, PropagatorRespectsSubMaximalRange) {
        prop.propagate(b0, p, 50);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(0, pos.Z(), 1E-5);
    }


    TEST_FIXTURE(SetupFixture, ParticleInsideWillPropagateToSurface) {
        p.setPosition({0,0,0});

        prop.propagate(b0, p);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(50, pos.Z(), 1E-5);
    }

    TEST_FIXTURE(SetupFixture, ParticleInsideWillPropagateToSurfaceEvenIfNotOrientedAlongAxis) {
        p.setPosition({0,0,0});
        p.setFourMomentumLab(1, {1,1,1});

        prop.propagate(b0, p);

        auto& pos = p.getPosition();
        CHECK_CLOSE(50, pos.X(), 1E-5);
        CHECK_CLOSE(50, pos.Y(), 1E-5);
        CHECK_CLOSE(50, pos.Z(), 1E-5);
    }
}
