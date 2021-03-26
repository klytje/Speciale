
#include <unittest++/UnitTest++.h>

#include <simX/Particle.h>
#include <simX/Layer.h>
#include "simX/propagator/Vacuum.h"

#include <ausa/constants/Mass.h>
#include <ausa/util/memory>
#include <ausa/geometry/Box.h>
#include <ausa/eloss/Material.h>

using namespace simX;
using namespace simX::propagator;
using namespace AUSA::EnergyLoss;
using AUSA::Geometry::Box;

SUITE(VacuumTest) {

    class SetupFixture {
    public:
        SetupFixture()
                : p(Ion(0,0)),
                  b0(Material::predefined("Silicon"),std::make_unique<Box>(100.,100.,100.,TVector3(2.,-3.,0.),TVector3(0,0,-1),TVector3(0,1,0)), false),
                  b1(Material::predefined("Silicon"),std::make_unique<Box>(100.,100.,100.,TVector3(2.,-3.,200),TVector3(0,0,-1),TVector3(0,1,0)), false)
        {
            p.setPosition({0,0,-100});
            p.setFourMomentumLab(1E3, {0,0,1});

            layers.push_back(&b0);
            layers.push_back(&b1);

            nothing = (size_t) -1;
        }

        Particle p;
        std::vector<const Layer*> layers;
        Layer b0, b1;

        size_t nothing;
    };

    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST_FIXTURE(SetupFixture, NoNextLayerIfEmptyVector) {
        std::vector<const Layer*> layers;

        auto next = propagateInVacuum(layers, p);

        CHECK_EQUAL(nothing, next);
    }

    TEST_FIXTURE(SetupFixture, ParticleWillHitBox) {
        auto next = propagateInVacuum(layers, p);

        CHECK_EQUAL(0, next);
    }

    TEST_FIXTURE(SetupFixture, ParticleMovingAwayFromBoxWillNotHitBox) {
        p.setFourMomentumLab(1E3, {0,0,-1});

        auto next = propagateInVacuum(layers, p);

        CHECK_EQUAL(nothing, next);
    }

    TEST_FIXTURE(SetupFixture, ParticlePositionHaveBeenUpdatedAfterPropagation) {
        propagateInVacuum(layers, p);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(-50, pos.Z(), 1E-5);
    }

    TEST_FIXTURE(SetupFixture, ParticleWithXOffset_HasXOffsetAfterProp) {
        p.setPosition({25,0,-100});

        propagateInVacuum(layers, p);

        auto& pos = p.getPosition();
        CHECK_CLOSE(25, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(-50, pos.Z(), 1E-5);
    }

    TEST_FIXTURE(SetupFixture, ParticleWillPropagateToNearestLayer) {
        layers[0] = &b1;
        layers[1] = &b0;

        auto next = propagateInVacuum(layers, p);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(-50, pos.Z(), 1E-5);

        CHECK_EQUAL(1, next);
    }
}
