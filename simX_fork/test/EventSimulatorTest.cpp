
#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <TVector3.h>
#include <iostream>
#include <simX/Beam.h>
#include <ausa/util/memory>
#include <ausa/constants/Mass.h>
#include <simX/Target.h>
#include <simX/EventSimulator.h>
#include <simX/propagator/NonIonizingPropagator.h>
#include <ausa/geometry/Box.h>

#include "simX/Particle.h"
#include "simX/parser/ReactionParser.h"

using namespace simX;
using namespace std;
using namespace simX::parser;
using namespace AUSA::EnergyLoss;

SUITE(EventSimulatorTest) {

    Target buildTarget(int n) {
        std::vector<Layer> v;
        auto silicon = Material::predefined("Silicon");
        for (int i = 1; i <= n; i++) {
//            auto box = make_unique<Layer>();
            v.emplace_back(silicon, make_unique<AUSA::Geometry::Box>(100.,100.,100.,TVector3(2.,-3,(i-1)*100.),TVector3(0,0,-1),TVector3(0,1,0)), i==n);
        }

        return Target{move(v), TVector3(2.,-3,100.)};
    }
    
    void setProp(ProcessChain& c) {
        auto prop = std::make_shared<propagator::NonIonizingPropagator>();
        c.getBeam().setPropagator(prop);

        for (auto& i : c.getTree()) {
            for (auto& j : i->getDaughters()) {
                j->setPropagator(prop);
            }
        }
    }

    class SetupFixture {
        public:
            SetupFixture()
                : beam(1E-4, 0, 0, -10, 0, 0), target(buildTarget(1))
            {
                prop = std::make_shared<propagator::NonIonizingPropagator>();
            }

            ReactionParser parser;
            Beam beam;
            Target target;
            std::shared_ptr<propagator::NonIonizingPropagator> prop;
    };

    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST_FIXTURE(SetupFixture, EventSimulatorUpdatesDirection) {
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                AD: FIXED(theta="90" phi="90")
                He4
                Be8
            }
            )");

        setProp(*c);

        auto i = begin(*c);

        EventSimulator sim(beam, target, *c);
        sim.run();

        ++i;
        // Access alpha particle
        NuclearProcess& proc = *i;
        auto& d = proc.getDaughters();
        auto d1 = d[0];
        auto& dir1 = d1 -> getDirectionLab();
        CHECK_CLOSE( 0., dir1.X(), 1E-3 );
        CHECK_CLOSE( 1., dir1.Y(), 1E-3 );
        CHECK_CLOSE( 0., dir1.Z(), 1E-3 );
    }

    TEST_FIXTURE(SetupFixture, BeamParticleIsPropagatedIntoTarget) {
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                AD: FIXED(theta="90" phi="90")
                He4
                Be8
            }
            )");

        auto& b = c->getBeam();
        setProp(*c);

        EventSimulator sim(beam, target, *c);
        sim.run();


        auto& pos = b.getPosition();
        CHECK_CLOSE( 0., pos.X(), 1E-3 );
        CHECK_CLOSE( 0, pos.Y(), 1E-3 );

        CHECK(pos.Z() < 50 && pos.Z() > -50);
    }

    TEST_FIXTURE(SetupFixture, BeamParticleIsPropagatedIntoActiveLayer) {
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            )");

        auto& b = c->getBeam();
        setProp(*c);

        target = buildTarget(2); // Two layer target

        EventSimulator sim(beam, target, *c);
        sim.run();


        auto& pos = b.getPosition();
        CHECK_CLOSE( 0., pos.X(), 1E-3 );
        CHECK_CLOSE( 0, pos.Y(), 1E-3 );

        CHECK(pos.Z() < 150 && pos.Z() > 50);
    }

    TEST_FIXTURE(SetupFixture, PhysicsEventContainsOnlyFinalStates) {
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                He4
                Be8
            }
            )");

        setProp(*c);

        EventSimulator sim(beam, target, *c);
        auto event = sim.run();

        CHECK_EQUAL(2, event.size());
    }

    TEST_FIXTURE(SetupFixture, TheEventContains3Alphas) {
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                He4
                Be8 -> {
He4
He4
}
            }
            )");

        setProp(*c);

        EventSimulator sim(beam, target, *c);
        auto event = sim.run();

        CHECK_EQUAL(3, event.size());

        for (auto p : event) {
            CHECK_EQUAL("He4", p->getName());
        }
    }

    TEST_FIXTURE(SetupFixture, FinalStatesArePropagatedToSurfaceOfSingleLayerTarget) {
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                AD: FIXED(theta="90" phi="90")
                He4
                Be8
            }
            )");

        setProp(*c);

        EventSimulator sim(beam, target, *c);
        sim.run();

        auto& e = c->findDecayOf("C12")->getDaughters();

        {
            auto he = *begin(e);

            CHECK_CLOSE(1, he->getDirectionLab().Y(), 1E-4);

            auto& pos = he->getPosition();
            CHECK_CLOSE( 0., pos.X(), 1E-3 );
            CHECK_CLOSE( -3+50, pos.Y(), 1E-3 );

            CHECK(pos.Z() < 50 && pos.Z() > -50);

        }

        {
            auto be = *(++begin(e));

            CHECK_CLOSE(-1, be->getDirectionLab().Y(), 1E-4);

            auto& pos = be->getPosition();
            CHECK_CLOSE( 0., pos.X(), 1E-3 );
            CHECK_CLOSE( -3-50, pos.Y(), 1E-3 );

            CHECK(pos.Z() < 50 && pos.Z() > -50);
        }
    }

    TEST_FIXTURE(SetupFixture, FinalStatesArePropagatedToSurfaceOfSingleLayerTargetWhenMovingInZ) {
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                AD: FIXED(theta="0" phi="0")
                He4
                Be8
            }
            )");

        setProp(*c);

        EventSimulator sim(beam, target, *c);
        sim.run();

        auto& e = c->findDecayOf("C12")->getDaughters();

        {
            auto p = *begin(e);

            CHECK_EQUAL("He4", p->getName());

            auto dir = p->getDirectionLab();
            CHECK_CLOSE( 0., dir.X(), 1E-3 );
            CHECK_CLOSE( 0, dir.Y(), 1E-3 );
            CHECK_CLOSE( 1, dir.Z(), 1E-3 );

            auto& pos = p->getPosition();
            CHECK_CLOSE( 0., pos.X(), 1E-3 );
            CHECK_CLOSE( 0, pos.Y(), 1E-3 );
            CHECK_CLOSE( 50, pos.Z(), 1E-3 );
        }

        {
            auto p = *(++begin(e));

            CHECK_EQUAL("Be8", p->getName());

            auto dir = p->getDirectionLab();
            CHECK_CLOSE( 0., dir.X(), 1E-3 );
            CHECK_CLOSE( 0, dir.Y(), 1E-3 );
            CHECK_CLOSE( -1, dir.Z(), 1E-3 );

            auto& pos = p->getPosition();
            CHECK_CLOSE( 0., pos.X(), 1E-3 );
            CHECK_CLOSE( 0, pos.Y(), 1E-3 );
            CHECK_CLOSE( -50, pos.Z(), 1E-3 );
        }
    }

    TEST_FIXTURE(SetupFixture, FinalStatesArePropagatedToSurfaceOf_Double_LayerTargetWhenMovingInZ) {
        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                AD: FIXED(theta="0" phi="0")
                He4
                Be8
            }
            )");

        setProp(*c);

        target = buildTarget(2);

        EventSimulator sim(beam, target, *c);
        sim.run();

        auto& e = c->findDecayOf("C12")->getDaughters();

        {
            auto p = *begin(e);

            auto& pos = p->getPosition();
            CHECK_CLOSE( 0., pos.X(), 1E-3 );
            CHECK_CLOSE( 0, pos.Y(), 1E-3 );
            CHECK_CLOSE( 150, pos.Z(), 1E-3 );
        }

        {
            auto p = *(++begin(e));

            auto& pos = p->getPosition();
            CHECK_CLOSE( 0., pos.X(), 1E-3 );
            CHECK_CLOSE( 0, pos.Y(), 1E-3 );
            CHECK_CLOSE( -50, pos.Z(), 1E-3 );
        }
    }


    TEST_FIXTURE(SetupFixture, IfNominalBeamCannotHitAllLayersThrowException) {
        Beam beam(1000, 10000, 10000, -100, 0, 0, nullptr, nullptr, nullptr, nullptr);

        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                AD: FIXED(theta="90" phi="90")
                He3
                Be9
            }
            )");

        CHECK_THROW(
                EventSimulator(beam, target, *c),
                std::runtime_error
        );
    }

    TEST_FIXTURE(SetupFixture, TargetMiss) {
        auto f = std::make_unique<TF1>("", "1", 150, 151);

        Beam beam(1000, 2, -3, -100, 0, 0, nullptr, move(f), nullptr, nullptr);

        auto c = parser.parseString(R"(
            beam: He3
            target: Be9
            -> {
                AD: FIXED(theta="90" phi="90")
                He3
                Be9
            }
            )");

        setProp(*c);
        EventSimulator simulator(beam, target, *c);

        auto e = simulator.run();

        CHECK_EQUAL(0, e.size());
    }
}
