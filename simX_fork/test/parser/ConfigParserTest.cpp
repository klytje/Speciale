//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <simX/propagator/NonIonizingPropagator.h>
#include "simX/CompoundFormation.h"

#include "simX/parser/ReactionParser.h"
#include "simX/parser/ConfigParser.h"

#include "simX/propagator/IonizingPropagator.h"
#include "simX/propagator/GaussianStragglingPropagator.h"
#include "simX/propagator/MCStragglingPropagator.h"

#include <TF1.h>

using namespace simX::parser;
using namespace simX;
using std::string;

SUITE(ConfigParserTest) {
    struct Fixture {
        ConfigParser cParser;
        ReactionParser rParser;
    };

    TEST_FIXTURE(Fixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(Fixture, ELOSS_type_gives_IonizingPropagator) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"({"target_propagator": {"type" : "ELOSS"}})");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto ionizing = std::dynamic_pointer_cast<propagator::IonizingPropagator>(p);
        CHECK(ionizing != nullptr);
    }

    TEST_FIXTURE(Fixture, ELOSS_string_gives_IonizingPropagator) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"({"target_propagator": "ELOSS"})");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto ionizing = std::dynamic_pointer_cast<propagator::IonizingPropagator>(p);
        CHECK(ionizing != nullptr);
    }

    TEST_FIXTURE(Fixture, NOLOSS_tag_gives_NonIonizingPropagator) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"( {
            "target_propagator": {"type" : "NOLOSS"}} )");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto cast = std::dynamic_pointer_cast<propagator::NonIonizingPropagator>(p);
        CHECK(cast != nullptr);
    }

    TEST_FIXTURE(Fixture, TargetPropagatorDefaultsToEloss) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"( {
            "detection_propagator": {"type" : "ELOSS"}})");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto ionizing = std::dynamic_pointer_cast<propagator::IonizingPropagator>(p);
        CHECK(ionizing != nullptr);
    }

    TEST_FIXTURE(Fixture, DetectionPropagtionIsAlsoEnabled) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"( {
            "detection_propagator": {"type" : "NOLOSS"}})");

        auto b = chain->getBeam();

        config.applyDetectionPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto cast = std::dynamic_pointer_cast<propagator::NonIonizingPropagator>(p);
        CHECK(cast != nullptr);
    }

    TEST_FIXTURE(Fixture, GammaIsAlwaysNonIonizing) {
        auto chain = rParser.parseString("beam: g");
        auto config = cParser.parse(R"( {
            "detection_propagator": {"type" : "NOLOSS"}})");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto cast = std::dynamic_pointer_cast<propagator::NonIonizingPropagator>(p);
        CHECK(cast != nullptr);
    }

    TEST_FIXTURE(Fixture, ConfigCanContainTargetFile) {
        auto res = cParser.parse(R"({"target": "target.json"})");

        CHECK(res.hasTargetFile());
    }

    TEST_FIXTURE(Fixture, BeamPropagatorCanBeSetIndependently) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"({"beam_propagator": {"type" : "NOLOSS"}})");

        auto b = chain->getBeam();

        config.applyBeamPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto ionizing = std::dynamic_pointer_cast<propagator::NonIonizingPropagator>(p);
        CHECK(ionizing != nullptr);
    }

    TEST_FIXTURE(Fixture, GAUSSSTRAG_string_gives_GaussianStragglingPropagator) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"( {
            "target_propagator": "GAUSSSTRAG"})");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto gaussstrag = std::dynamic_pointer_cast<propagator::GaussianStragglingPropagator>(p);
        CHECK(gaussstrag != nullptr);
    }

    TEST_FIXTURE(Fixture, MCSTRAG_string_gives_MCStragglingPropagator) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"( {
            "target_propagator": "MCSTRAG"})");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto mcstrag = std::dynamic_pointer_cast<propagator::MCStragglingPropagator>(p);
        CHECK(mcstrag != nullptr);
    }

    TEST_FIXTURE(Fixture, Can_set_MCSTRAG_with_default_multiplier) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"({
            "target_propagator": {"type": "MCSTRAG"}})");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto mcstrag = std::dynamic_pointer_cast<propagator::MCStragglingPropagator>(p);
        CHECK(mcstrag != nullptr);
        CHECK_EQUAL(0.04, mcstrag->getCrossSectionMultiplier());
    }

    TEST_FIXTURE(Fixture, Can_set_MCSTRAG_with_multiplier) {
        auto chain = rParser.parseString("beam: C12");
        auto config = cParser.parse(R"({
            "target_propagator": {"type": "MCSTRAG", "multiplier":0.77}})");

        auto b = chain->getBeam();

        config.applyTargetPropagator(b);

        auto p = b.getPropagator();
        CHECK(p != nullptr);

        auto mcstrag = std::dynamic_pointer_cast<propagator::MCStragglingPropagator>(p);
        CHECK(mcstrag != nullptr);
        CHECK_CLOSE(0.77, mcstrag->getCrossSectionMultiplier(), 0.001);
    }
}
