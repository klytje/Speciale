//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <simX/propagator/NonIonizingPropagator.h>
#include "simX/CompoundFormation.h"

#include "simX/parser/BeamParser.h"

#include "simX/propagator/IonizingPropagator.h"

using namespace simX::parser;
using namespace simX;
using std::string;

SUITE(BeamParserTest) {
    struct Fixture {
        BeamParser cParser;
    };

    TEST_FIXTURE(Fixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(Fixture, EnergyDistributionCanBeGaussian) {
        cParser.parse(R"(energy : 2000keV center: (0mm, 0mm, 0mm) dist_energy : GAUSSIAN(sigma = "10keV"))");
    }

    TEST_FIXTURE(Fixture, EnergyPerNucleon) {
        auto p = cParser.parse(R"(energy : 2000keV/A center: (0mm, 0mm, 0mm) dist_energy : GAUSSIAN(sigma = "10keV"))");
        CHECK(p.isPerNucleon());
    }

    TEST_FIXTURE(Fixture, EnergyNOTPerNucleon) {
        auto p = cParser.parse(R"(energy : 2000keV center: (0mm, 0mm, 0mm) dist_energy : GAUSSIAN(sigma = "10keV"))");
        CHECK(!p.isPerNucleon());
    }

    TEST_FIXTURE(Fixture, AngleDistributionCanBeGaussian) {
        cParser.parse(R"(energy : 2000keV center: (0mm, 0mm, 0mm) dist_theta: GAUSSIAN(sigma = "10deg"))");
    }

    TEST_FIXTURE(Fixture, AngleDistributionCanBeUniform) {
        cParser.parse(R"(energy : 2000keV center: (0mm, 0mm, 0mm) dist_theta: UNIFORM(low = "10deg" high="30deg"))");
    }

    TEST_FIXTURE(Fixture, AngleDistributionCanBePointDistribution) {
        cParser.parse(R"(energy : 2000keV center: (0mm, 0mm, 0mm) dist_theta: POINT(low = "10deg" high="30deg"))");
    }

    TEST_FIXTURE(Fixture, AngleDistributionCanBeTF1) {
        cParser.parse(R"(energy : 2000keV center: (0mm, 0mm, 0mm) dist_theta: POINT(low = "10deg" high="30deg" expression="1+x"))");
    }

    TEST_FIXTURE(Fixture, XYDistributionCanBeTF1) {
        cParser.parse(R"(energy : 2000keV center: (0mm, 0mm, 0mm) dist_xy: UNIFORM(radius = "10mm"))");
    }
}
