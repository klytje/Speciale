//
// Created by jesper on 4/1/16.
//

#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <simX/propagator/NonIonizingPropagator.h>
#include <simX/parser/BeamParser.h>
#include "simX/CompoundFormation.h"
#include <ausa/json/JSONUtil.h>
#include <ausa/eloss/ListOfMaterials.h>
#include <ausa/eloss/Material.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>

#include "simX/parser/ReactionParser.h"
#include "simX/parser/ConfigParser.h"

#include "simX/propagator/IonizingPropagator.h"
#include "simX/parser/ElossFactory.h"


SUITE(ElossFactoryTest) {
    struct Fixture {
        simX::parser::ReactionParser rParser;
    };

    TEST_FIXTURE(Fixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }


    TEST_FIXTURE(Fixture, TabulationTypeCanBeSRIM) {
        auto chain = rParser.parseString("beam: C12");
        rapidjson::Document ref;
        ref.Parse(R"({"tabulation": "SRIM13","target_propagator": "ELOSS"})");

        AUSA::EnergyLoss::Layer layer(AUSA::EnergyLoss::Material::predefined("Silicon"), nullptr, true);

        auto calc = simX::parser::elossCalculator(ref);
        std::unique_ptr<AUSA::EnergyLoss::EnergyLossCalculator> p = calc(layer, chain->getBeam());


        auto cast = dynamic_cast<AUSA::EnergyLoss::EnergyLossRangeInverter*>(p.get());
        CHECK(cast != nullptr);
    }

    TEST_FIXTURE(Fixture, TabulationTypeCanBeGEANT) {
        auto chain = rParser.parseString("beam: C12");
        rapidjson::Document ref;
        ref.Parse(R"({"tabulation": "GEANT","target_propagator": "ELOSS"})");

        AUSA::EnergyLoss::Layer layer(AUSA::EnergyLoss::Material::predefined("Silicon"), nullptr, true);

        auto calc = simX::parser::elossCalculator(ref);
        std::unique_ptr<AUSA::EnergyLoss::EnergyLossCalculator> p = calc(layer, chain->getBeam());


        auto cast = dynamic_cast<AUSA::EnergyLoss::EnergyLossRangeInverter*>(p.get());
        CHECK(cast == nullptr);
    }

    TEST_FIXTURE(Fixture, TabulationTypeIsUnsupported) {
        auto chain = rParser.parseString("beam: C12");
        rapidjson::Document ref;
        ref.Parse(R"({"tabulation": "MythicTabulation","target_propagator": "ELOSS"})");

        AUSA::EnergyLoss::Layer layer(AUSA::EnergyLoss::Material::predefined("Silicon"), nullptr, true);

        auto calc = simX::parser::elossCalculator(ref);
        CHECK_THROW(calc(layer, chain->getBeam()), std::invalid_argument);


    }
}