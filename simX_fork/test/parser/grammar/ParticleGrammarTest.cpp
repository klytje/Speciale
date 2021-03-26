#include <unittest++/UnitTest++.h>
#include <memory>

#include "simX/parser/grammar/ParticleGrammar.h"

#include "ausa/constants/Mass.h"
#include "ausa/constants/Constants.h"

using namespace simX::parser;
using AUSA::EnergyLoss::Ion;
using std::string;

SUITE(ParticleParserTest) {

    Particle parseString(const string& s) {
        using boost::spirit::qi::phrase_parse;
        using boost::spirit::ascii::space;

        Particle particle;
        bool r = phrase_parse(s.begin(), s.end(), ParticleGrammar{}, space, particle);
        if (!r) throw std::invalid_argument("failed to parse '" + s + "'");
        return particle;
    }

    TEST(IonIsParticleWithNoExcitation) {
        auto particle = parseString("Be9  ");

        CHECK_EQUAL(0., particle.ex());
        CHECK_EQUAL("Be9", particle.ion.getName());
        CHECK_CLOSE(AUSA::Constants::isotopeMass("Be9"), particle.ion.getMass(), 1E-2);
    }

    TEST(ParticleCanBeExcited) {
        auto particle = parseString("Be9 Ex: 200keV");

        CHECK_EQUAL(200., particle.ex());
        CHECK_EQUAL("Be9", particle.ion.getName());
        CHECK_CLOSE(AUSA::Constants::isotopeMass("Be9"), particle.ion.getMass(), 1E-2);
    }

    TEST(ParticleCanHaveAWidth) {
        auto particle = parseString("Be9 G0: 200keV");

        CHECK_EQUAL(200., particle.G0());
        CHECK_EQUAL("Be9", particle.ion.getName());
        CHECK_CLOSE(AUSA::Constants::isotopeMass("Be9"), particle.ion.getMass(), 1E-2);
    }


    TEST(ParticleCanHaveWidthAndExcitation) {
        auto particle = parseString("Be9 Ex: 900keV G0: 200keV");

        CHECK_EQUAL(900., particle.ex());
        CHECK_EQUAL(200., particle.G0());
    }

    TEST(ParticleCanHaveWidthInTime) {
        auto expect = AUSA::HBAR / 10E-9;
        auto particle = parseString("Be9 G0: 10ns");

        CHECK_EQUAL(expect, particle.G0());
    }

    TEST(OrderOfExcitationAndWidthCanBeSwitched) {
        auto particle = parseString("Be9 G0: 900keV Ex: 200keV");

        CHECK_EQUAL(900., particle.G0());
        CHECK_EQUAL(200., particle.ex());
    }

    TEST(WidthAndExcitationCanBeLowercase) {
        auto particle = parseString("Be9 g0: 900keV ex: 200keV");

        CHECK_EQUAL(900., particle.G0());
        CHECK_EQUAL(200., particle.ex());
    }

    TEST(ParticleWithoutNoTrackIsTracked) {
        auto particle = parseString("Be9 g0: 900keV ex: 200keV");

        CHECK(particle.doTracking());
    }

    TEST(ParticleWithNoTrackIs__NOT__Tracked) {
        auto particle = parseString("Be9 notrack");

        CHECK(!particle.doTracking());
    }

    TEST(NoTrackWithCompatibleWithOtherOptions) {
        auto particle = parseString("Be9 g0: 900keV notrack");

        CHECK(!particle.doTracking());

        particle = parseString("Be9 ex: 200keV notrack");

        CHECK(!particle.doTracking());

        particle = parseString("Be9 ex: 200keV g0: 900keV notrack");

        CHECK(!particle.doTracking());

        particle = parseString("Be9 g0: 900keV notrack ex: 200keV");

        CHECK(!particle.doTracking());
        CHECK_CLOSE(200, particle.ex(), 1E-5);
        CHECK_CLOSE(900, particle.G0(), 1E-5);

        particle = parseString("Be9 notrack ex: 200keV g0: 900keV");

        CHECK(!particle.doTracking());
        CHECK_CLOSE(200, particle.ex(), 1E-5);
        CHECK_CLOSE(900, particle.G0(), 1E-5);

        particle = parseString("Be9 g0: 900keV ex:200keV notrack ");

        CHECK(!particle.doTracking());
    }

    TEST(NoTrackIsNotCaseSensitive) {
        auto particle = parseString("Be9 g0: 900keV ex:200keV noTrack ");

        CHECK(!particle.doTracking());
    }
}
