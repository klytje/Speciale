#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>

#include "simX/parser/grammar/BeamGrammar.h"

using namespace simX::parser;
using std::string;
using namespace std;

SUITE(BeamParserTest) {

    BeamPrototype parseString(const string& s) {
        using boost::spirit::qi::phrase_parse;
        using boost::spirit::ascii::space;

        BeamPrototype particle;
        bool r = phrase_parse(s.begin(), s.end(), BeamGrammar{}, space, particle);
        if (!r) throw std::invalid_argument("failed to parse '" + s + "'");
        return particle;
    }


    TEST(BeamHaveEnergy_2) {
        auto r = parseString(" energy: 2keV center: (2mm, 4mm, 10mm)");

        CHECK_EQUAL(2, r.energy.get().first);
    }

    TEST(BeamHaveCenter_3_4) {
        auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 10mm))");

        CHECK_EQUAL(3, get<0>(r.center.get()));
        CHECK_EQUAL(4, get<1>(r.center.get()));
    }

    TEST(AngleIsOptional) {
        auto r = parseString(" energy: 2keV center: (2mm, 4mm, 5mm)");

        CHECK(!r.theta.is_initialized());
    }

    TEST(BeamHaveAngle_3) {
        auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 5mm)
                                     theta: 3rad
)");

        CHECK_EQUAL(3, r.theta.get());
    }

    TEST(BeamHaveUniformAngleDist) {
      auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 5mm)
                                     theta: 3rad
                                     dist_theta: UNIFORM
)");

      CHECK_EQUAL("UNIFORM", r.thetaDist->first);

      CHECK(!r.thetaDist->second.is_initialized());
    }

    TEST(BeamHaveUniformXYDist) {
        auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 5mm)
                                     theta: 3rad
                                     dist_xy: UNIFORM
)");

        CHECK(r.xyDist.is_initialized());
        CHECK_EQUAL("UNIFORM", r.xyDist->first);

        CHECK(!r.xyDist->second.is_initialized());
    }


    TEST(OptionalArgumentsCanBeProvidedInParathesis) {
      auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 5mm)
                                     theta: 3rad
                                     dist_theta: UNIFORM(width = "10rad")
)");

      CHECK(r.thetaDist->second.is_initialized());

      auto& o = r.thetaDist->second.get();

      CHECK_EQUAL("10rad", o["width"]);
    }

    TEST(BeamHaveUniformEnergyDist) {
      auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 5mm)
                                     theta: 3rad
                                     dist_energy : UNIFORM
)");

      CHECK_EQUAL("UNIFORM", r.energyDist->first);

      CHECK(!r.energyDist->second.is_initialized());
    }

    TEST(OptionalArgumentsCanBeProvidedInParathesisToEnergy) {
      auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 5mm)
                                     theta: 3rad
                                     dist_energy : UNIFORM(width = "10keV")
)");

      CHECK(r.energyDist->second.is_initialized());

      auto& o = r.energyDist->second.get();

      CHECK_EQUAL("10keV", o["width"]);
    }

    TEST(BeamHavePhi_3) {
        auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 5mm)
                                     phi: 3rad
)");

        CHECK_EQUAL(3, r.phi.get());
    }


    TEST(BeamHavePhiDistributionPoint) {
        auto r = parseString(R"( energy: 2keV
                                     center: (3mm, 4mm, 5mm)
                                     dist_phi: POINT
)");

        CHECK(r.phiDist.is_initialized());

        auto name = r.phiDist->first;
        CHECK_EQUAL("POINT", name);
    }
}
