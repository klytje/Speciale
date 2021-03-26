#include <unittest++/UnitTest++.h>
#include <memory>

#include "simX/parser/grammar/IonGrammar.h"

using namespace simX::parser;
using AUSA::EnergyLoss::Ion;
using std::string;

SUITE(IonParserTest) {

    TEST(BeryliumHasA4) {
        ElementSymbolTable e;
        CHECK_EQUAL(4, e.at("Be"));
    }

    TEST(RoengtiumHasA111) {
        ElementSymbolTable e;
        CHECK_EQUAL(111, e.at("Rg"));
    }

    TEST(Be9HasA9) {
        Ion ion;
        AzeIsotopeGrammar grammar;

        const string s = "Be9";
        bool r = parse(s.begin(), s.end(), grammar, ion);

        CHECK(r);
        CHECK_EQUAL(9, ion.getA());
    }

    TEST(Be9HasZ4) {
        Ion ion;
        AzeIsotopeGrammar grammar;

        const string s = "Be9";
        bool r = parse(s.begin(), s.end(), grammar, ion);

        CHECK(r);
        CHECK_EQUAL(4, ion.getZ());
    }

    TEST(GammaIsZ0A0) {
        PredefinedIons e;
        auto i = e.at("g");
        CHECK_EQUAL(0, i.getZ());
        CHECK_EQUAL(0, i.getA());
    }

    TEST(ProtonIsZ1A1) {
        PredefinedIons e;
        auto i = e.at("p");
        CHECK_EQUAL(1, i.getZ());
        CHECK_EQUAL(1, i.getA());
    }

    TEST(DeutronIsZ1A2) {
        PredefinedIons e;
        auto i = e.at("d");
        CHECK_EQUAL(1, i.getZ());
        CHECK_EQUAL(2, i.getA());
    }

    TEST(TritonIsZ1A3) {
        PredefinedIons e;
        auto i = e.at("t");
        CHECK_EQUAL(1, i.getZ());
        CHECK_EQUAL(3, i.getA());
    }

    TEST(AlphaIsZ2A4) {
        PredefinedIons e;
        auto i = e.at("a");
        CHECK_EQUAL(2, i.getZ());
        CHECK_EQUAL(4, i.getA());
    }

    TEST(NeutronIsZ0A1) {
        PredefinedIons e;
        auto i = e.at("n");
        CHECK_EQUAL(0, i.getZ());
        CHECK_EQUAL(1, i.getA());
    }

    //
    // Combined parser
    //

    TEST(CombinedAlsoGivesBe9HasZ4) {
        Ion ion;
        IonGrammar grammar;

        const string s = "Be9";
        bool r = parse(s.begin(), s.end(), grammar, ion);

        CHECK(r);
        CHECK_EQUAL(4, ion.getZ());
    }

    TEST(CombinedAlsoGivesBe9HasA9) {
        Ion ion;
        IonGrammar grammar;

        const string s = "Be9";
        bool r = parse(s.begin(), s.end(), grammar, ion);

        CHECK(r);
        CHECK_EQUAL(9, ion.getA());
    }

    TEST(CombinedParserUnderstandsProton) {
        Ion ion;
        IonGrammar grammar;

        const string s = "p";
        bool r = parse(s.begin(), s.end(), grammar, ion);

        CHECK(r);
        CHECK_EQUAL(1, ion.getA());
        CHECK_EQUAL(1, ion.getZ());
    }
}
