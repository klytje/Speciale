//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>
#include "simX/parser/DAQTriggerParser.h"

using namespace simX::parser;
using namespace simX;
using std::string;

SUITE(DAQTriggerParserTest) {
        struct SetupFixture {
            SetupFixture() {
                names = {"A", "B", "C", "D"};
                set.resize(4);
                set.reset();
            }

            std::vector<string> names;
            boost::dynamic_bitset<size_t> set;
        };

        TEST_FIXTURE(SetupFixture, Sanity) {
            CHECK_EQUAL(1,1);
        }

    TEST_FIXTURE(SetupFixture, CheckSimpleExpressionIsTrue) {
        set[0] = true;
        auto f = parseTriggerFunction("A", names);
        auto b = f(set);

        CHECK(b);
    }

    TEST_FIXTURE(SetupFixture, CheckSimpleExpressionIsFalse) {
        set[0] = false;
        auto f = parseTriggerFunction("A", names);
        auto b = f(set);

        CHECK(!b);
    }

    TEST_FIXTURE(SetupFixture, CheckThatAndOperatorWorks) {
        auto f = parseTriggerFunction("A & B", names);

        CHECK_EQUAL(true, (set[0] = set[1] = true, f(set)));
    }

    TEST_FIXTURE(SetupFixture, CheckThatOrOperatorWorks) {
        auto f = parseTriggerFunction("A | B", names);
        CHECK_EQUAL(false, (set[0] = set[1] = false, f(set)));
        CHECK_EQUAL(true, (set[0] = true, set[1] = false, f(set)));
        CHECK_EQUAL(true, (set[0] = false, set[1] = true, f(set)));
    }

    TEST_FIXTURE(SetupFixture, CheckThatNotOperatorWorks) {
        auto f = parseTriggerFunction("!A", names);

        CHECK_EQUAL(true, (set[0] = false, f(set)));
    }


    TEST_FIXTURE(SetupFixture, DownscaledFactor2IsTrueFirstTimeNotSecondTime) {
        set[0] = true;
        auto f = parseTriggerFunction("A(2)", names);
        CHECK(!f(set));
        CHECK(f(set));
    }

    TEST_FIXTURE(SetupFixture, DownscaledFactor4IsTrueFourthTime) {
        set[0] = true;
        auto f = parseTriggerFunction("A(4)", names);
        CHECK(!f(set));
        CHECK(!f(set));
        CHECK(!f(set));
        CHECK(f(set));
    }

    TEST_FIXTURE(SetupFixture, DownscaledFactor4IsTrueEigthTime) {
        set[0] = true;
        auto f = parseTriggerFunction("A(4)", names);
                CHECK(!f(set));
                CHECK(!f(set));
                CHECK(!f(set));
                CHECK(f(set));

                CHECK(!f(set));
                CHECK(!f(set));
                CHECK(!f(set));
                CHECK(f(set));
    }

    TEST_FIXTURE(SetupFixture, NotDownscaledAndExpressionIsTrueEveryTime) {
        set[0] = true;
        set[1] = true;
        auto f = parseTriggerFunction("A & B", names);
        CHECK(f(set));
        CHECK(f(set));
        CHECK(f(set));
        CHECK(f(set));
        CHECK(f(set));
    }

    TEST_FIXTURE(SetupFixture, DownscaledAndExpressionTriggersWhenBothAreTrue) {
        set[0] = true;
        set[1] = true;
        auto f = parseTriggerFunction("A & B(2)", names);
                CHECK(!f(set));
                CHECK(f(set));
                CHECK(!f(set));
                CHECK(f(set));
                CHECK(!f(set));
                CHECK(f(set));
    }

    TEST_FIXTURE(SetupFixture, DownscaledOrExpressionFiresWhenFirstIsTrue) {
        set[1] = true;
        auto bAlone = parseTriggerFunction("B(3)", names);
        CHECK(!bAlone(set));
        CHECK(!bAlone(set));
        CHECK(bAlone(set));

        CHECK(!bAlone(set));
        CHECK(!bAlone(set));
        CHECK(bAlone(set));


        set[0] = true;
        set[1] = true;
        auto or_ = parseTriggerFunction("B(3) | A(2)", names);

        CHECK(!or_(set));
        CHECK(or_(set)); // A
        CHECK(or_(set)); // B

        CHECK(or_(set)); // A
        CHECK(!or_(set));
        CHECK(or_(set)); // B
    }

    TEST_FIXTURE(SetupFixture, AndExpressionTriggersEvenWhenDownscaledDoesnt) {
        set[0] = true;
        set[1] = true;
        auto f = parseTriggerFunction("A(3) | B(3) | A & B", names);
        CHECK(f(set));
        CHECK(f(set));
        CHECK(f(set));

        set[0] = true;
        set[1] = false;
        CHECK(!f(set));
        CHECK(!f(set));
        CHECK(f(set));
    }

    TEST_FIXTURE(SetupFixture, AndAndParses) {
        parseTriggerFunction("A && B", names);
    }
}
