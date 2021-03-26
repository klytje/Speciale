
#include <unittest++/UnitTest++.h>
#include <memory>
#include <iostream>
#include <math.h>

#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include "ausa/util/memory"

#include "simX/ProcessChain.h"
#include "simX/ExcitationSampler.h"
#include "simX/parser/ReactionParser.h"

using namespace std;
using namespace simX;
using namespace AUSA::EnergyLoss;
using namespace simX::parser;

SUITE(ExcitationSamplerTest) {

    class SetupFixture {
        public:
            SetupFixture() {
            
                // Construct Process chain
                chain = px.parseString(R"(
                beam: He3
                target: B11
                -> {
                    AD: ISO
                    Li6 Ex: 2500keV G0: 250keV
                    -> {
                        AD: ISO
                        d
                        He4                
                    }
                    Be8 Ex: 3000keV G0: 1500keV
                    -> {
                        AD: ISO
                        He4
                        He4                
                    }
                }
                )");

            }

            ReactionParser px;
            unique_ptr<ProcessChain> chain;    
    };

    TEST_FIXTURE(SetupFixture, ExcitationSamplerAcceptsSimpleInput) {
        ExcitationSampler exs( *chain );
    }

    TEST_FIXTURE(SetupFixture, ReturnsCorrectDimension) {
        ExcitationSampler exs( *chain );
        int DIM=2;
        CHECK_EQUAL( DIM, exs.getDimension() );
    }
    
    TEST_FIXTURE(SetupFixture, CheckLimitsAreCorrectlySet) {
        auto first = chain -> getFirstProcess();
        // simulate compound formation
        first -> runProcess();
        // set nominal excitation energy of compound
        double ex = first -> getDaughters()[0] -> getExcitationEnergy();
        first -> getDaughters()[0] -> setNominalExcitationEnergy( ex );
        // create excitation sampler
        ExcitationSampler exs( *chain );
        for (auto& p : *chain) { // loop over processes
            auto& d = p.getDaughters();
            for (auto& daughter : d) { // loop over daughters
                double exmin = daughter -> getMinExcitationEnergy();
                double exmax = daughter -> getMaxExcitationEnergy();
                if (daughter->getName()=="Li6") {
                    CHECK_CLOSE(1473.76,exmin,1e-1);
                    CHECK_CLOSE(4662.39,exmax,1e-1);
                }
            } 
        }
                
    }

    TEST_FIXTURE(SetupFixture, CatchNonEnergyConservingDecay) {
        ReactionParser p;
        auto c = p.parseString("beam: C12 -> {p B11}");
        CHECK_THROW( ExcitationSampler exs(*c), std::invalid_argument );
    }

    TEST_FIXTURE(SetupFixture, CatchNonEnergyConservingBinaryReaction) {
        ReactionParser p;
        auto c = p.parseString("beam: p target: C12 -> {d C11}");
        CHECK_THROW( ExcitationSampler exs(*c), std::invalid_argument );
    }

    TEST_FIXTURE(SetupFixture, FinalStateParticlesMustHaveZeroWidth) {
        ReactionParser p;
        auto c = p.parseString("beam: He3 target: Be9 -> {a Be8 Ex: 0keV G0:100keV}");
        CHECK_THROW( ExcitationSampler exs(*c), std::invalid_argument );
    }

    TEST_FIXTURE(SetupFixture, CatchNonEnergyConservingThreeAlphaDecay) {
        ReactionParser p;
        auto c = p.parseString("beam: C12 -> {a a a}");
        CHECK_THROW( ExcitationSampler exs(*c), std::invalid_argument );
    }

    TEST_FIXTURE(SetupFixture, EnergyConservationEnsuredByAddingExcitationEnergy) {
        ReactionParser p;
        auto c = p.parseString("beam: C12 Ex:10000keV -> {a a a}");
        ExcitationSampler exs(*c);
    }

    TEST_FIXTURE(SetupFixture, CatchNonEnergyConservingTwoStepReaction) {
        ReactionParser p;
        auto c = p.parseString("beam: He3 target: Be9 -> {a Be8 Ex: -1000keV -> {a a} }");
        CHECK_THROW( ExcitationSampler exs(*c), std::invalid_argument );
    }

    TEST_FIXTURE(SetupFixture, FiniteWidthMakesItPossibleToSatisfyEnergyConservation) {
        ReactionParser p;
        auto c = p.parseString("beam: He3 target: Be9 -> {a Be8 Ex: -1000keV G0: 100keV -> {a a} }");
        ExcitationSampler exs(*c);
    }

    TEST_FIXTURE(SetupFixture, ThreeAlphaDecayEnergeticallyAllowedIfC12HasSufficientlyLargeWidth) {
        ReactionParser p;
        auto c = p.parseString("beam: C12 G0: 1500keV -> {a a a}");
        // default sampling range is +-5*width, so by setting the width
        // to 1.5 MeV we sample up to 7.5 MeV which is above the 3a threshold.
        ExcitationSampler exs(*c);
    }

/*
    // this test is rather slow (31 sec)
    TEST_FIXTURE(SetupFixture, CheckSamplerIsInitializedAndSamplingWorks) {
        ExcitationSampler exs( *chain );
        // sample
        exs.sampleExcitationEnergies();
        // check that generated energies are within prescribed limits
        double ex1=0;
        for (auto& p : *chain) { 
            auto& d = p.getDaughters();
            for (auto& daughter : d) { // loop over daughters
                ex1 = daughter -> getExcitationEnergy();
                double exmin = daughter -> getMinExcitationEnergy();
                double exmax = daughter -> getMaxExcitationEnergy();
                if (daughter->getName()=="Li6")
                    CHECK(ex1>=exmin && ex1<=exmax);
            } 
        }
        // sample again
        exs.sampleExcitationEnergies();
        // check that new value is generated which is different from the previous
        double ex2=0;
        for (auto& p : *chain) { 
            auto& d = p.getDaughters();
            for (auto& daughter : d) { // loop over daughters
                ex2 = daughter -> getExcitationEnergy();
                if (daughter->getName()=="Li6")
                    CHECK(abs(ex2-ex1)>0.1);
            } 
        }
    }   
*/

    TEST_FIXTURE(SetupFixture, CheckDepthIs3) {
        auto& t = chain->getTree();

        CHECK_EQUAL(0, ProcessChain::Tree::depth(t.begin())); // Root is at depth 0

        // Hence depth is 0 based
        // The tree is then [compound formation] -> [compound decay] -> [C12 decay]
        CHECK_EQUAL(2, t.max_depth());

        // Iterator to compound decay
        auto b = t.begin_fixed(t.begin(), 1); // t.begin() given iterator to root
           //b = t.begin(1)                   // shortcut that assumes you want iterator from root
        CHECK(!b.isLeaf());

        /* How one would normally do, but apparent this is broke */
        //auto e = t.end_fixed(t.begin(), 1);

        // Does the iterator point to something meaningful
        CHECK(t.is_valid(b));

        int count = 0;
        while (t.is_valid(b)) {
            ProcessChain::DecayPtr& p = *b;

            p->getDaughters();

            ++count;
            ++b; // Increment to next element
        }

        // There is only 1 item @ depth 1
        CHECK_EQUAL(1, count);

        b = t.begin_fixed(2);
        CHECK(b.isLeaf());
    }

    TEST_FIXTURE(SetupFixture, CheckPossibleToIterateBothWays) {
        auto& t = chain->getTree();

        auto N = t.max_depth();

        // Bottom up
        for (auto n = N; n >= 0; --n) {
            auto i = t.begin_fixed(n);

            // Loop from left to right
            while (t.is_valid(i)) {
                ProcessChain::DecayPtr& proc = *i;

                proc->getDaughters();

                ++i;
            }
        }

        // Top down
        for (auto n = 0; n <= N; n++) {
            auto i = t.begin_fixed(n);

            while (t.is_valid(i)) {
                ProcessChain::DecayPtr& proc = *i;

                proc->getDaughters();

                ++i;
            }
        }
    }
    
}
