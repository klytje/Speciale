
#include <unittest++/UnitTest++.h>
#include <simX/samplers/RootSampleSet.h>
#include <simX/Particle.h>


using namespace std;
using namespace simX;
using namespace simX::samplers;

SUITE(RootSampleSetTest) {

    TEST(Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST(SimpleFileParses) {
        RootSampleSet{"test/_res/samplers/sample_with_weight.root", 6};
    }

    TEST(FirstWeightIs10) {
        RootSampleSet input{"test/_res/samplers/sample_with_weight.root", 6};
        input.next();

        CHECK_EQUAL(10, input.getWeight());
    }

    TEST(SecondWeightIs9) {
        RootSampleSet input{"test/_res/samplers/sample_with_weight.root", 6};
        input.next();
        input.next();

        CHECK_EQUAL(9, input.getWeight());
    }

    TEST(SimpleInputIsNotSampled) {
        RootSampleSet input{"test/_res/samplers/sample_with_weight.root", 6};

        CHECK(!input.isSampled());
    }

    TEST(IfInputHas3TimesParticleColumnsItIsSampled) {
        RootSampleSet input{"test/_res/samplers/sample_without_weight.root", 6};

        CHECK(input.isSampled());
    }

    TEST(FirstSampleIs000000) {
        RootSampleSet input{"test/_res/samplers/sample_without_weight.root", 6};
        input.next();

        auto& sample = input.getSample();
        vector<double> expected = {0,0,0, 0,0,0};

        CHECK_ARRAY_EQUAL(expected, sample, 6);
    }

    TEST(SecondSampleIs0_1E3_0_00_1E3) {
        RootSampleSet input{"test/_res/samplers/sample_without_weight.root", 6};
        input.next();
        input.next();

        auto& sample = input.getSample();
        vector<double> expected = {0,1E3,0, 0,0,1E3};

        CHECK_ARRAY_EQUAL(expected, sample, 6);
    }

    TEST(HasNextDoesNotModifySample) {
        RootSampleSet input{"test/_res/samplers/sample_without_weight.root", 6};

        auto& sample = input.getSample();
        vector<double> expected = {0,0,0, 0,0,0};

        CHECK_ARRAY_EQUAL(expected, sample, 6);

        input.hasNext();
        auto& sample6 = input.getSample();
        CHECK_ARRAY_EQUAL(expected, sample6, 6);
    }

    TEST(HasNextIndicatesWhetherMoreInputAvaiable) {
        RootSampleSet input{"test/_res/samplers/sample_without_weight.root", 6};

        CHECK(input.hasNext());
        input.next();
        CHECK(input.hasNext());
        input.next();
        CHECK(!input.hasNext());
    }

    TEST(ResetResets) {
        RootSampleSet input{"test/_res/samplers/sample_without_weight.root", 6};

        input.next();
        input.next();
        CHECK(!input.hasNext());
        input.reset();
        CHECK(input.hasNext());
    }
}
