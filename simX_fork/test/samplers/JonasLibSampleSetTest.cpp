
#include <unittest++/UnitTest++.h>
#include <simX/samplers/JonasLibSampleSet.h>


using namespace std;
using namespace simX;
using namespace simX::samplers;

SUITE(JonasLibSampleSetTest) {

    class SetupFixture {
    public:
        SetupFixture()
                : sampler("test/_res/samplers/jonas_lib_sample.root", 3)
        {
        }
        JonasLibSampleSet sampler;
    };
    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1,1);
    }


    TEST_FIXTURE(SetupFixture, SamplerCanLoadFirstEntry) {
        sampler.next();
    }

    TEST_FIXTURE(SetupFixture, JonasLibFilesAreNotSampled) {
        CHECK(!sampler.isSampled());
    }

    TEST_FIXTURE(SetupFixture, FileHasTenEntries) {
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();

        CHECK(!sampler.hasNext());
    }

    TEST_FIXTURE(SetupFixture, FirstEntryFirstVectorMatches) {
        sampler.next();

        double expected[] = {-109451.0, 77442.924, -50395.46,
                             109315.16, -75080.33, 49221.841,
                             135.86813, -2362.586, 1173.6200};
        CHECK_ARRAY_CLOSE(expected, sampler.getSample(), 9, 1E-1);
    }

    TEST_FIXTURE(SetupFixture, LastEntryFirstVectorMatches) {
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();
        sampler.next();

        double expected[] = {108465.50, 63809.016, -87882.08,
                             -59193.27, -74007.09, 15681.589,
                             -49272.22, 10198.081, 72200.496};
        CHECK_ARRAY_CLOSE(expected, sampler.getSample(), 9, 1E-1);
    }

    TEST(WrongNumberOfParticlesIsDetected) {
        CHECK_THROW(JonasLibSampleSet("test/_res/samplers/jonas_lib_sample.root", 2), std::runtime_error);
    }

}
