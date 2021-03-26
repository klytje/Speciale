
#include <unittest++/UnitTest++.h>
#include <simX/samplers/SampleSampler.h>
#include "DummySampler.h"
#include <simX/Particle.h>
#include <ausa/util/memory>
#include <simX/samplers/MultiSampler.h>


using namespace std;
using namespace simX;
using namespace simX::samplers;

SUITE(MultiSamplerTest) {

    initializer_list<initializer_list<double>> init(initializer_list<initializer_list<double>> list) {
        return list;
    }

    TEST(Sanity) {
        CHECK_EQUAL(1,1);
    }


    TEST(IfOnlySingleElementThatIsReturned) {
        auto d = std::make_unique<DummySampler>(init({{1,1,1}}));
        MultiSampler::SampleCollection samples;
        samples.push_back(move(d));
        MultiSampler multi{move(samples)};

        CHECK(multi.hasNext());

        multi.next();

        auto& sample = multi.getSample();
        CHECK_ARRAY_EQUAL(vector<double>({1, 1, 1}), sample, 3);
        CHECK(!multi.hasNext());
    }

    TEST(IfTwoElementsProvidedTheyWillBeSampled) {
        auto d = std::make_unique<DummySampler>(init({{1,1,1}, {2,2,2}}));
        MultiSampler::SampleCollection samples;
        samples.push_back(move(d));
        MultiSampler multi{move(samples)};

        CHECK(multi.hasNext());

        multi.next();
        auto& sample = multi.getSample();
        CHECK_ARRAY_EQUAL(vector<double>({1, 1, 1}), sample, 3);

        multi.next();

        auto& sample2 = multi.getSample();
        CHECK_ARRAY_EQUAL(vector<double>({2, 2, 2}), sample2, 3);
        CHECK(!multi.hasNext());
    }

    TEST(IfMultipleSamplesProvidedSampleAll) {
        auto d1 = std::make_unique<DummySampler>(init({{1,1,1}}));
        auto d2 = std::make_unique<DummySampler>(init({{9,9,9}}));
        MultiSampler::SampleCollection samples;
        samples.push_back(move(d1));
        samples.push_back(move(d2));
        MultiSampler multi{move(samples)};

        CHECK(multi.hasNext());
        multi.next();
        CHECK_ARRAY_EQUAL(vector<double>({1, 1, 1}), multi.getSample(), 3);

        CHECK(multi.hasNext());
        multi.next();
        CHECK_ARRAY_EQUAL(vector<double>({9, 9, 9}), multi.getSample(), 3);

        CHECK(!multi.hasNext());
    }

    TEST(SampleSkipsMiddleEmptySets) {
        auto d1 = std::make_unique<DummySampler>(init({{1,1,1}}));
        auto d2 = std::make_unique<DummySampler>(init({}));
        auto d3 = std::make_unique<DummySampler>(init({{9,9,9}}));
        MultiSampler::SampleCollection samples;
        samples.push_back(move(d1));
        samples.push_back(move(d2));
        samples.push_back(move(d3));
        MultiSampler multi{move(samples)};

        CHECK(multi.hasNext());
        multi.next();
        CHECK_ARRAY_EQUAL(vector<double>({1, 1, 1}), multi.getSample(), 3);

        CHECK(multi.hasNext());
        multi.next();
        CHECK_ARRAY_EQUAL(vector<double>({9, 9, 9}), multi.getSample(), 3);

        CHECK(!multi.hasNext());
    }

    TEST(SampleSkipsInitialEmptySets) {
        auto d1 = std::make_unique<DummySampler>(init({{1,1,1}}));
        auto d2 = std::make_unique<DummySampler>(init({}));
        auto d3 = std::make_unique<DummySampler>(init({{9,9,9}}));
        MultiSampler::SampleCollection samples;
        samples.push_back(move(d2));
        samples.push_back(move(d1));
        samples.push_back(move(d3));
        MultiSampler multi{move(samples)};

        CHECK(multi.hasNext());
        multi.next();
        CHECK_ARRAY_EQUAL(vector<double>({1, 1, 1}), multi.getSample(), 3);

        CHECK(multi.hasNext());
        multi.next();
        CHECK_ARRAY_EQUAL(vector<double>({9, 9, 9}), multi.getSample(), 3);

        CHECK(!multi.hasNext());
    }
}
