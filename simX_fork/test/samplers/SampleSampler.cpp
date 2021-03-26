
#include <unittest++/UnitTest++.h>
#include <simX/samplers/SampleSampler.h>
#include "DummySampler.h"
#include <simX/Particle.h>
#include <ausa/util/memory>


using namespace std;
using namespace simX;
using namespace simX::samplers;

SUITE(SampleSamplerTest) {

    initializer_list<initializer_list<double>> init(initializer_list<initializer_list<double>> list) {
        return list;
    }

    TEST(Sanity) {
        CHECK_EQUAL(1,1);
    }


    TEST(IfOnlySingleElementThatIsReturned) {
        auto d = std::make_unique<DummySampler>(init({{1,1,1}}));
        SampleSampler sampler(move(d));

        sampler.next();

        auto& sample = sampler.getSample();
        CHECK_ARRAY_EQUAL(vector<double>({1, 1, 1}), sample, 3);
    }

    TEST(IfTwoElementsProvidedTheyWillBeSampled) {
        auto d = std::make_unique<DummySampler>(init({{1,1,1}, {2,2,2}}));
        d->weight = {1,1}; // Equal weights
        SampleSampler sampler(move(d));


        sampler.next();
        auto sample = sampler.getSample();
        CHECK_ARRAY_EQUAL(vector<double>({2, 2, 2}), sample, 3);

        sampler.next();
        sample = sampler.getSample();
        CHECK_ARRAY_EQUAL(vector<double>({1, 1, 1}), sample, 3);

        sampler.next();
        sample = sampler.getSample();
        CHECK_ARRAY_EQUAL(vector<double>({1, 1, 1}), sample, 3);
    }
}
