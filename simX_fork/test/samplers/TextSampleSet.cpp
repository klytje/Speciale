
#include <unittest++/UnitTest++.h>
#include <simX/samplers/TextSampleSet.h>
#include <simX/Particle.h>


using namespace std;
using namespace simX;
using namespace simX::samplers;

SUITE(TextSampleSetTest) {

    TEST(Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST(SimpleFileParses) {
        TextSampleSet{"test/_res/samplers/input_with_weight.txt", 6};
    }

    TEST(FirstWeightIs10) {
        TextSampleSet input{"test/_res/samplers/input_with_weight.txt", 6};
        input.next();

        CHECK_EQUAL(10, input.getWeight());
    }

    TEST(SecondWeightIs20) {
        TextSampleSet input{"test/_res/samplers/input_with_weight.txt", 6};
        input.next();
        input.next();

        CHECK_EQUAL(9, input.getWeight());
    }

    TEST(SimpleInputIsNotSampled) {
        TextSampleSet input{"test/_res/samplers/input_with_weight.txt", 6};
        input.next();

        CHECK(!input.isSampled());
    }

    TEST(IfInputHas3TimesParticleColumnsItIsSampled) {
        TextSampleSet input{"test/_res/samplers/input_without_weight.txt", 6};
        input.next();

        CHECK(input.isSampled());
    }

    TEST(FirstSampleIs163123) {
        TextSampleSet input{"test/_res/samplers/input_with_weight.txt", 6};
        input.next();

        auto& sample = input.getSample();
        vector<double> expected = {1,2,3, 1,2,3};

        CHECK_ARRAY_EQUAL(expected, sample, 6);
    }

    TEST(SecondSampleIs456456) {
        TextSampleSet input{"test/_res/samplers/input_with_weight.txt", 6};
        input.next();
        input.next();

        auto& sample = input.getSample();
        vector<double> expected = {4,5,6,4,5,6};

        CHECK_ARRAY_EQUAL(expected, sample, 6);
    }

    TEST(HasNextDoesNotModifySample) {
        TextSampleSet input{"test/_res/samplers/input_with_weight.txt", 6};

        auto& sample = input.getSample();
        vector<double> expected = {1,2,3, 1,2,3};

        CHECK_ARRAY_EQUAL(expected, sample, 6);

        input.hasNext();
        auto& sample2 = input.getSample();
        CHECK_ARRAY_EQUAL(expected, sample2, 6);
    }

    TEST(HasNextIndicatesWhetherMoreInputAvaiable) {
        TextSampleSet input{"test/_res/samplers/input_with_weight.txt", 6};

        CHECK(input.hasNext());
        input.next();
        CHECK(input.hasNext());
        input.next();
        CHECK(!input.hasNext());
    }

    TEST(FlowerDataIsNotSampled) {
        TextSampleSet input{"test/_res/samplers/flower.dat", 9};
        CHECK(!input.isSampled());
    }

    TEST(FirstFlowerLineIsCorrect) {
        TextSampleSet input{"test/_res/samplers/flower.dat", 9};

        double expected[] = {
                16.2407, -42.9347, -19.4811, 77.1825, -44.9068, 91.8585, -93.4232, 87.8415, -72.3775,
        };

        input.next();
        CHECK_ARRAY_CLOSE(expected, input.getSample(), 9, 1E-3);
        CHECK_CLOSE(3.55877e-08, input.getWeight(), 1E-11);


        input.next();
        double expected2[] = {
                74.3061, -57.1667, 1.47314, 16.1366, -54.7682, 0.504639, -90.4427, 111.935, 3.42261
        };
        CHECK_ARRAY_CLOSE(expected2, input.getSample(), 9, 1E-3);
        CHECK_CLOSE(2, input.getWeight(), 1E-11);
    }

    TEST(SecondFlowerLineIsCorrect) {
        TextSampleSet input{"test/_res/samplers/flower.dat", 9};

        double expected[] = {
                74.3061, -57.1667, 1.47314, 16.1366, -54.7682, 0.504639, -90.4427, 111.935, 3.42261
        };

        input.next();
        input.next();
        CHECK_ARRAY_CLOSE(expected, input.getSample(), 9, 1E-3);
        CHECK_CLOSE(2, input.getWeight(), 1E-11);
    }

    TEST(FlowerLine5000IsCorrect) {
        TextSampleSet input{"test/_res/samplers/flower.dat", 9};

        double expected[] = {
                -62.3323, 126.653, 2.86377, 28.8217, -14.8525, 0.656982, 33.5106, -111.801, 1.87964
        };

        for (size_t i = 0; i < 5000; ++i) {
            input.next();
        }

        CHECK_ARRAY_CLOSE(expected, input.getSample(), 9, 1E-3);
        CHECK_CLOSE(5357, input.getWeight(), 1E-11);
    }
}
