//
// Created by munk on 6/16/16.
//

#ifndef SIMX_TEXTINPUT_H_H
#define SIMX_TEXTINPUT_H_H

#include "SampleSet.h"
#include <string>
#include <fstream>

namespace simX {
    namespace samplers {

        /**
         * A specilization of SampleSet, where the samples are stored in a text file.
         *
         * Syntax:
         * E0 E1 ... En [weight]
         * ie. N elements with an optional weight in the end.
         *
         * This input is considered sampled if the weight is not provided.
         *
         */
        class TextSampleSet : public SampleSet {
        public:

            /**
             * Create a TextInput from a file
             *
             * @param path Path of input file
             * @param elements Number of particles
             */
            TextSampleSet(std::string path, size_t elements);

            virtual const std::vector<double> &getSample() override;

            virtual const double getWeight() override;

            virtual bool isSampled() const override;

            virtual void next() override;

            virtual bool hasNext() override;

            virtual void reset() override;


        private:
            std::fstream stream;
            bool sampled;
            std::vector<double> sample;
            std::string line, path;
            size_t elements, current;
            double weight;
            bool read;

            void readSample();

        };
    }
}

#endif //SIMX_TEXTINPUT_H_H
