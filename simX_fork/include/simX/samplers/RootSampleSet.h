//
// Created by munk on 6/19/16.
//

#ifndef SIMX_SAMPLERS_ROOTINPUT_H
#define SIMX_SAMPLERS_ROOTINPUT_H

#include "SampleSet.h"
#include <string>

class TTree;
class TFile;

namespace simX {
    namespace samplers {

        /**
         * A specilization of SampleSet, where the samples are stored in a ROOT file.
         *
         * The ROOT file for N elements must contain a TTree named 'sample' containing
         *   1. A TParameter<int> named N_ELM which states the number of elements (should be N)
         *   2. A branch named p, which is a double array of size N
         *   3. Optional. A branch named w containing the weight. This must be a double.
         */
        class RootSampleSet : public SampleSet {
        public:

            /**
             * Creates a RootInput from the path.
             *
             * @param path The path of the root file.
             * @param nElements Number of expected particles
             * @throws invalid_argument If file cannot be found or branch is missing.
             */
            RootSampleSet(std::string file, size_t nElements);
            virtual ~RootSampleSet();

            virtual const std::vector<double> &getSample() override;

            virtual const double getWeight() override;

            virtual bool isSampled() const override;

            virtual void next() override;

            virtual bool hasNext() override;

            virtual void reset() override;

        private:
            TFile* file;
            TTree* tree;
            ssize_t current;

            std::vector<double> result;
            double w;
            bool sampled;
        };
    }
}
#endif //SIMX_SAMPLERS_ROOTINPUT_H
