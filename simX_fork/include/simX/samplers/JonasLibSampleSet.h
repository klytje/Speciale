//
// Created by munk on 6/19/16.
//

#ifndef SIMX_SAMPLERS_JONASLIBSAMPLESET_H
#define SIMX_SAMPLERS_JONASLIBSAMPLESET_H

#include "SampleSet.h"
#include <string>
#include <memory>
#include <TVector3.h>

class TChain;
class TClonesArray;

namespace simX {
    namespace samplers {

        /**
         * A specilization of SampleSet, where the samples from Jonas are stored in a ROOT file.
         */
        class JonasLibSampleSet : public SampleSet {
        public:

            /**
             * Creates a RootInput from the path.
             *
             * @param path The path of the root file.
             * @param nParticles Number of expected particles
             * @throws invalid_argument If file cannot be found or branch is missing.
             */
            JonasLibSampleSet(std::string file, size_t nParticles);
            virtual ~JonasLibSampleSet();

            virtual const std::vector<double> &getSample() override;

            virtual const double getWeight() override;

            virtual bool isSampled() const override;

            virtual void next() override;

            virtual bool hasNext() override;

            virtual void reset() override;

        private:
            std::unique_ptr<TChain> tree;
            ssize_t current;

            TClonesArray *arr;
            std::vector<double> result;
            double w;
        };
    }
}
#endif //SIMX_SAMPLERS_JONASLIBSAMPLESET_H
