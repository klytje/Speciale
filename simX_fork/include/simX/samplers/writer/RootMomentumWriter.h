//
// Created by munk on 6/22/16.
//

#ifndef SIMX_ROOTMOMENTUMWRITER_H
#define SIMX_ROOTMOMENTUMWRITER_H

#include "MomentumSampleWrite.h"

class TFile;
class TTree;

namespace simX {
    namespace samplers {
        namespace writer {
            class RootMomentumWriter : public MomentumSampleWriter {
            public:
                RootMomentumWriter(std::string path, bool hasWeight, size_t nElements);

                virtual ~RootMomentumWriter();

                virtual void addSample(const std::vector<TVector3> &sample, double w) override;

            private:
                TFile* file;
                TTree* tree;

                std::vector<double> sample;
                double weight;
            };
        }
    }
}
#endif //SIMX_ROOTMOMENTUMWRITER_H
