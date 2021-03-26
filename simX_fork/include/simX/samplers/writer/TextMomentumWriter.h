//
// Created by munk on 6/22/16.
//

#ifndef SIMX_TEXTMOMENTUMWRITER_H
#define SIMX_TEXTMOMENTUMWRITER_H

#include <fstream>
#include "MomentumSampleWrite.h"

namespace simX {
    namespace samplers {
        namespace writer {
            class TextMomentumWriter : public MomentumSampleWriter {
            public:
                TextMomentumWriter(std::string path, bool weight);

                virtual void addSample(const std::vector<TVector3> &sample, double w) override;

            private:
                std::ofstream stream;
                bool weight;
            };
        }
    }
}
#endif //SIMX_TEXTMOMENTUMWRITER_H
