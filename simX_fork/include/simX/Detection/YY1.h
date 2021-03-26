//
// Created by munk on 27-03-17.
//

#ifndef SIMX_YY1_H
#define SIMX_YY1_H

#include "SegmentedDetector.h"

namespace simX {
    namespace Detection {
        class YY1 : public SegmentedDetector {
        public:
            YY1(const std::string& name, const TVector3& position, const TVector3& normal, const TVector3& orientation,
                double activeThickness, double dlFront, double dlBack);

        protected:
            bool strikesGrid(double u) override;

            int gridIndex() override;

            struct Args;

        };
    }
}

#endif //SIMX_YY1_H
