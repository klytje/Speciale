//
// Created by munk on 11-10-16.
//

#ifndef EVENTANALYSIS_MULSPECTRUM_H
#define EVENTANALYSIS_MULSPECTRUM_H

#include <ausa/event/id_cut/AbstractIdCutter.h>
#include <TH1I.h>

namespace AUSA {
    namespace Event {
        class MulSpectrum : public AbstractIdCutter{
        public:
            MulSpectrum();
            virtual std::vector<bool> cut(const IdVector &ids) const override;

            mutable TH1I mul;
        };
    }
}

#endif //EVENTANALYSIS_MULSPECTRUM_H
