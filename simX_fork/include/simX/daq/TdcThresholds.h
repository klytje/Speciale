//
// Created by munk on 09-11-16.
//

#ifndef SIMX_TDCTHRESHOLDS_H
#define SIMX_TDCTHRESHOLDS_H

#include <vector>
#include <TFile.h>
#include <memory>
#include <TH1D.h>
#include "simX/Random.h"

namespace simX {
    namespace daq {
        class TdcThresholds {
        public:
            TdcThresholds(std::vector<TH1D> hists) : hists(hists) {};

            bool isTriggering(int channel, double E) const {
                auto& h = hists[channel];
                auto n = h.FindFixBin(E);
                if (n >= h.GetNbinsX()) return true;

                auto y = h.GetBinContent(n);
                return y == 1 || rnd() >= y;
            }

        private:
            std::shared_ptr<TFile> file;
            std::vector<TH1D> hists;
        };
    }
}
#endif //SIMX_TDCTHRESHOLDS_H
