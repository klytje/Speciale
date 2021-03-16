
#ifndef ANTI_PROTON_CUTTER
#define ANTI_PROTON_CUTTER

#include <memory>
#include <TH1I.h>
#include "ausa/event/cut/AbstractEventCutter.h"
#include "ausa/eloss/EnergyLossCalculator.h"

namespace AUSA {
    namespace Event {
        class AntiProtonCutter : public AbstractEventCutter {

        public:
            AntiProtonCutter(double cutoff);

            virtual ~AntiProtonCutter() = default;

            virtual bool cut(const PhysicsEvent& physEvent) const override;

            mutable size_t accepted, rejected;

        private:
            double cutoff;
            std::unique_ptr<EnergyLoss::EnergyLossCalculator> calc;
        };
    }
}

#endif // ANTI_PROTON_CUTTER
