
#include <iostream>
#include "AntiProtonCutter.h"
#include "ausa/identify/Identification.h"
#include "ausa/identify/Hit.h"
#include "ausa/eloss/Default.h"
#include "ausa/constants/Mass.h"

using namespace std;
using namespace AUSA;
using namespace AUSA::Event;

namespace {
    auto const ALPHA = EnergyLoss::Ion{"He4"};
}

AntiProtonCutter::AntiProtonCutter(double cutoff)
 : AbstractEventCutter("AntiProtonCutter"), cutoff(cutoff),
   accepted(0), rejected(0)
{
    calc = EnergyLoss::defaultRangeInverter("p", "Silicon");


}


bool AntiProtonCutter::cut(const PhysicsEvent& event) const
{
    for (size_t i = 0; i < event.getDetectionMultiplicity(); ++i) {
        auto id = event.getDetectionEvent()->getIdentification(i);
        auto hit = id->getDssdHit();

        if (hit->getDetectorIndex() == 1
                && id->getKineticEnergy() < 400
                && id->getDeadLayerLoss() < 220
                && id->getParticleType()->ion == ALPHA
           ) {
            rejected++;
            return false;
        }
    }
    accepted++;
    return true;
}
