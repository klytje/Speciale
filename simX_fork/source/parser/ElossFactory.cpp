#include <simX/propagator/NonIonizingPropagator.h>
#include <simX/propagator/NoTrackPropagator.h>
#include <ausa/eloss/EnergyLossCalculator.h>
#include <simX/propagator/IonizingPropagator.h>
#include <simX/propagator/GaussianStragglingPropagator.h>
#include <simX/propagator/MCStragglingPropagator.h>
#include "simX/parser/ElossFactory.h"

#include <ausa/eloss/Ion.h>
#include <ausa/eloss/Default.h>
#include <ausa/json/JSONUtil.h>

#include <ausa/util/Resource.h>
#include <ausa/util/memory>
#include "ausa/eloss/EnergyLossIntegrator.h"
#include "ausa/eloss/SRIMTabulation.h"
#include "ausa/eloss/SRIM13Tabulation.h"
#include "ausa/eloss/StoppingPowerInterpolator.h"
#include "ausa/eloss/RangeSplineFitter.h"
#include "ausa/eloss/RangeInterpolator.h"
#include "ausa/eloss/Material.h"

namespace simX {
    namespace parser {
        Config::PropagatorFactory noLossFactory(ConfigParser::Options& o) {
            return [](simX::Particle& p) -> simX::Particle::Propagator {
                if (p.getTracking())
                    return std::make_shared<propagator::NonIonizingPropagator>();
                else
                    return std::make_shared<propagator::NoTrackPropagator>();
            };
        }

        Config::PropagatorFactory elossFactory(ConfigParser::Options& o) {
            auto fInner = elossCalculator(o);
            auto nonIon = noLossFactory(o);
            return [=](simX::Particle& p) -> simX::Particle::Propagator {
                if (!p.getTracking())
                    return std::make_shared<propagator::NoTrackPropagator>();
                else if (p.getZ() > 0) {
                    return std::make_shared<propagator::IonizingPropagator>(fInner);
                }
                else
                    return nonIon(p);
            };
        }

        Config::PropagatorFactory gaussianStragglingFactory(ConfigParser::Options& o) {
            auto fInner = elossCalculator(o);
            auto nonIon = noLossFactory(o);
            return [=](simX::Particle& p) -> simX::Particle::Propagator {
                if (!p.getTracking())
                    return std::make_shared<propagator::NoTrackPropagator>();
                else if (p.getZ() > 0) {
                    auto ionprop = std::make_shared<propagator::IonizingPropagator>(fInner);
                    return std::make_shared<propagator::GaussianStragglingPropagator>(ionprop);
                }
                else
                    return nonIon(p);
            };
        }

        Config::PropagatorFactory mcStragglingFactory(ConfigParser::Options& o) {
            auto mul = AUSA::JSON::readDoubleOrDefault(o, "multiplier", 0.04);

            auto nonIon = noLossFactory(o);
            return [=](simX::Particle& p) -> simX::Particle::Propagator {
                if (!p.getTracking())
                    return std::make_shared<propagator::NoTrackPropagator>();
                else if (p.getZ() > 0) {
                    auto prop = std::make_shared<propagator::MCStragglingPropagator>(lossFactory, rangeFactory);
                    prop->setCrossSectionMultiplier(mul);
                    return prop;
                }
                else
                    return nonIon(p);
            };
        }

        simX::propagator::IonizingPropagator::CalculatorFactory elossCalculator(ConfigParser::Options& o) {
            std::string tab = AUSA::JSON::readStringOrDefault(o, "tabulation", "GEANT");
            return [tab](const simX::Layer& layer, const simX::Particle& p) -> std::unique_ptr<AUSA::EnergyLoss::EnergyLossCalculator> {
                if(tab == "SRIM13") return defaultRangeInverter(AUSA::EnergyLoss::Ion{p.getZ(), p.getA()}, layer.getMaterial());
                else if(tab == "GEANT") return eLossDefault(AUSA::EnergyLoss::Ion{p.getZ(), p.getA()}, layer.getMaterial());
                else throw std::invalid_argument("'"  + tab + "' is not a valid tabulation!");
            };
        }

        simX::propagator::MCStragglingPropagator::LossCalc lossFactory(const Layer& layer, const Particle& p)
        {
            AUSA::EnergyLoss::SRIM13Tabulation tabulation{p.getZ(), p.getA(), layer.getMaterial()};
            auto interpolator = std::make_unique<AUSA::EnergyLoss::StoppingPowerInterpolator>(tabulation);
            return std::make_unique<AUSA::EnergyLoss::EnergyLossIntegrator>(std::move(interpolator));
        }

        simX::propagator::MCStragglingPropagator::RangeCalc rangeFactory(const Layer& layer, const Particle& p)
        {
            AUSA::EnergyLoss::SRIM13Tabulation tabulation{p.getZ(), p.getA(), layer.getMaterial()};
            return std::make_unique<AUSA::EnergyLoss::RangeSplineFitter>(tabulation);
        }

    }
}
