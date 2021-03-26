#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <ausa/setup/Setup.h>
#include <ausa/setup/RoundDSSD.h>
#include <TVector3.h>

#include "simX/Particle.h"
#include "simX/Detection/Detector.h"
#include "simX/Detection/SegmentedDetector.h"
#include <simX/Detection/S3.h>
#include <simX/parser/DetectionSystemParser.h>
#include "simX/Layer.h"
#include "simX/Random.h"
#include <simX/propagator/ParticlePropagator.h>
#include <simX/propagator/IonizingPropagator.h>

#include <ausa/eloss/SRIMTabulation.h>
#include <ausa/eloss/RangeInterpolator.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>
#include <ausa/util/Resource.h>
#include <ausa/util/memory>

#include <math.h>
#include <TMath.h>
#include <tuple>
#include <iostream>
#include <ausa/json/IO.h>

using namespace simX;
using namespace AUSA;
using namespace AUSA::EnergyLoss;
using namespace simX::propagator;
using namespace simX::parser;
using namespace simX::detection;

SUITE(S3SetupTest) {

    // use SRIM stopping powers for alpha in Silicon
    IonizingPropagator::LossCalc factory(const Layer &layer, const Particle &p) {
        std::string suffix = "/SRIM08/He4_Si.dat";
        SRIMTabulation tab(AUSA::getResourceDirectory() + suffix);
        return std::make_unique<EnergyLossRangeInverter>(std::make_unique<RangeInterpolator>(tab));
    }

    class SetupFixture {
    public:
        SetupFixture()
                : alpha(Ion(2, 4)),
                  ionProp(std::make_shared<IonizingPropagator>(IonizingPropagator(factory))),
                  setup(AUSA::JSON::readSetupFromJSON("test/_res/hit_setup/setup.json")) {
            alpha.setPropagator(ionProp);
            aSB = std::dynamic_pointer_cast<RoundDSSD>(setup->getDSSD("SB"));
            aST = std::dynamic_pointer_cast<RoundDSSD>(setup->getDSSD("ST"));

            DetectionSystemParser parser{false};
            detSys = parser.buildSimXSetup(*setup);

            for (auto &d : detSys.getDetectors()) {
                if (d->getName() == "SB") sSB = std::dynamic_pointer_cast<S3>(d);
                if (d->getName() == "ST") sST = std::dynamic_pointer_cast<S3>(d);
            }
        }

        Particle alpha;
        std::shared_ptr<IonizingPropagator> ionProp;
        std::shared_ptr<Setup> setup;
        std::shared_ptr<RoundDSSD> aSB, aST;
        DetectionSystem detSys;
        std::shared_ptr<S3> sSB, sST;
    };


    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1, 1);

        CHECK(sST);
        CHECK(aST);
        CHECK(sSB);
        CHECK(aSB);
    }

    TEST_FIXTURE(SetupFixture, TestSimX_SB_IsEqualToAUSAlibVersion) {
        CHECK_EQUAL(24, aSB->backStripCount());
        CHECK_EQUAL(24, sSB->getNRings());

        CHECK_EQUAL(32, sSB->spokeCount());
        CHECK_EQUAL(32, sSB->getNSpokes());

        CHECK_CLOSE(1038E-3, aSB->getThickness(), 1E-6);
//        CHECK_CLOSE(1038E-3, sSB->getThickness(), 1E-6);

        CHECK_CLOSE(0.5, aSB->getRingPitch(), 1E-3);
        CHECK_CLOSE(0.5, aSB->getRingWidth(), 1E-3);

        CHECK_CLOSE(0.5, sSB->getRingWidth(), 1E-3);
        CHECK_CLOSE(0.0, sSB->getRingGap(), 1E-3);

        CHECK_CLOSE(11, aSB->getInnerRadius(), 1E-3);
        CHECK_CLOSE(11, sSB->getInnerRadius(), 1E-3);
    }



    TEST_FIXTURE(SetupFixture, MiddleOfSBPixel11IsReportedAsPixel11) {
        auto p = aSB->getPixelPosition(1, 1);

        auto n = aSB->getNormal();
        alpha.setPosition(p + 100 * n);
        alpha.setFourMomentumLab(5E3, -n);

        CHECK(sSB->getDetectorVolume().isIntersecting(alpha.getPosition(), alpha.getDirectionLab()));

        auto out = sSB->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(0, out[0].channel);
        CHECK_EQUAL(32, out[1].channel);
    }

    TEST_FIXTURE(SetupFixture, ItIsPossibleToHitPixbel11From000) {
        auto p = aSB->getPixelPosition(1, 1);

        alpha.setFourMomentumLab(5E3, p.Unit());

        auto out = sSB->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(0, out[0].channel);
        CHECK_EQUAL(32, out[1].channel);
    }

    TEST_FIXTURE(SetupFixture, SBRingHalfShouldBeRadius) {
        auto p = aSB->getContinuousPixelPosition(1, 0.5001);

        auto n = aSB->getNormal();
        alpha.setPosition(p+100*n);
        alpha.setFourMomentumLab(5E3, -n);

        auto out = sSB->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(0, out[0].channel);
        CHECK_EQUAL(32, out[1].channel);
    }

    TEST_FIXTURE(SetupFixture, STRingHalfShouldBeRadius) {
        auto p = aST->getContinuousPixelPosition(1, 0.5 + 1E-6);

        auto n = aST->getNormal();
        alpha.setPosition(p+100*n);
        alpha.setFourMomentumLab(5E3, -n);

        auto out = sST->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(0, out[0].channel);
        CHECK_EQUAL(32, out[1].channel);
    }

    TEST_FIXTURE(SetupFixture, CenterOfSBPixel2_23_CanBeHit) {
        auto p = aST->getPixelPosition(2, 23);

        auto n = aST->getNormal();
        alpha.setPosition(p+100*n);
        alpha.setFourMomentumLab(5E3, -n);

        auto out = sST->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(1, out[0].channel);
        CHECK_EQUAL(31+23, out[1].channel);
    }

    TEST_FIXTURE(SetupFixture, CenterOfSBPixel6_23_CanBeHit) {
        auto p = aST->getPixelPosition(6, 23);

        auto n = aST->getNormal();
        alpha.setPosition(p+100*n);
        alpha.setFourMomentumLab(5E3, -n);

        auto out = sST->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(5, out[0].channel);
        CHECK_EQUAL(31+23, out[1].channel);
    }

    TEST_FIXTURE(SetupFixture, CenterOfSBPixel7_23_CanBeHit) {
        auto p = aST->getPixelPosition(7, 23);

        auto n = aST->getNormal();
        alpha.setPosition(p+100*n);
        alpha.setFourMomentumLab(5E3, -n);

        auto out = sST->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(6, out[0].channel);
        CHECK_EQUAL(31+23, out[1].channel);
    }

    TEST_FIXTURE(SetupFixture, CenterOfSBPixel8_23_CanBeHit) {
        auto p = aST->getPixelPosition(8, 23);

        auto n = aST->getNormal();
        alpha.setPosition(p+100*n);
        alpha.setFourMomentumLab(5E3, -n);

        auto out = sST->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(7, out[0].channel);
        CHECK_EQUAL(31+23, out[1].channel);
    }


    TEST_FIXTURE(SetupFixture, CenterOfSBPixel32_23_CanBeHit) {
        auto p = aST->getPixelPosition(32, 23);

        auto n = aST->getNormal();
        alpha.setPosition(p+100*n);
        alpha.setFourMomentumLab(5E3, -n);

        auto out = sST->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(31, out[0].channel);
        CHECK_EQUAL(31+23, out[1].channel);
    }

    TEST_FIXTURE(SetupFixture, EdgeOfPixel32_23_CanBeHit) {
        auto p = aST->getContinuousPixelPosition(32.5-1E-6, 23);

        auto n = aST->getNormal();
        alpha.setPosition(p+100*n);
        alpha.setFourMomentumLab(5E3, -n);

        auto out = sST->detect(alpha);

        CHECK_EQUAL(2, out.size());
        CHECK_EQUAL(31, out[0].channel);
        CHECK_EQUAL(31+23, out[1].channel);
    }
}
