
#include <unittest++/UnitTest++.h>
#include <ausa/eloss/Ion.h>
#include <TVector3.h>
#include "simX/Particle.h"
#include "simX/Detection/Detector.h"
#include "simX/Detection/W1.h"
#include "simX/Detection/SegmentedDetector.h"
#include "simX/Layer.h"

#include <ausa/util/memory>
#include <simX/daq/SimpleDAQ.h>
#include <iostream>

using namespace simX;
using namespace std;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Calibration;
using namespace simX::propagator;
using namespace simX::daq;

SUITE(SimpleDAQTest) {

    bool trigger(int d, int c, const DAQ::BufferItem& b) {
        return b.energy > 0;
    }

    DetectorCalibration buildDummyCalibration(size_t n, double slope, double offset) {
        return DetectorCalibration{vector<LinearCalibration>{n, LinearCalibration{offset, slope}}};
    }

    class SetupFixture {
        public:
            SetupFixture() 
                  : detectors({std::make_shared<W1>("W1", 16, 16,
                                                    TVector3{0, 0, 1}, TVector3{0, 0, 1}, TVector3{0, 1, 0},
                                                    false, 60E-3, 0.1E-3, 0.4E-3, 0.5E-3, 0.2E-3,
                                                    3.120, 3.120, 0.1, 0.1, 90E-3
            )}),
                    daq(detectors, trigger, trigger)
            {}
        std::vector<std::shared_ptr<Detector>> detectors;
        SimpleDAQ daq;
    };


    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(SetupFixture, WithNoInputOutputIsEmpty) {
        DAQ::Output o;
        daq.getData(o);
        CHECK_EQUAL(0, o.size());
    }

    TEST_FIXTURE(SetupFixture, WithSingleInputOutputIsEqualToThat) {
        Detector::DetectorOutput detectorOutput = {{4, 600, 800}}; // channel 4, energy 600, time 800
        daq.feed(0, detectorOutput);
        DAQ::Output o;
        daq.getData(o);
        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(0, o[0].detector);
        CHECK_EQUAL(4, o[0].channel);
        CHECK_EQUAL(600, o[0].energy);
    }

    TEST_FIXTURE(SetupFixture, TwoDetectionsWithoutClearWillJustSum) {
        Detector::DetectorOutput detectorOutput = {{4, 600, 800}};
        daq.feed(0, detectorOutput);
        daq.feed(0, detectorOutput);
        DAQ::Output o;
        daq.getData(o);
        CHECK_EQUAL(1, o.size());
        CHECK_CLOSE(1200, o[0].energy, 1E-6);
    }

    TEST_FIXTURE(SetupFixture, AfterClearingDataIsGone) {
        Detector::DetectorOutput detectorOutput = {{4, 600, 800}};
        daq.feed(0, detectorOutput);
        daq.clear();
        DAQ::Output o;
        daq.getData(o);
        CHECK_EQUAL(0, o.size());
    }

    TEST_FIXTURE(SetupFixture, IfEnergyIsAboveTDCThresholdDAQTriggers) {
        DAQ::Output o;
        Detector::DetectorOutput detectorOutput = {{4, 600, 800}};
        daq.feed(0, detectorOutput);
        CHECK(daq.getData(o));
        daq.clear();

        daq.setCommonTDCThreshold(500.);
        daq.feed(0, detectorOutput);
        CHECK(daq.getData(o));  // should still trigger
        daq.clear();

        daq.setCommonTDCThreshold(700.);
        daq.feed(0, detectorOutput);
        CHECK(!daq.getData(o)); // should no longer trigger
    }

    TEST_FIXTURE(SetupFixture, IfEnergyIsAboveADCThresholdItIsReadOut) {
        DAQ::Output o;
        Detector::DetectorOutput detectorOutput = {{4, 600, 800},{5, 2000, 800}};
        daq.feed(0, detectorOutput);
        daq.getData(o);
        CHECK_EQUAL(2, o.size());
        daq.setCommonADCThreshold(500.);
        daq.feed(0, detectorOutput);
        daq.getData(o);
        CHECK_EQUAL(2, o.size());   // should still be ok
        daq.clear();

        daq.setCommonADCThreshold(700.);
        daq.feed(0, detectorOutput);
        daq.getData(o);
        CHECK_EQUAL(1, o.size()); // no longer registered
    }


    TEST_FIXTURE(SetupFixture, InverseCalibrationIsAppliedToData) {
        detectors[0]->setCalibration(buildDummyCalibration(32, 1, 100));

        Detector::DetectorOutput detectorOutput = {{4, 600, 800}}; // channel 4, energy 600, time 800
        daq.feed(0, detectorOutput);
        DAQ::Output o;
        daq.getData(o);
        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(0, o[0].detector);
        CHECK_EQUAL(4, o[0].channel);
        CHECK_EQUAL(500, o[0].energy);
    }

    TEST_FIXTURE(SetupFixture, InverseCalibrationIsAppliedToData2) {
        detectors[0]->setCalibration(buildDummyCalibration(32, 100, 0));

        Detector::DetectorOutput detectorOutput = {{4, 600, 800}}; // channel 4, energy 600, time 800
        daq.feed(0, detectorOutput);
        DAQ::Output o;
        daq.getData(o);
        CHECK_EQUAL(1, o.size());
        CHECK_EQUAL(0, o[0].detector);
        CHECK_EQUAL(4, o[0].channel);
        CHECK_EQUAL(6, o[0].energy);
    }
}
