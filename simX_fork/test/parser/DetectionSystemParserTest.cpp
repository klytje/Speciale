//
// Created by munk on 28-06-15.
//


#include <unittest++/UnitTest++.h>

#include "simX/parser/DetectionSystemParser.h"

#include <typeinfo>
#include <simX/Detection/S3.h>
#include <simX/Detection/W1.h>
#include <simX/Detection/SiPad.h>
#include <simX/Detection/Scaler.h>

using namespace simX::parser;
using namespace simX::detection;
using namespace simX;
using std::string;
using BufferItem = simX::daq::DAQ::BufferItem;

SUITE(DetectionSystemParserTest) {
    struct Fixture {
        DetectionSystemParser parser;
        DetectionSystem system;

        Fixture()
            : parser(false)
        {
            system = parser.parseString(R"( {"ausalib_setup" : "test/_res/setup/setup.json"})").system;
        }
    };

    std::shared_ptr<Detector> findDet(const std::vector<std::shared_ptr<Detector>>& v, const std::string& s) {
        auto i = std::find_if(begin(v), end(v), [&](const std::shared_ptr<Detector> d) {
            return d->getName() == s;
        });
        return *i;
    }

    TEST_FIXTURE(Fixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(Fixture, DetectionSystemHas6Detectors) {
        CHECK_EQUAL(6, system.size());
    }


    TEST_FIXTURE(Fixture, FirstDetectorIsS3) {
        auto d = system.getDetectors()[0];

        CHECK(dynamic_cast<S3*>(d.get()) != nullptr);
    }

    TEST_FIXTURE(Fixture, ThirdDetectorIsW1) {
        auto d = system.getDetectors()[2];

        CHECK(dynamic_cast<W1*>(d.get()) != nullptr);
    }

    TEST_FIXTURE(Fixture, FifthDetectorIsPad) {
        auto d = system.getDetectors()[4];

        CHECK(dynamic_cast<SiPad*>(d.get()) != nullptr);
    }

    TEST_FIXTURE(Fixture, NumberOfStripsIsPropagatedToW1) {
        auto& dets = system.getDetectors();
        auto d = findDet(dets, "Det2");
        auto& w = dynamic_cast<W1&>(*d);

        CHECK_EQUAL(24, w.frontStripCount());
        CHECK_EQUAL(32, w.backStripCount());
    }

    TEST_FIXTURE(Fixture, ParserAcceptsSimpleOneDetectorTrigger) {
        parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU"
})");
    }

    TEST_FIXTURE(Fixture, S3CanHaveSharing) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "SU" : {
          "sharing" : true
         }
}
})");

        auto& dets = res.system.getDetectors();
        auto d = findDet(dets, "SU");
        auto& w = dynamic_cast<S3&>(*d);

        CHECK(w.isSharingEnabled());
    }

    TEST_FIXTURE(Fixture, W1CanHaveSharing) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "Det1" : {
          "sharing" : true
         }
}
})");

        auto& dets = res.system.getDetectors();
        auto d = findDet(dets, "Det1");
        auto& w = dynamic_cast<W1&>(*d);

        CHECK(w.isSharingEnabled());
    }

    TEST_FIXTURE(Fixture, __Global__SharingIsPropagated) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "__GLOBAL__" : {
          "sharing" : true
         }
}
})");

        auto& dets = res.system.getDetectors();
        auto d = findDet(dets, "Det1");
        auto& w = dynamic_cast<W1&>(*d);

        CHECK(w.isSharingEnabled());
    }

    TEST_FIXTURE(Fixture, LocalSettingOverrideGlobals) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "__GLOBAL__" : {
          "sharing" : true
         },
  "Det1" : {
          "sharing" : false
         }
}
})");

        auto& dets = res.system.getDetectors();
        auto d = findDet(dets, "Det1");
        auto& w = dynamic_cast<W1&>(*d);

        CHECK(!w.isSharingEnabled());
    }

    TEST_FIXTURE(Fixture, LocalAndGlobalSettingsAreMerged) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "__GLOBAL__" : {
          "sharing" : true
         },
  "Det1" : {
         }
}
})");

        auto& dets = res.system.getDetectors();
        auto d = findDet(dets, "Det1");
        auto& w = dynamic_cast<W1&>(*d);

        CHECK(w.isSharingEnabled());
    }

    TEST_FIXTURE(Fixture, SUHaveAdcThreshold100keV) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "SU" : {
          "trigger_threshold" : "100keV"
         }
}
})");

        auto& dets = res.system;

        auto& f = dets.getTDCTrigger();

        BufferItem b;

        CHECK(f(0,0, (b.energy = 200, b)));
        CHECK(!f(0,0, (b.energy = 99, b)));


    }


    TEST_FIXTURE(Fixture, DetectionSystemHave4Scalers) {
        CHECK_EQUAL(4, system.getScalers().size());
    }

    TEST_FIXTURE(Fixture, FirstScalerNamed) {
        CHECK_EQUAL("CLOCK", system.getScalers()[0]->getName());
    }

    TEST_FIXTURE(Fixture, SuCanHaveResolutionOf100keV) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "SU" : {
          "trigger_threshold" : "100keV",
          "resolution" : {"type": "GAUSS",
                          "fwhm":
                             {  "front":"100keV",
                                "back":"50keV"
                             }
                         }
         }
}
})");
    }

    TEST_FIXTURE(Fixture, W1CanDetectOnlyIonizingEnergy) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "Det1" : {
          "ionizing" : true
         }
}
})");

        auto& dets = res.system.getDetectors();
        auto d = findDet(dets, "Det1");
        auto& w = dynamic_cast<W1&>(*d);

        CHECK(w.isDetectOnlyIonizingEnergyEnabled());
    }

    TEST_FIXTURE(Fixture, S3CanDetectOnlyIonizingEnergy) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "SU" : {
          "ionizing" : true
         }
}
})");

        auto& dets = res.system.getDetectors();
        auto d = findDet(dets, "SU");
        auto& s3 = dynamic_cast<S3&>(*d);

        CHECK(s3.isDetectOnlyIonizingEnergyEnabled());
    }

    TEST_FIXTURE(Fixture, IonizingCanBeSetGlobally) {
        auto res = parser.parseString(R"( {
"ausalib_setup" : "test/_res/setup/setup.json",
"daq_trigger" : "SU & Det2",
"detectors" : {
  "__GLOBAL__" : {
          "ionizing" : true
         },
  "Det1" : {
         }
}
})");

        auto& dets = res.system.getDetectors();
        auto dw = findDet(dets, "Det1");
        auto& w = dynamic_cast<W1&>(*dw);
        auto ds3 = findDet(dets, "SU");
        auto& s3 = dynamic_cast<S3&>(*ds3);

        CHECK(w.isDetectOnlyIonizingEnergyEnabled());
        CHECK(s3.isDetectOnlyIonizingEnergyEnabled());
    }
}
