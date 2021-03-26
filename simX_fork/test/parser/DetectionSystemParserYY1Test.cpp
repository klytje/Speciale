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

SUITE(DetectionSystemParserYY1Test) {
    struct Fixture {
        DetectionSystemParser parser;
        DetectionSystem system;

        Fixture()
            : parser(false)
        {
            system = parser.parseString(R"( {"ausalib_setup" : "test/_res/setupYY1/setup.json"})").system;
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
        CHECK_EQUAL(12, system.size());
    }

}
