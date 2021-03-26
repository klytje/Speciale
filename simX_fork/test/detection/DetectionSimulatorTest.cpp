
#include <unittest++/UnitTest++.h>

#include <simX/propagator/IonizingPropagator.h>
#include <simX/Particle.h>
#include <simX/PhysicsEvent.h>
#include <simX/Detection/DetectionSimulator.h>

#include <ausa/eloss/SRIMTabulation.h>
#include <ausa/util/Resource.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>
#include <ausa/util/memory>
#include <ausa/eloss/RangeInterpolator.h>
#include <ausa/eloss/Ion.h>
#include <simX/Detection/SegmentedDetector.h>
#include <simX/Detection/W1.h>

using namespace simX;
using namespace std;
using namespace AUSA::EnergyLoss;
using namespace simX::propagator;

SUITE(DetectionSystemTest) {

    // use SRIM stopping powers for alpha in Silicon
    IonizingPropagator::LossCalc factory(const Layer& layer, const Particle& p) {
        std::string suffix = "/SRIM08/He4_Si.dat";
        SRIMTabulation tab(AUSA::getResourceDirectory() + suffix);
        return std::make_unique<EnergyLossRangeInverter>(std::make_unique<RangeInterpolator>(tab));
    }

    class SetupFixture {
        public:
            SetupFixture()
                    : alpha(Ion(2,4)),
                      particles({&alpha}),
                      evt(particles),
                      sim({std::make_shared<W1>("W1", 16,16,TVector3(0,0,30), TVector3(0,0,-1), TVector3(0,1,0), false, 60E-3, 0.1E-3, 0.4E-3, 0.5E-3, 0.2E-3,
                                                3.120, 3.120, 0.1, 0.1, 90E-3)}),
                      ionProp( std::make_shared<IonizingPropagator>(IonizingPropagator(factory)) )
            {
                alpha.setFourMomentumLab( 1E3, {0,0,1} );
                alpha.setPropagator( ionProp );
            }

        Particle alpha;
        std::vector<Particle*> particles;
        PhysicsEvent evt;
        DetectionSimulator sim;
        std::shared_ptr<IonizingPropagator> ionProp;
    };


    TEST_FIXTURE(SetupFixture, Sanity) {
        CHECK_EQUAL(1, 1);
    }

    TEST_FIXTURE(SetupFixture, SimCanRun) {
        sim.run(evt);
    }

    TEST_FIXTURE(SetupFixture, AlphaIsDetectedInFrontAndBackOfDetector) {
        auto& out = sim.run(evt);

        // We detected a single particle.
        CHECK_EQUAL(1, out.size());

        auto& d = out[0];
        auto& detection = d.output;

        CHECK_EQUAL(2, detection.size());

        auto& d0 = detection[0];
        auto& d1 = detection[1];
        CHECK_CLOSE(969., d0.energy, 1);
        CHECK_CLOSE(969., d1.energy, 1);
    }
}
