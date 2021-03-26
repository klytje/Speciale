
#include <unittest++/UnitTest++.h>

#include <ausa/eloss/Ion.h>
#include <ausa/constants/Mass.h>
#include <simX/propagator/IonizingPropagator.h>

#include "simX/Particle.h"
#include "simX/Layer.h"

#include <ausa/eloss/SRIMTabulation.h>
#include <ausa/eloss/RangeInterpolator.h>
#include <ausa/eloss/EnergyLossRangeInverter.h>
#include <ausa/util/Resource.h>
#include <ausa/util/memory>
#include <ausa/geometry/Box.h>

using namespace simX;
using namespace simX::propagator;
using namespace AUSA::EnergyLoss;

using AUSA::Geometry::Box;

SUITE(IonizingTest) {

    IonizingPropagator::LossCalc factory(const Layer& layer, const Particle& p) {
        std::string suffix = "/SRIM08/";
        if (layer.getMaterial().getName() == "Silicon")
            suffix += "He4_Si.dat";
        else
            suffix += "He4_Fe.dat";

        SRIMTabulation tab(AUSA::getResourceDirectory() + suffix);

        return std::make_unique<EnergyLossRangeInverter>(std::make_unique<RangeInterpolator>(tab));
    }

    class SetupFixture {
    public:
        SetupFixture()
                : silicon("Silicon", 2.3290, 14, AUSA::Constants::atomicMass(14)),
                  iron("Iron", 7.874, 26, AUSA::Constants::atomicMass(26)),
                  p(Ion(2,4)),
                  b0(Material::predefined("Silicon"),std::make_unique<Box>(100.,100.,100.E-6,TVector3(0.,0,0.),TVector3(0,0,-1),TVector3(0,1,0)), false),
                  b1(Material::predefined("Carbon"),std::make_unique<Box>(100.,100.,100.E-6,TVector3(0.,0,0.),TVector3(0,0,-1),TVector3(0,1,0)), false),
                  prop(factory)
        {
            p.setPosition({0,0,-50});
            p.setFourMomentumLab(1E3, {0,0,1});
        }

        Material silicon, iron;
        Particle p;
        Layer b0, b1;
        IonizingPropagator prop;
    };

    TEST_FIXTURE(SetupFixture, Sanity) {
        factory(b0, p); // Check this does not crash

        CHECK_EQUAL(1,1);
    }

    TEST_FIXTURE(SetupFixture, ThicknessAreCorrectlyConvertedFrom_mm_to_nm) {
        auto c = factory(b0, p);
        auto dist = 100E-6; // 100nm
        auto loss = c -> getTotalEnergyLoss(p.getKineticEnergyLab(), dist);

        auto rest = p.getKineticEnergyLab() - loss;

        prop.propagate(b0, p, dist);

        auto E = p.getKineticEnergyLab();

        CHECK_CLOSE(rest, E, 1E-3);
    }

    TEST_FIXTURE(SetupFixture, ThicknessIsNotHardcoded) {
        auto c = factory(b0, p);
        auto dist = 200E-6; // 100nm
        auto loss = c -> getTotalEnergyLoss(p.getKineticEnergyLab(), dist);

        auto rest = p.getKineticEnergyLab() - loss;

        prop.propagate(b0, p, dist);

        auto E = p.getKineticEnergyLab();

        CHECK_CLOSE(rest, E, 1E-3);
    }

    TEST_FIXTURE(SetupFixture, ParticlePositionIsUpdated) {
        p.setFourMomentumLab(1E5, {0,0,1});
        prop.propagate(b0, p, 2);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(-50+2, pos.Z(), 1E-5);
    }

    TEST_FIXTURE(SetupFixture, IfNoRangeSpecifiedPropagateThrough) {
        p.setFourMomentumLab(1E5, {0,0,1});

        Layer b(silicon, std::make_unique<Box>(100,100,2,TVector3(0,0.,0.),TVector3(0,0,-1),TVector3(0,1,0)), false);

        prop.propagate(b, p);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(-50+2, pos.Z(), 1E-5);
    }

    TEST_FIXTURE(SetupFixture, IonizingPropagatorUsesDifferentCalcForDifferentLayers) {
        auto cSi = factory(b0, p);
        auto cFe = factory(b1, p);

        auto E = 1E3;

        auto eSi = cSi -> getTotalEnergyLoss(E, 100E-6);
        auto eFe = cFe -> getTotalEnergyLoss(E, 100E-6);

        p.setFourMomentumLab(E, {0,0,1});

        prop.propagate(b0, p, 100E-6);

        CHECK_CLOSE(E-eSi, p.getKineticEnergyLab(), 1E-3);

        p.setFourMomentumLab(E, {0,0,1});
        prop.propagate(b1, p, 100E-6);

        CHECK_CLOSE(E-eFe, p.getKineticEnergyLab(), 1E-3);

        CHECK(eSi != eFe);
    }

    TEST_FIXTURE(SetupFixture, IfParticleIsStoppedPositionIsSetAccordingly) {
        p.setFourMomentumLab(500, {0,0,1});
        auto c = factory(b0, p);
        double R;
        c -> getTotalEnergyLoss(p.getKineticEnergyLab(), 3000E-6, R);

        prop.propagate(b0, p, 3000E-6);

        auto& pos = p.getPosition();
        CHECK_CLOSE(0, pos.X(), 1E-5);
        CHECK_CLOSE(0, pos.Y(), 1E-5);
        CHECK_CLOSE(-50+R, pos.Z(), 1E-3);
    }
}
