//
// Created by munk on 26-03-17.
//

#include <unittest++/UnitTest++.h>

#include <simX/Detection/VolumeFactory.h>
#include <simX/Detection/VolumeFactory.h>
#include <TVector3.h>
#include <ausa/geometry/Box.h>
#include <ausa/util/memory>
#include <ausa/geometry/HollowCylinder.h>

using namespace simX;
using namespace simX::Detection;
using namespace AUSA::Geometry;

SUITE(VolumeFactoryTest) {
    struct Fixture {

        Fixture()
                : Si(AUSA::EnergyLoss::Material::predefined("Silicon")),
                  boxFac(boxFactory(100, 100, TVector3(0, 0, 1), TVector3(0, 1, 0))),
                  cylFac(hollowCylinderFactory(10, 100, TVector3(0, 0, 1)))
        {

        }

        VolumeFactory boxFac, cylFac;
        AUSA::EnergyLoss::Material Si;
    };

    TEST_FIXTURE(Fixture,Sanity) {
        CHECK_EQUAL(1,1);
    }

    TEST_FIXTURE(Fixture,EmptyVectorsGivesEmptyResult) {
        auto res = stackLayers(TVector3(0,0,1), {}, {}, boxFac, 0);

        CHECK_EQUAL(0, res.size());
    }

    TEST_FIXTURE(Fixture,StarkIsOneLayer) {
        auto res = stackLayers(TVector3(0,0,1), {1}, {Si}, boxFac, 0);

        CHECK_EQUAL(1, res.size());
    }

    TEST_FIXTURE(Fixture,FirstMaterialIsSilicon) {
        auto res = stackLayers(TVector3(0,0,1), {1}, {Si}, boxFac, 0);

        CHECK_EQUAL(Si.getName(), res[0]->getMaterial().getName());
    }

    TEST_FIXTURE(Fixture,FirstLayerIsBox) {
        auto res = stackLayers(TVector3(0,0,1), {1}, {Si}, boxFac, 0);

        auto& first = res.front();
        auto& vol = first->getVolume();

        CHECK(dynamic_cast<const Box*>(&vol) != nullptr);
    }

    TEST_FIXTURE(Fixture,FirstLayerIsActive) {
        auto res = stackLayers(TVector3(0,0,1), {1}, {Si}, boxFac, 0);

        auto& first = res.front();
        CHECK(first->isActive());
    }

    TEST_FIXTURE(Fixture,SecondLayerIsActive) {
        auto res = stackLayers(TVector3(0,0,1), {1, 1}, {Si, Si}, boxFac, 1);

        CHECK(res[1]->isActive());
    }

    TEST_FIXTURE(Fixture,BoxIs1mmThick) {
        auto res = stackLayers(TVector3(0,0,1), {1}, {Si}, boxFac, 0);

        auto& first = res.front();
        auto& vol = first->getVolume();

        auto box = dynamic_cast<const Box*>(&vol);
        CHECK_CLOSE(1, box->getNDim(), 1E-5);
    }

    TEST_FIXTURE(Fixture,FirstBoxIsPlacedAtHalfThickness) {
        auto res = stackLayers(TVector3(0,0,1), {1}, {Si}, boxFac, 0);

        auto& first = res.front();
        auto& vol = first->getVolume();

        CHECK_ARRAY_CLOSE(TVector3(0, 0, -0.5), vol.getCenter(), 3, 1E-5);
    }

    TEST_FIXTURE(Fixture,SecondBoxIsPlacedAfterFirst) {
        auto res = stackLayers(TVector3(0,0,1), {1, 1}, {Si, Si}, boxFac, 0);

        auto& first = res[1];
        auto& vol = first->getVolume();

        CHECK_ARRAY_CLOSE(TVector3(0, 0, -1.5), vol.getCenter(), 3, 1E-5);
    }

    TEST_FIXTURE(Fixture,LayersAreShiftedOppositeNormal) {
        auto res = stackLayers(TVector3(0,0,-1), {1}, {Si}, boxFac, 0);

        auto& first = res.front();
        auto& vol = first->getVolume();

        CHECK_ARRAY_CLOSE(TVector3(0, 0, 0.5), vol.getCenter(), 3, 1E-5);
    }

    TEST_FIXTURE(Fixture,LayersWith0ThicknessIsSkipped) {
        auto res = stackLayers(TVector3(0,0,-1), {0, 1}, {Si, Si}, boxFac, 1);

        CHECK_EQUAL(1, res.size());
    }

    TEST_FIXTURE(Fixture,IfSomeLayerIs0ThickActiveIsStillCorrect) {
        auto res = stackLayers(TVector3(0,0,-1), {0, 1}, {Si, Si}, boxFac, 1);

        CHECK(res[0]->isActive());
    }

    TEST_FIXTURE(Fixture,PossibleToMakeHollowCylinder) {
        auto res = stackLayers(TVector3(0,0,-1), {0, 1}, {Si, Si}, cylFac, 1);

        CHECK(dynamic_cast<const HollowCylinder*>(&res.front()->getVolume()));
    }
}