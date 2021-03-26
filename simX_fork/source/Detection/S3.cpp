#include "simX/Detection/S3.h"
#include "simX/Detection/VolumeFactory.h"

#include <ausa/util/memory>
#include <ausa/util/stream.h>
#include <sstream>
#include <ausa/geometry/HollowCylinder.h>
#include <iostream>

using namespace std;
using namespace simX;
using AUSA::EnergyLoss::Material;

struct S3::Args : public SegmentedDetector::Args {
    Args(const TVector3 &center, const TVector3 &normal, const TVector3 &up, UInt_t nSpokes, UInt_t nRings,
         double innerRadius, double ringWidth, double ringGap, double spokeGap, bool sharing, double activeThickness,
         double deadlayerFrontThickness, double deadlayerBackThickness, bool reverseSpokesOrdering,
         bool reverseRingOrdering) : center(center), normal(normal), up(up), nSpokes(nSpokes), nRings(nRings),
                                     innerRadius(innerRadius), ringWidth(ringWidth), ringGap(ringGap),
                                     spokeGap(spokeGap), sharing(sharing),
                                     reverseSpokesOrdering(reverseSpokesOrdering),
                                     reverseRingOrdering(reverseRingOrdering) {

        thickness = {deadlayerFrontThickness, activeThickness, deadlayerBackThickness};
        totalThickness = 0;
        for (auto& t : thickness) totalThickness += t;

        outerRadius = innerRadius + nRings*ringWidth;
    }

    vector<LayerPtr> buildLayers() override {
        auto Si = Material::predefined("Silicon");

        vector<Material> material{Si, Si, Si};
        auto factory = simX::Detection::hollowCylinderFactory(innerRadius, outerRadius, normal);
        return simX::Detection::stackLayers(normal, thickness, material, factory, 1);
    }

    unique_ptr<Volume> buildDetectionVolume() override {
        // NOTE: detVol gets placed such that the center of its surface is at (0,0,0)
        auto pos = -0.5*totalThickness*normal;
        return make_unique<AUSA::Geometry::HollowCylinder>(innerRadius, outerRadius, totalThickness, TVector3{0,0,0}, normal);
    }

    TVector3 buildCenter() override {
        return center;
    }

    TVector3 buildNormal() override {
        return normal;
    }

    TVector3 buildOrientation() override {
        return up;
    }

    vector<SegmentGroup> buildSegmentation() override {
    return std::vector<SegmentGroup>({
                                             SegmentGroup{2.*TMath::Pi() / nSpokes, 0, nSpokes, reverseSpokesOrdering, spokeCoordinateFunction(), modulateSpokes(nSpokes), spokeGapFunction()},
                                             SegmentGroup{ringWidth, innerRadius, nRings, reverseRingOrdering, ringCoordinateFunction, modulateRings, ringGapFunction()},
                                      });
    }

    bool doSharing() override {
        return sharing;
    }

    bool doIonizing() override {
        return false;
    }

    UInt_t nChannels() override {
        return nSpokes+nRings;
    }

    static double ringCoordinateFunction(const TVector2& uv) {
        return uv.Mod();
    }

    SegmentedDetector::SegmentGroup::CoordinateFunction spokeCoordinateFunction() {
        const auto spokeAngle = 2.*TMath::Pi() / nSpokes;
        return [=](const TVector2& uv) {
            return fmod(uv.Phi() + 3*TMath::Pi()/2. + 0.5*spokeAngle,  2*TMath::Pi());
        };
    }

    static int modulateRings(int i) {
        return i;
    }

    SegmentedDetector::SegmentGroup::ModulationFunction modulateSpokes(UInt_t nSpokes) {
        return [=](int i) {
            return i % nSpokes;
        };
    }

    SegmentedDetector::SegmentGroup::GapFunction ringGapFunction() {
        auto local = ringGap;
        return [=](const TVector2&) {
            return local;
        };
    }

    SegmentedDetector::SegmentGroup::GapFunction spokeGapFunction() {
        auto local = spokeGap;
        return [=](const TVector2& uv) {
            return local / uv.Mod();
        };
    }

    TVector3 center; TVector3 normal; TVector3 up;
    UInt_t nSpokes; UInt_t nRings;
    double innerRadius; double ringWidth; double ringGap; double spokeGap;
    bool sharing;
    bool reverseSpokesOrdering; bool reverseRingOrdering;

    vector<double> thickness;
    double totalThickness;
    double outerRadius;
};

S3::S3( std::string name, TVector3 center, TVector3 normal, TVector3 up,
        UInt_t nSpokes, UInt_t nRings, double innerRadius, double ringWidth, double ringGap, double spokeGap,
        bool sharing, double activeThickness, double deadlayerFrontThickness,
        double deadlayerBackThickness, bool reverseSpokesOrdering, bool reverseRingOrdering )
    : SegmentedDetector( name,
                         make_unique<Args>(center, normal, up,
                                           nSpokes, nRings,
                                           innerRadius, ringWidth, ringGap, spokeGap,
                                           sharing,
                                           activeThickness, deadlayerFrontThickness, deadlayerBackThickness,
                                           reverseSpokesOrdering, reverseRingOrdering
                         )
),
      nSpokes(nSpokes), nRings(nRings), innerRadius(innerRadius), ringWidth(ringWidth), ringGap(ringGap), spokeGap(spokeGap)
{

}

S3::~S3() 
{
    // Nothing
}

bool S3::strikesGrid(double u) {
    return true;
}

// S3 does not have a grid of electrodes, just a single effective deadlayer
// that covers the entire surface
int S3::gridIndex() {
    return -1;
}

string S3::description() {
    stringstream out;
    string prefix = " - ";
    out << "Name: " << getName() << endl;
    out << prefix << "Type: " << "S3" << endl;
    out << prefix << "Strips: " << spokeCount() << "x" << ringCount() << endl;
    out << prefix << "Position: " << getPosition() << endl;
    out << prefix << "Direction: " << getDirection() << endl;
    out << prefix << "Up: " << getUp() << endl;
    return out.str();
}

UInt_t S3::spokeCount() {
    return nSegments(0);
}

UInt_t S3::ringCount() {
    return nSegments(1);
}


