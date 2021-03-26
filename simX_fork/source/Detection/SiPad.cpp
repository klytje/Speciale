#include "simX/Detection/SiPad.h"
#include <ausa/util/memory>
#include <ausa/util/stream.h>
#include <ausa/geometry/Box.h>
#include <sstream>
#include <simX/Logger.h>
#include <simX/Detection/VolumeFactory.h>

using namespace AUSA::EnergyLoss;
using namespace std;
using namespace simX;

using AUSA::EnergyLoss::Material;
using AUSA::Geometry::Box;

struct SiPad::Args : public SegmentedDetector::Args {
    Args(const TVector3 &center, const TVector3 &normal, const TVector3 &up, double activeVolumeThickness,
              double deadlayerFrontThickness, double deadlayerBackThickness, double lengthHorizontal,
              double lengthVertical, bool doIonizing) : center(center), normal(normal), up(up),
                                                        lengthHorizontal(lengthHorizontal),
                                                        lengthVertical(lengthVertical), ionizing(doIonizing)
    {
        thickness = {deadlayerFrontThickness, activeVolumeThickness, deadlayerBackThickness};
        totalThickness = 0;
        for (auto& t : thickness) totalThickness += t;
    }

    vector<SegmentedDetector::LayerPtr> buildLayers() override {
        auto Si = Material::predefined("Silicon");
        vector<Material> material{Si, Si, Si};

        auto factory = Detection::boxFactory(lengthHorizontal, lengthVertical, normal, up);
        return Detection::stackLayers(normal, thickness, material, factory, 1);
    }

    unique_ptr<Volume> buildDetectionVolume() override {
        // NOTE: detVol gets placed such that the center of its surface is at (0,0,0)
        auto pos = -0.5*totalThickness*normal;
        return make_unique<AUSA::Geometry::Box>(lengthHorizontal, lengthVertical, totalThickness, pos, normal, up);
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

    vector<SegmentedDetector::SegmentGroup> buildSegmentation() override {
        return std::vector<SegmentGroup>({
                                                 SegmentGroup{lengthHorizontal, -lengthHorizontal/2, 1, false, frontCoordinateFunction, modulate, gap}
                                         });
    }

    bool doSharing() override {
        return false;
    }

    bool doIonizing() override {
        return ionizing;
    }

    UInt_t nChannels() override {
        return 1;
    }


    static SegmentedDetector::SegmentGroup::CoordinateFunction::result_type frontCoordinateFunction(const TVector2& uv) {
        return uv.X();
    }

    static int modulate(int i) {
        return i;
    }

    static double gap(const TVector2&) {
        return 0;
    }

    TVector3 center; TVector3 normal; TVector3 up;
    double lengthHorizontal; double lengthVertical;
    bool ionizing;

    vector<double> thickness;
    double totalThickness;
};

// Constructor
SiPad::SiPad( std::string name, TVector3 center, TVector3 normal, TVector3 up,
            double activeVolumeThickness, double deadlayerFrontThickness, double deadlayerBackThickness,
            double lengthHorizontal, double lengthVertical, bool doIonizing)
    : SegmentedDetector(name,
            std::make_unique<SiPad::Args>(
                    center, normal, up,
                    activeVolumeThickness, deadlayerFrontThickness, deadlayerBackThickness,
                    lengthHorizontal, lengthVertical,
                    doIonizing
            ))
{
}

string SiPad::description() {
    stringstream out;
    string prefix = " - ";
    out << "Name: " << getName() << endl;
    out << prefix << "Type: " << "SiPad" << endl;
    out << prefix << "Position: " << getPosition() << endl;
    out << prefix << "Direction: " << getDirection() << endl;
    return out.str();
}

bool SiPad::strikesGrid(double u) {
    return false;
}

int SiPad::gridIndex() {
    return -1;
}
