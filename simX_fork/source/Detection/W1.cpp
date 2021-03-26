#include "simX/Detection/W1.h"
#include "simX/Detection/VolumeFactory.h"

#include <ausa/util/memory>
#include <ausa/util/stream.h>
#include <ostream>
#include <sstream>
#include <ausa/geometry/Box.h>

using namespace std;
using namespace simX;

using AUSA::EnergyLoss::Material;

struct W1::Args : public SegmentedDetector::Args {

    Args(UInt_t nFrontStrip, UInt_t nBackStrip, const TVector3 &center, const TVector3 &normal, const TVector3 &up,
             bool sharing, double activeVolumeThickness, double deadlayerFrontThickness, double deadlayerBackThickness,
             double gridThickness, double contactBackThick, double stripFrontWidth, double stripBackWidth,
             double stripFrontGap, double stripBackGap, bool reverseOrderingFront, bool reverseOrderingBack)
            : nFrontStrip(nFrontStrip), nBackStrip(nBackStrip), center(center), normal(normal),
              up(up), sharing(sharing),
              stripFrontWidth(stripFrontWidth),
              stripBackWidth(stripBackWidth), stripFrontGap(stripFrontGap),
              stripBackGap(stripBackGap),
              reverseOrderingFront(reverseOrderingFront),
              reverseOrderingBack(reverseOrderingBack) {

        thickness = {gridThickness, deadlayerFrontThickness, activeVolumeThickness, deadlayerBackThickness, contactBackThick};
        totalThickness = 0;
        for (auto& t : thickness) totalThickness += t;

        frontWidth = nFrontStrip*stripFrontWidth;
        backWidth = nBackStrip*stripBackWidth;
    }

    vector<LayerPtr> buildLayers() override {
        auto Al = Material::predefined("Aluminum");
        auto Si = Material::predefined("Silicon");

        vector<Material> material{Al, Si, Si, Si, Al};
        size_t activeIndex = 2;

        auto factory = Detection::boxFactory(frontWidth, backWidth, normal, up);

        return Detection::stackLayers(normal, thickness, material, factory, activeIndex);
    }

    unique_ptr<Volume> buildDetectionVolume() override {
        // NOTE: detVol gets placed such that the center of its surface is at (0,0,0)
        auto pos = -0.5*totalThickness*normal;
        return make_unique<AUSA::Geometry::Box>(frontWidth, backWidth, totalThickness, pos, normal, up);
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

    bool doSharing() override {
        return sharing;
    }

    bool doIonizing() override {
        return false;
    }

    UInt_t nChannels() override {
        return nFrontStrip+nBackStrip;
    }

    vector<SegmentGroup> buildSegmentation() override {
        return std::vector<SegmentGroup>({
                                          SegmentGroup{stripFrontWidth, -frontWidth/2, nFrontStrip, reverseOrderingFront, frontCoordinateFunction, modulate, gapFunction(stripFrontGap)},
                                          SegmentGroup{stripBackWidth, -backWidth/2, nBackStrip, !reverseOrderingBack, backCoordinateFunction, modulate, gapFunction(stripBackGap)},
                                  });
    }

    static SegmentedDetector::SegmentGroup::CoordinateFunction::result_type frontCoordinateFunction(const TVector2& uv) {
        return uv.X();
    }

    static SegmentedDetector::SegmentGroup::CoordinateFunction::result_type backCoordinateFunction(const TVector2& uv) {
        return uv.Y();
    }

    static int modulate(int i) {
        return i;
    }

    SegmentedDetector::SegmentGroup::GapFunction gapFunction(double gap) {
        return [=](const TVector2&) {
            return gap;
        };
    }


    UInt_t nFrontStrip; UInt_t nBackStrip;
    TVector3 center; TVector3 normal; TVector3 up; bool sharing;
    double stripFrontWidth; double stripBackWidth; double stripFrontGap; double stripBackGap;
    bool reverseOrderingFront; bool reverseOrderingBack;

    vector<double> thickness;
    double totalThickness;

    double frontWidth, backWidth;
};


// Constructor
W1::W1(std::string name, UInt_t nFrontStrip, UInt_t nBackStrip,
       TVector3 center, TVector3 normal, TVector3 up, bool sharing, double activeVolumeThickness,
       double deadlayerFrontThickness, double deadlayerBackThickness, double gridThickness, double contactBackThick,
       double stripFrontWidth, double stripBackWidth, double stripFrontGap, double stripBackGap, double gridWidth,
        bool reverseOrderingFront, bool reverseOrderingBack )
    :
        SegmentedDetector(name, std::make_unique<W1::Args>(nFrontStrip, nBackStrip,
                                                           center, normal, up,
                                                           sharing,
                                                           activeVolumeThickness, deadlayerFrontThickness, deadlayerBackThickness, gridThickness, contactBackThick,
                                                           stripFrontWidth, stripBackWidth, stripFrontGap, stripBackGap,
                                                           reverseOrderingFront, reverseOrderingBack )),

      sideLength(nFrontStrip*stripFrontWidth), gridWidth(gridWidth),
      stripWidthFront(stripFrontWidth), stripWidthBack(stripBackWidth),
      stripGapFront(stripFrontGap), stripGapBack(stripBackGap),
      nBackStrip(nBackStrip), nFrontStrip(nFrontStrip)

{
    gridI = gridThickness != 0 ? 0 : -1;
}


// To speed up the code, we assume one electrode on each strip
// instead of three. Should't affect the physics.
bool W1::strikesGrid(double u) {
    double u0 = -0.5*sideLength;
    double w = stripWidthFront;
    double halfg = 0.5*gridWidth;

    // Which strip has been hit?
    int i = (u-u0)/w;

    // Are we hitting the electrode grid?
    double du = u - (u0+i*w) - 0.5*w;    
    return abs(du)<halfg;
}

int W1::gridIndex() {
    return gridI;
}

string W1::description() {
    stringstream out;
    string prefix = " - ";
    out << "Name: " << getName() << endl;
    out << prefix << "Type: " << "W1" << endl;
    out << prefix << "Strips: " << frontStripCount() << "x" << backStripCount() << endl;
    out << prefix << "Position: " << getPosition() << endl;
    out << prefix << "Direction: " << getDirection() << endl;
    out << prefix << "Up: " << getUp() << endl;
    out << prefix << "Front strip width/gap: " << stripWidthFront << "/" << stripGapFront << endl;
    out << prefix << "Back strip width/gap: " << stripWidthBack << "/" << stripGapBack << endl;
    return out.str();
}

void W1::reverseStripOrdering(size_t i, bool b) {
    if (i == 1) b = !b;

    SegmentedDetector::reverseStripOrdering(i, b);
}

UInt_t W1::frontStripCount() const {
    return nSegments(0);
}

UInt_t W1::backStripCount() const {
    return nSegments(1);
}
