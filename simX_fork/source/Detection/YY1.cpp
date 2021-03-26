//
// Created by munk on 27-03-17.
//

#include "simX/Detection/YY1.h"
#include "simX/Detection/VolumeFactory.h"
#include <vector>
#include <ausa/geometry/CylinderSegment.h>
#include <ausa/util/memory>
#include <iostream>
#include <boost/math/constants/constants.hpp>

using namespace AUSA::EnergyLoss;
using namespace simX;
using namespace simX::Detection;
using namespace std;

namespace {

    /**
     * This detector consist of 13 full size strips of with an opening angle of roughly 40 degrees.
     * These are follow by 3 strips of decreasing size.
     *
     * We model this as a wedge with the full opening angle.
     * If the particle hits the outer strips, then we calculate whether it is inside or outside the strips.
     * If outside we return an invalid coordinate.
     */

    constexpr auto PI = boost::math::constants::pi<double>();

    const auto N_SEGMENTS = 16;
    const auto N_SMALL_SEGMENTS = 3;
    const auto PITCH = 5.;
    const auto INNER_RADIUS = 50.;
    const auto OUTER_RADIUS = INNER_RADIUS+N_SEGMENTS*PITCH;
    const auto MEAN_RADIUS = (OUTER_RADIUS + INNER_RADIUS)/2;
    const auto FULL_OPENING = 42.;

    const auto STEP_IN = 0.3+0.7; // Distance between edge of silicon and active area
    const auto STEP_IN_14 = 5.266;
    const auto STEP_IN_15 = 7.961;
    const auto STEP_IN_16 = 11.807;

    double calcOpening() {
        auto R = INNER_RADIUS + (N_SEGMENTS-N_SMALL_SEGMENTS)*PITCH/2;
        auto reduction = (2*STEP_IN) / (2*PI*R) * 360;
        return FULL_OPENING - reduction;
    }

    array<double,3> calcReducedOpening() {
        auto step = STEP_IN;
        double steps[] = {STEP_IN_14, STEP_IN_15, STEP_IN_16};

        auto fullOpeningRad = FULL_OPENING*PI/180.;

        array<double, 3> result;

        for (int i = 0; i < 3; ++i) {
            step += steps[i];
            auto strip = 14+i;
            auto R = INNER_RADIUS+(strip-0.5)*PITCH;

            result[i] = (fullOpeningRad - (2*step)/R)/2.;
        }
        return result;
    };

    const auto OPENING = calcOpening();
    const auto HALF_REDUCED_OPENINGS = calcReducedOpening();

    double coordinateFunction(const TVector2& uv) {
        /*
         * uv is position relative to center of YY1.
         * y is along orientation vector
         * x is along cross vector
         */

        const auto x = uv.X();
        const auto y = MEAN_RADIUS + uv.Y();
        auto r = sqrt(x*x + y*y);

        if (r == OUTER_RADIUS) r -= 1E-9;

        if (r <= OUTER_RADIUS - 3*PITCH) return r;

        // Do a small displacement so boundary belongs to previous strip
        const auto deg = abs(std::atan2(y, -x) - TMath::PiOver2());

        if (r < OUTER_RADIUS - 2*PITCH) {
            if (deg > HALF_REDUCED_OPENINGS[0]) return -1;
        }
        else if (r < OUTER_RADIUS - 1*PITCH) {
            if (deg > HALF_REDUCED_OPENINGS[1]) return -1;
        }
        else if (r < OUTER_RADIUS - 0*PITCH) {
            if (deg > HALF_REDUCED_OPENINGS[2]) return -1;
        }

        return r;
    }

    double modulationFunction(int i) {
        return i;
    }

    double gapFunction(const TVector2& uv) {
        return 0;
    }
}

struct YY1::Args : public SegmentedDetector::Args {

    Args(const TVector3 &position, const TVector3 &normal, const TVector3 &orientation, double activeThickness,
         double dlFront, double dlBack) : position(position), normal(normal), orientation(orientation)
    {
        thickness = {dlFront, activeThickness, dlBack};
        totalThickness = 0;
        for (auto& t : thickness) totalThickness += t;

        outerRad = INNER_RADIUS + N_SEGMENTS*PITCH;
    }

    vector<LayerPtr> buildLayers() override {
        using namespace simX::Detection;

        auto Si = Material::predefined("Silicon");

        vector<Material> material{Si, Si, Si};

        auto factory = hollowCylinderSegmentFactory(INNER_RADIUS, outerRad, OPENING, normal, orientation);

        return stackLayers(normal, thickness, material, factory, 1);
    }

    unique_ptr<Volume> buildDetectionVolume() override {
        using namespace simX::Detection;

        auto Si = Material::predefined("Silicon");

        TVector3 pos = -0.5*totalThickness*normal;

        return std::make_unique<AUSA::Geometry::CylinderSegment>(pos, normal, orientation,
                                                                 INNER_RADIUS, outerRad, OPENING, totalThickness);
    }

    TVector3 buildCenter() override {
        return position;
    }

    TVector3 buildNormal() override {
        return normal;
    }

    TVector3 buildOrientation() override {
        return orientation;
    }

    vector<SegmentGroup> buildSegmentation() override {
        return std::vector<SegmentGroup>(
            {SegmentGroup{PITCH, INNER_RADIUS, N_SEGMENTS, false, coordinateFunction,modulationFunction,gapFunction}}
        );
    }

    bool doSharing() override {
        return false;
    }

    bool doIonizing() override {
        return false;
    }

    UInt_t nChannels() override {
        return N_SEGMENTS;
    }

    const TVector3 &position;
    const TVector3 &normal;
    const TVector3 &orientation;

    vector<double> thickness;
    double totalThickness, outerRad;
};


YY1::YY1(const std::string &name, const TVector3 &position, const TVector3 &normal, const TVector3 &orientation,
         double activeThickness, double dlFront, double dlBack)
    : SegmentedDetector(name, std::make_unique<Args>(position, normal, orientation, activeThickness, dlFront, dlBack))
{


}

bool YY1::strikesGrid(double) {
    return false;
}

int YY1::gridIndex() {
    return -1;
}
