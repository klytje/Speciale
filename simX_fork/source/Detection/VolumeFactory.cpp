//
// Created by munk on 26-03-17.
//

#include <TVector3.h>
#include <ausa/util/memory>
#include <ausa/geometry/Box.h>
#include <ausa/geometry/HollowCylinder.h>
#include <ausa/geometry/CylinderSegment.h>
#include "simX/Detection/VolumeFactory.h"

using namespace simX;
using namespace simX::Detection;
using namespace std;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Geometry;


namespace simX {
    namespace Detection {
        vector<unique_ptr<Layer>> stackLayers(const TVector3& normal,
                                              const vector<double>& thickness, const vector<Material> material, VolumeFactory factory, size_t active) {
            TVector3 pos(0,0,0);
            double last = 0;
            vector<unique_ptr<Layer>> result;
            for (size_t i = 0; i < thickness.size(); i++) {
                pos -= 0.5 * (last + thickness[i]) * normal;
                last = thickness[i];
                if (thickness[i] == 0) continue;
                result.emplace_back(make_unique<Layer>(material[i], factory(pos, thickness[i]), active==i));
            }
            return result;
        }

        VolumeFactory boxFactory(double xDim, double upDim, const TVector3 &normal, const TVector3 &up) {
            return [=](const TVector3& position, double thickness) {
                return std::make_unique<AUSA::Geometry::Box>(xDim, upDim, thickness, position, normal, up);
            };
        }

        VolumeFactory hollowCylinderFactory(double innerRadius, double outerRadius, const TVector3 &normal) {
            return [=](const TVector3& position, double thickness) {
                return std::make_unique<AUSA::Geometry::HollowCylinder>(innerRadius, outerRadius, thickness, position, normal);
            };
        }

        VolumeFactory
        hollowCylinderSegmentFactory(double innerRadius, double outerRadius, double opening, const TVector3 &normal,
                                     const TVector3 &orientation) {
            return [=](const TVector3& position, double thickness) {
                return std::make_unique<AUSA::Geometry::CylinderSegment>(position, normal, orientation,
                                                                         innerRadius, outerRadius, opening, thickness);
            };
        }
    }
}