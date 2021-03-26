//
// Created by munk on 26-03-17.
//

#ifndef SIMX_VOLUMEFACTORY_H
#define SIMX_VOLUMEFACTORY_H

#include <vector>
#include <memory>
#include "simX/Layer.h"
#include "simX/Volume.h"
#include "ausa/eloss/Material.h"
#include <functional>

class TVector3;

namespace simX {
    namespace Detection {

        using VolumeFactory = std::function<std::unique_ptr<Volume>(const TVector3&, double)>;
        std::vector<std::unique_ptr<Layer>> stackLayers(const TVector3& normal, const std::vector<double>& thickness,
                                              const std::vector<AUSA::EnergyLoss::Material> material, VolumeFactory factory, size_t active);

        VolumeFactory boxFactory(double xDim, double upDim, const TVector3& normal, const TVector3& up);
        VolumeFactory hollowCylinderFactory(double innerRadius, double outerRadius, const TVector3& normal);
        VolumeFactory hollowCylinderSegmentFactory(double innerRadius, double outerRadius, double opening, const TVector3& normal, const TVector3& orientation);
    }
}
#endif //SIMX_VOLUMEFACTORY_H
