#ifndef SIMX_LAYER_H
#define SIMX_LAYER_H

#include <ausa/eloss/Layer.h>

namespace simX {

    /**
     * Tolerance used for determining if a particle is inside a layer
     */
    const double VOLUME_TOLERANCE = 1E-9;

    using AUSA::EnergyLoss::Layer;

}

#endif	/* SIMX_LAYER_H */
