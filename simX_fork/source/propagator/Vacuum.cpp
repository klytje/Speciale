//
// Created by munk on 06-08-15.
//

#include "simX/propagator/Vacuum.h"
#include "simX/Layer.h"
#include "simX/Particle.h"

#include <limits>

size_t simX::propagator::propagateInVacuum(const std::vector<const Layer*>& layers, simX::Particle & part){
    auto pos = part.getPosition();
    auto& dir = part.getDirectionLab();

    TVector3 intersection;
    double thickness, distance;

    double best = std::numeric_limits<double>::infinity();
    size_t result = (size_t) -1;

    for (size_t i = 0; i < layers.size(); ++i) {
        auto& layer = layers[i];
        if (layer == nullptr) continue;

        auto hit = layer->getIntersection(pos, dir, intersection, thickness, distance);

        if (hit && distance < best) {
            best = distance;
            part.setPosition(intersection);
            result = i;
        }
    }

    return result;
}
