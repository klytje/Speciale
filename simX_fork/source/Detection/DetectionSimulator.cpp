
#include "simX/Detection/DetectionSimulator.h"
#include "simX/propagator/Vacuum.h"
#include "simX/Detection/SegmentedDetector.h"
#include "simX/Logger.h"

#include "ausa/util/stream.h"


#include <TVector3.h>
#include <numeric>

using namespace std;
using namespace simX;

DetectionSimulator::DetectionSimulator(const vector<DetectionSimulator::DetectorPtr> detectors)
 :  detectors(detectors), indexes(detectors.size()) {

    // Detector volumes
    std::iota(begin(indexes), end(indexes), 0);
}


DetectionSimulator::~DetectionSimulator() 
{}


const DetectionSimulator::Detections& DetectionSimulator::run(PhysicsEvent& physicsEvent) {
    static auto log = log::getLogger("DetectionSimulator");
    detections.clear();

    // Loop over final-state particles
    for (auto& part : physicsEvent) {

        size_t max = detectors.size();
//        log->debug("START sim of {}", part->getName());

        // Propagate particle while it still have energy.
        double E;
        while ((E = part->getKineticEnergyLab())>0.) {
//            log->debug("{} have {} keV", part->getName(), E);

            // Propagate to next detector
            auto i = propagateToNext(*part, max);


            // If particle does not hit any of the detectors, skip to next particle
            if (i==-1) break;

            // Get index of detector
            auto detectorIndex = indexes[i];
//            log->debug("{} hit detector {} = {}", part->getName(), detectorIndex, detectors[detectorIndex]->getName());

            // Propagate particle in detector and simulate response
            auto& o = detectors[detectorIndex] -> detect(*part);

            // Save result
            detections.push_back({detectorIndex, o});

            // Swap the detector we hit to the end, such that we cannot hit it again.
            indexes[i] = indexes[max-1];
            indexes[max-1] = detectorIndex;
            --max;
        }

//        log->debug("STOP sim of {}", part->getName());
    }
    return detections;
}

size_t DetectionSimulator::propagateToNext(Particle& part, size_t max) {
    auto pos = part.getPosition();
    auto& dir = part.getDirectionLab();

    TVector3 intersection;
    double thickness, distance;

    double best = std::numeric_limits<double>::infinity();
    size_t result = (size_t) -1;

    for (size_t i = 0; i < max; ++i) {
        auto& layer = detectors[indexes[i]]->getDetectorVolume();

        auto hit = layer.getIntersection(pos, dir, intersection, thickness, distance);

        if (hit && distance < best) {
            best = distance;
            part.setPosition(intersection);
            result = i;
        }
    }

    return result;
}
