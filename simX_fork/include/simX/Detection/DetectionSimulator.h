
#ifndef DETECTIONSIMULATOR_h
#define DETECTIONSIMULATOR_h 1

// simX header files
#include "simX/Layer.h"
#include "simX/PhysicsEvent.h"
#include "simX/Detection/Detector.h"

// C++ and ROOT header files
#include <vector>
#include <memory>


namespace simX {

    /**
    * Class for simulating the response of the detection system.
    * @param physEvent Final-state particles ejected from target 
    */
    class DetectionSimulator {
        public:
            struct DetectionTuple {
                size_t detectorID;
                Detector::DetectorOutput output;
            };

            using Detections = std::vector<DetectionTuple>;
            using DetectorPtr = std::shared_ptr<Detector>;

            DetectionSimulator(const std::vector<DetectorPtr> detectors);
            ~DetectionSimulator();

            const Detections& run(PhysicsEvent &physEvent);

        private:
            const std::vector<DetectorPtr> detectors;
            std::vector<size_t> indexes;
            Detections detections;

            size_t propagateToNext(Particle& part, size_t max);
    };
}

#endif	/* DETECTIONSIMULATOR_H */
