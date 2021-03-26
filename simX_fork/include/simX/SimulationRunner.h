//
// Created by munk on 19-11-15.
//

#ifndef SIMX_SIMULATIONRUNNER_H
#define SIMX_SIMULATIONRUNNER_H

#include <string>
#include <memory>

namespace simX {
    class SimulationRunner {
    public:
        SimulationRunner(std::string outputFile, std::string configFile,
                         std::string targetFile, std::string detectionFile, std::string reactionFile,
                         std::string beamFile,
                         double energyOverride = -1,
                         bool inverseCalibrate = false);
        ~SimulationRunner();

        void run(int i, bool progress = false);

    private:
        struct Impl;
        std::unique_ptr<Impl> pimpl;
    };
}

#endif //SIMX_SIMULATIONRUNNER_H
