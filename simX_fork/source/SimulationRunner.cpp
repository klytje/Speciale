//
// Created by munk on 19-11-15.
//

#include "simX/SimulationRunner.h"
#include "simX/parser/DetectionSystemParser.h"
#include "simX/parser/ConfigParser.h"
#include "simX/parser/BeamParser.h"
#include "simX/parser/ReactionParser.h"
#include "simX/EventSimulator.h"
#include "simX/PhysicsContainer.h"
#include "simX/Detection/DetectionSimulator.h"
#include "simX/daq/SimpleDAQ.h"
#include "ausa/util/FileUtil.h"
#include "rapidjson/document.h"
#include "ausa/json/JSONUtil.h"
#include "ausa/json/MalformedJsonException.h"

#include <ausa/json/IO.h>
#include <simX/daq/DaqTreeWriter.h>
#include <TObjString.h>

using namespace simX;
using namespace simX::detection;
using namespace simX::daq;
using namespace AUSA;

using simX::parser::DetectionSystemParser;
using simX::parser::ConfigParser;
using simX::parser::BeamParser;
using simX::parser::ReactionParser;

namespace {
    bool trigger(int detector, int channel, const DAQ::BufferItem& b) {
        return b.energy > 0;
    }
}

struct SimulationRunner::Impl {
    Impl(std::string outputFile, std::string configFile, std::string targetFile,
            std::string detectionFile, std::string reactionFile, std::string beamFile, double energyOverride, bool inverseCalibrate)
        : output(std::make_shared<TFileWrapper>(outputFile, "RECREATE")),
          target(JSON::readTargetFromJSON(targetFile)),
          config(ConfigParser().parseFile(configFile)),
          physicsContainer(output),
          beam(BeamParser().parseFile(beamFile))

    {
        // Check beam starts outside target
        double beamX, beamY, beamZ;
        beamZ = beam.getZ0();
        beam.nominalXY(beamX, beamY);
        TVector3 beamPosition(beamX, beamY, beamZ);
        for(const Layer& l : target.getLayers()) {
            if(l.isInside(beamPosition, VOLUME_TOLERANCE)) std::cerr << "The beam starts inside the target. This can cause unwanted effects!" << std::endl;
        }

        auto dRes = DetectionSystemParser(inverseCalibrate).parseFile(detectionFile);
        detectionSystem = dRes.system;

        if (energyOverride != -1) beam.setNominalEnergy(energyOverride);


        detectionSimulator = std::make_unique<DetectionSimulator>(detectionSystem.getDetectors());

        daq = std::make_unique<SimpleDAQ>(detectionSystem.getDetectors(), dRes.system.getTDCTrigger(), dRes.system.getADCTrigger());
        daq->setTrigger(detectionSystem.getTriggerFunction());
        daqTree = std::make_unique<DaqTreeWriter>(output, detectionSystem, dRes.ausalibFile);
        chains = ReactionParser().parseMultiple(reactionFile);

        logConfig(configFile, targetFile, detectionFile, reactionFile, beamFile);

        TObjString calibrated(inverseCalibrate ? "FALSE" : "TRUE");
        output->writeWithKey(calibrated, "H101_CALIBRATED");
    }

    void setupEventSimulator() {
        // Set propagators for target/physics
        config.applyBeamPropagator(chain->getBeam());
        for (auto& r : *chain) {
            for (auto& d : r.getDaughters()) {
                config.applyTargetPropagator(*d);
            }
        }
        eventSimulator = std::make_unique<EventSimulator>(beam, target, *chain);

        // Save final state propagators for fast switching later.
        finalStates = chain->findFinalStates();

        for (auto p : finalStates) {
            targetProp.push_back(p->getPropagator());
            config.applyDetectionPropagator(*p);
            detProp.push_back(p->getPropagator());
        }

    }

    void run(int N, bool progress = false) {

        auto t = N/200;
        int j = 0;
        int total = 0;
        for(ReactionParser::ReactionBranch& b : chains) {
            chain = std::move(b.chain);
            setupEventSimulator();
            if(progress) std::cerr << b.name << std::endl;
            int iterations = N * b.ratio;
            for (int i = 0; i < iterations; i++, j++, total++) {

                runSimulationStep();

                if (progress && j == t) {
                    std::cerr << "\rProcessing entry " << total << " (" << round(100.*total/(double) N)  << "%)" << std::flush;
                    j = 0;
                }
            }
            if(progress)  {
                std::cerr << "\rProcessing entry " << total << " (" << round(100.*total/(double) N)  << "%)" << std::flush;
                std::cerr << std::endl;
            }
            b.chain = std::move(chain); // My job is done here, please put me back!
        }

    }

    void runSimulationStep() {
        // Set target propagators
        for (size_t i = 0; i < finalStates.size(); ++i) {
            finalStates[i]->setPropagator(targetProp[i]);
        }
        // Simulate physics
        auto pEvent = eventSimulator->run();
        // Fill the physics tree
        physicsContainer.fill(pEvent);

        // Set detections propagators
        for (size_t i = 0; i < finalStates.size(); ++i) {
            finalStates[i]->setPropagator(detProp[i]);
        }

        // Simulate detection
        auto& dEvent = detectionSimulator->run(pEvent);

        // Simulate daq
        daq->feed(dEvent);

        // Retrieve data.
        daq->getData(daqOutput);

        // Write to TTree
        daqTree->feed(daqOutput);

        // Clear daq so we are ready for next round
        daq->clear();
    }

    void logConfig(const std::string& configFile, const std::string& targetFile,
                 const std::string& detectionFile, const std::string& reactionFile, const std::string& beamFile) {

        logFile("CONFIG", configFile);
        logFile("BEAM", beamFile);
        logFile("TARGET", targetFile);
        logFile("DETECTION_SYSTEM", detectionFile);

        TObjArray c;
        for (auto& p : chains) {
            auto s = asciiFileToString(p.path);

            auto str = new TObjString(s.c_str());
            c.Add(str);
        }
        output->writeWithKey(c, "REACTIONS", "SIMX_CONFIG");
    }

    void logFile(const std::string& key, const std::string& path) {
        auto s = asciiFileToString(path);

        TObjString str(s.c_str());
        output->writeWithKey(str, key, "SIMX_CONFIG");
    }



    std::shared_ptr<TFileWrapper> output;
    Target target;
    DetectionSystem detectionSystem;
    Config config;
    Beam beam;
    ReactionParser::ChainPtr chain;
    std::vector<ReactionParser::ReactionBranch> chains;

    std::unique_ptr<EventSimulator> eventSimulator;
    PhysicsContainer physicsContainer;

    std::unique_ptr<DetectionSimulator> detectionSimulator;

    std::unique_ptr<SimpleDAQ> daq;
    std::unique_ptr<DaqTreeWriter> daqTree;
    DAQ::Output daqOutput;

    std::vector<Particle*> finalStates;
    std::vector<Particle::Propagator> targetProp, detProp;
};

SimulationRunner::SimulationRunner(std::string outputFile, std::string configFile, std::string targetFile,
                                   std::string detectionFile, std::string reaction, std::string beamFile, double energyOverride, bool inverseCalibrate)
    : pimpl(std::make_unique<Impl>(outputFile, configFile, targetFile, detectionFile, reaction, beamFile, energyOverride, inverseCalibrate))
{
 /* Nothing */
}

SimulationRunner::~SimulationRunner() {
    /*Here because we use pimpl*/
}

void SimulationRunner::run(int i, bool progress) {
    pimpl -> run(i, progress);
}