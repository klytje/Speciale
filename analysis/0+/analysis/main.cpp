
#include "ausa/event/builder/DefaultEventBuilder.h"
#include "ausa/event/EventRunner.h"
#include "ausa/event/ParticleConfiguration.h"
#include "ausa/event/cut/TotalEnergyCutter.h"

// user includes
#include "StandardAnalyzer.h"

#include "ausa/identify/IdReader.h"
#include <ausa/sort/SortedReader.h>
#include <ausa/identify/IdReader.h>
#include <ausa/json/IO.h>
#include <ausa/eloss/Ion.h>
#include <ausa/util/FileUtil.h>
#include <ausa/TFileWrapper.h>


#include <iostream>
#include <fstream>
#include <utility>   // std::pair

using namespace std;
using namespace AUSA;
using namespace AUSA::Identify;
using namespace AUSA::EnergyLoss;
using namespace AUSA::JSON;
using namespace AUSA::Event;

// simple script which reads the energy from the beam.simX file
double readBeamEnergyFromJSON(string beamfile) {
    ifstream file(beamfile);
    string line;
    getline(file, line);

    // remove everything but the last word
    line.erase(0, line.find_last_of(' ')+1);

    // determine number of digits
    size_t i = 0;
    for (; i < line.length(); i++ ){ 
        if (!isdigit(line[i])) {
            break;
        }
    }

    // keep only the digits, removing everything else (the unit)
    line.erase(i, line.size());
    return atof(line.c_str());
}

int main(int argc, char* argv[]) {
    // Beam energy, setup and target
    string beamFile = "beam.simX";
    string setupFile  = "setup/setup.json";
    string targetFile = "setup/target.json";

    // Get the beam energy from the beam.simX file
    double beamEnergy = readBeamEnergyFromJSON("beam.simX");

    // Load AUSAlib setup file
    auto setup = readSetupFromJSON(setupFile);

    // Load Target file
    auto target = readTargetFromJSON(targetFile);

    // Create event builder
    auto builder = make_unique<DefaultEventBuilder>(target, beamEnergy, Ion("H1"), Ion("B11"));
    builder->setBailOutMultiplity(14);
    builder->setAllowBackground(true);

    // Set final state
    ParticleConfiguration fs;
    fs.push_back(make_pair(Ion("He4"), 3));
    builder -> setFinalState(fs);

    // Create runner
    EventRunner runner(setup, move(builder));

    // Attach analyzer
    auto output = argv[1];
    TFile f{output, "RECREATE"};
    auto a = std::make_shared<StandardAnalyzer>();
    runner.setVerbose(true);
    runner.attach(a);

    // Add sorted and identified files
    vector<string> files;
    for (int i = 2; i < argc; ++i) findFilesMatchingWildcard(argv[i], files);
    std::cout << "CHECKPOINT" << std::endl;
    runner.addFiles(files);
    std::cout << "CHECKPOINT" << std::endl;

    // Run analysis
    runner.run();

    // Save to ROOT file
    f.Close();
    runner.saveToRootFile(output, "UPDATE");
    return 0;
}   

