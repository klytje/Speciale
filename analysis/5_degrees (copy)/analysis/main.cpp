
#include "ausa/event/builder/DefaultEventBuilder.h"
#include "ausa/event/EventRunner.h"
#include "ausa/event/ParticleConfiguration.h"
#include "ausa/event/id_cut/TdcCutter.h"
#include "ausa/event/cut/TotalEnergyCutter.h"

// user includes
#include "IFA006Analyzer.h"
#include "MulSpectrum.h"

#include "ausa/identify/IdReader.h"
#include "AntiProtonCutter.h"
#include <ausa/sort/SortedReader.h>
#include <ausa/identify/IdReader.h>
#include <ausa/json/IO.h>
#include <ausa/eloss/Ion.h>
#include <ausa/util/FileUtil.h>
#include <ausa/TFileWrapper.h>


#include <iostream>
#include <utility>   // std::pair

using namespace std;
using namespace AUSA;
using namespace AUSA::Identify;
using namespace AUSA::EnergyLoss;
using namespace AUSA::JSON;
using namespace AUSA::Event;


int main(int argc, char* argv[]) {

//    // Beam energy, setup and target
    double beamEnergy = 669; // keV
    string setupFile  = "setup/setup.json";
    string targetFile = "setup/target.json";

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

    // TDC cutter
    string name     = "Tdc";
    double cutMin   = -15E3; // ps
    double cutMax   =  15E3; // ps
    double cutThres = 0;   // keV
    auto tdcCutter = make_shared<TdcCutter> (name, cutMin, cutMax, cutThres);

    auto mulPlot = std::make_shared<MulSpectrum>();

    auto antiProtonCutter = make_shared<AntiProtonCutter>(400);

    // Attach cutters
	//builder -> attachCutter(tdcCutter);
    builder -> attachCutter(mulPlot);
//    builder -> attachCutter(antiProtonCutter);

    // Create runner
    EventRunner runner(setup, move(builder));

    // Attach analyzer
    auto output = argv[1];
    TFile f{output, "RECREATE"};
    auto a = std::make_shared<IFA006Analyzer>();
    runner.setVerbose(true);
    runner.attach(a);

    // Add sorted and identified files
    vector<string> files;
    for (int i = 2; i < argc; ++i) findFilesMatchingWildcard(argv[i], files);
    runner.addFiles(files);

    // Run analysis
    runner.run();

    // Save to ROOT file
//    auto output = argv[1];
    f.WriteTObject(&tdcCutter->getHistAll());
    f.WriteTObject(&tdcCutter->getHistSurvivor());
    f.WriteTObject(&mulPlot->mul);

    f.Close();

    runner.saveToRootFile(output, "UPDATE");

    cout << antiProtonCutter->rejected / ((double) antiProtonCutter->accepted) << endl;

    return 0;
}   

