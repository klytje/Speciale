//
// Created by munk on 24-09-15.
//

#include <simX/SimulationRunner.h>
#include <simX/git.h>
#include <iostream>
#include <cmath>
#include <unistd.h>
#include <getopt.h>
#include <simX/parser/ConfigParser.h>
#include <ausa/util/FileUtil.h>
#include <ausa/parser/UnitParser.h>
#include <simX/Logger.h>
#include <simX/Random.h>

using namespace std;
using namespace simX;
using namespace AUSA::Parser;

// Functions for parsing commandline options
void parseCommandline(int argc, char* argv[]);
void printUsage(char* argv[]);
void printLongUsage(char* argv[]);


// Settings
std::string outputFile, configFile, targetFile, detectionFile, reactionFile, beamFile;
double energy = -1;
size_t N = 100000;
int progress = true;
int inverseCalibrate = true;

int main(int argc, char* argv[]) {

    parseCommandline(argc, argv);

//

    SimulationRunner runner(outputFile, configFile, targetFile, detectionFile, reactionFile, beamFile, energy, inverseCalibrate);
    runner.run(N, progress != 0);

}

void parseCommandline(int argc, char* argv[]) {
    auto logger = simX::log::getLogger("Main");

    if (argc==1) printUsage(argv);
    int c;
    string optName;
    while (1)
    {
        struct option long_options[] =
                {
                        {"calibrate",           no_argument, &inverseCalibrate, 0},
                        {"help",                no_argument, 0, 'h'},
                        {"debug",               optional_argument, 0, 0},
                        {"detection-system",    required_argument, 0, 'd'},
                        {"target",              required_argument, 0, 't'},
                        {"output",              required_argument, 0, 'o'},
                        {"seed",                required_argument, 0, 's'},
                        {"reaction",            required_argument, 0, 'r'},
                        {"beam",                required_argument, 0, 'b'},
                        {"config",              required_argument, 0, 'c'},
                        {"energy",              required_argument, 0, 'e'},
                        {"version",             no_argument, 0, 'v'},
                        {"progress",            no_argument, 0, 'P'},
                        {"quiet",               no_argument, 0, 'Q'},
                        {0,                     required_argument, 0, 'N'},
                        {0, 0, 0, 0}
                };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "d:t:o:s:r:c:e:b:hN:PQv",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        unsigned hash = 5381;

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;

                optName = long_options[option_index].name;
                if (optName == "debug") {
                    spdlog::level::level_enum l = spdlog::level::debug;
                    if(optarg != nullptr) {
                        l = (spdlog::level::level_enum) stoi(optarg);
                    }
                    simX::log::setLogLevel(l);
                    logger->set_level(l);
                    logger->info("Loglevel {} => {}", l, spdlog::level::level_names[l]);
                }
                break;

            case 'd':
                detectionFile = optarg;
                break;

            case 't':
                targetFile = optarg;
                break;

            case 'o':
                outputFile = optarg;
                break;

            case 's':
                for (size_t i = 0; i < strlen(optarg); ++i)
                    hash = 33 * hash + (unsigned char)optarg[i];
                simX::setSeed(hash);
                break;

            case 'r':
                reactionFile = optarg;
                break;

            case 'c':
                configFile = optarg;
                break;

            case 'b':
                beamFile = optarg;
                break;

            case 'e':
                energy = UnitParser('k', "eV").parse(optarg);
                break;

            case 'N':
                N = stoull(optarg);
                break;

            case 'P':
                progress = true;
                break;

            case 'v':
                cout << "simX " << GIT_HASH << endl;
                exit(0);

            case 'Q':
                progress = false;
                break;

            case 'h':
                printLongUsage(argv);
                exit(0);

            case '?':
                /* getopt_long already printed an error message. */
                printUsage(argv);
                break;

            default:
                abort ();
        }
    }

    if (configFile.empty())
        throw std::runtime_error("You must provide a config file!");

    if (outputFile.empty())
        throw std::runtime_error("You must provide a output file!");

    auto config = parser::ConfigParser().parseFile(configFile);
    auto prefix = AUSA::extractDirectory(configFile) + "/";

    if (detectionFile.empty()) {
        if (config.hasDetectionSystemFile())
            detectionFile = prefix + config.getDetectionSystemFile();
        else
            throw std::runtime_error("No detection system file provided!");
    }

    if (reactionFile.empty()) {
        if (config.hasReactionFile())
            reactionFile = prefix + config.getReactionFile();
        else
            throw std::runtime_error("No reaction file provided!");
    }

    if (targetFile.empty()) {
        if (config.hasTargetFile())
            targetFile = prefix + config.getTargetFile();
        else
            throw std::runtime_error("No target file provided!");
    }

    if (beamFile.empty()) {
        if (config.hasBeamFile())
            beamFile = prefix + config.getBeamFile();
        else
            throw std::runtime_error("No beam file provided!");
    }


}

void printUsage(char* argv[]) {
    cerr << "Usage " << argv[0] << " [options] -c config_file -o output.root" << endl;
    exit(1);
}

void printLongUsage(char* argv[]) {
    cerr << "Usage " << argv[0] << " [options] -c config_file -o output.root" << endl;
    cerr << "Options: " << endl;

    cerr << "  -c/config                       Config file                    (mandatory)" << endl;
    cerr << "  -o/output                       Output file                    (mandatory)" << endl;
    cerr << "  -s/seed                         Seed for TRandom               (mandatory)" << endl;
    cerr << "  -d/detection_system             File with detection system     (overrides config)" << endl;
    cerr << "  -t/target                       Target json file               (overrides config)" << endl;
    cerr << "  -r/reaction                     Reaction file                  (overrides config)" << endl;
    cerr << "  -b/beam                         Beam file                      (overrides config)" << endl;
    cerr << "  -e/energy                       Energy of beam particle        (overrides config)" << endl;
    cerr << "  -N                              Number of events to simulate   (default = 100000)" << endl;
    cerr << "  -P/progress                     Show progress                  (default)" << endl;
    cerr << "  --calibrate                     Write energy to the h101 tree  " << endl;
    cerr << "  --debug                         Show debug messages            (default = no debug info)" << endl;
    cerr << "  --v/version                     Print version" << endl;
    cerr << "  -Q/quiet                        Do not show progress           " << endl;
}
