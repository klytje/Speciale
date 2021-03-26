//
// Created by munk on 08-11-17.
//

#include "cml.h"
#include <cmath>
#include <getopt.h>
#include <ausa/parser/UnitParser.h>
#include <string>
#include <iostream>
#include <TRandom.h>

using namespace std;

std::pair<int,int> parseSpinParity(const std::string& arg) {
    size_t pos;

    auto J = stoi(arg, &pos);
    if (pos == arg.size()) {
        cerr << "Missing parity!" << endl;
        exit(1);
    }
    else if (pos < arg.size() - 1) {
        cerr << "To much input after spin!" << endl;
        exit(1);
    }

    int parity = 0;
    auto c = arg[pos];
    switch (c) {
        case '+': parity = 1; break;
        case '-': parity = -1; break;
        default:
            cerr << "Parity should be +/-. You specified: " << arg << endl;
            exit(1);
    }
    return make_pair(J, parity);
};

void printUsage(char** argv) {
    cerr << "Usage: " << argv[0] << " options..." << endl << endl;
    
    cerr << "--" << J1 << "                 Carbon spin" << endl;
    cerr << "--" << J2 <<  "                Berylium spin" << endl;
    cerr << "--" << R0C << "                r0 for the C breakup" << endl;
    cerr << "--" << R0Be << "               r0 for the Be breakup" << endl;
    cerr << "--" << gBe << "                reduced width of the Be system" << endl;
    cerr << "--" << ExC << "                C excitation energy" << endl;
    cerr << "--" << ExBe << "               Be excitation energy" << endl;
    cerr << "--" << SEED << "               seed for random number generator" << endl;
    cerr << "--" << COLOUMB_R << "     radius used for Coloumb corrections" << endl;
    cerr << "--" << OUTPUT << "/-o          Output file (.root or text)" << endl;
    cerr << "-N                   Number of events" << endl;
    cerr << "-h                   This text" << endl;

    cerr << "Example:" << endl;
    cerr << argv[0] << "--JC 1+ --JBe 2+ -L 2 --ExC 12.7MeV "
            "--ExBe 3030keV  -N 200000 -o out.dat" << endl;
}

CommandLineOptions parseCml(int argc, char **argv) {
    if (argc == 1) {
        printUsage(argv);
        exit(1);
    }

    CommandLineOptions cml{};
    AUSA::Parser::UnitParser eParser{'k', "eV"};
    AUSA::Parser::UnitParser rParser{'f', "m"};
    AUSA::Parser::UnitParser gParser{'k', "eV^(1/2)"};

    cml.coloumbCorrection = false;
    cml.N = 100000;

    cml.J1 = cml.J2 = cml.L = -1;
    cml.p1 = cml.p2 = 0;
    cml.r0p = 1.87;
    cml.r0s = 1.42;
    cml.gs = NAN;
    cml.ExBe = cml.ExC = NAN;

    int c;
    static struct option long_options[] = {
            {J1,               required_argument, 0, 0},
            {J2,        required_argument, 0, 0},
            {R0C,       required_argument, 0, 0},
            {R0Be,      required_argument, 0, 0},
            {gBe,       required_argument, 0, 0},
            {ExC,       required_argument, 0, 0},
            {ExBe,      required_argument, 0, 0},
            {SEED,      required_argument, 0, 0},
            {COLOUMB_R, required_argument, 0, 0},
            {OUTPUT,    required_argument, 0, 'o'},
            {HELP,      required_argument, 0, 'h'},
            {0,         required_argument, 0, 'N'},
            {0,         required_argument, 0, 'L'},
            {NULL, 0,  NULL,                  0}
    };
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "N:L:o:h",
                            long_options, &option_index)) != -1) {

        switch (c) {
            case 0: {
                string opt = long_options[option_index].name;
                if (opt == J1) {
                    std::tie(cml.J1, cml.p1) = parseSpinParity(optarg);
                }
                else if (opt == J2) {
                    std::tie(cml.J2, cml.p2) = parseSpinParity(optarg);
                }
                else if (opt == R0C) {
                    cml.r0p = rParser.parse(optarg);
                }
                else if (opt == R0Be) {
                    cml.r0s = rParser.parse(optarg);
                }
                else if (opt == gBe) {
                    cml.gs = gParser.parse(optarg);
                }
                else if (opt == ExC) {
                    cml.ExC = eParser.parse(optarg);
                }
                else if (opt == ExBe) {
                    cml.ExBe = eParser.parse(optarg);
                }
                else if (opt == SEED) {
                    gRandom->SetSeed(stoul(optarg));
                }
                else if (opt == COLOUMB_R) {
                    cml.coloumbRadius = rParser.parse(optarg);
                    cml.coloumbCorrection = true;
                }
                break;
            }
                /* ----------------------------------------*/
            case 'N':
                cml.N = stoull(optarg);
                break;
            case 'L':
                cml.L = stoi(optarg);
                break;
            case 'o':
                cml.output = optarg;
                break;
            case 'h':
                printUsage(argv);
                exit(0);
            case '?':
                break;
            default:
                printUsage(argv);
        }
    }

    if (cml.J1 < 0) {
        cerr << "C spin must be > 0!" << endl;
        exit(1);
    }
    if (cml.J2 < 0) {
        cerr << "Be spin must be > 0!" << endl;
        exit(1);
    }
    if (cml.L < 0) {
        cerr << "L must be > 0!" << endl;
        exit(1);
    }
    if (cml.p1 == 0) {
        cerr << "C parity must specified!" << endl;
        exit(1);
    }
    if (cml.p2 == 0) {
        cerr << "Be parity must specified!" << endl;
        exit(1);
    }

    if (std::isnan(cml.ExC)) {
        cerr << "Please specify C excitation energy" << endl;
        exit(1);
    }

    if (cml.ExC < 0) {
        cerr << "C excitation energy must be postive!" << endl;
        exit(1);
    }

    if (cml.ExBe < 0) {
        cerr << "Be8 excitation energy must be postive!" << endl;
        exit(1);
    }

    if (std::isnan(cml.ExBe)) {
        cerr << "Please specify Be excitation energy" << endl;
        exit(1);
    }
    if (cml.output.empty()) {
        cerr << "Please specify output file" << endl;
        exit(1);
    }

    if (std::isnan(cml.gs)) {
        if (cml.ExBe < 5) {
            cml.gs = 28.81;
            cerr << "No reduced width specified and E(Be) < 5keV" << endl
                 << "Assuming ground state g = 28.81 keV^(1/2)." << endl;
        }
        else if (cml.ExBe >= 2950 && cml.ExBe < 3100) {
            cml.gs = 32.787;
            cerr << "No reduced width specified and 2950keV < E(Be) < 3100keV" << endl
                 << "Assuming first excited g = 32.787 keV^(1/2)." << endl;
        }
        else {
            cerr << "No reduced width specified and unknown Be8 state. Specify gs!" << endl;
            exit(1);
        }
    }

    auto pi = cml.p1;
    auto pf = cml.p2 * pow(-1, cml.L);
    if (pi != pf) {
        cerr << "Parity violation! " << pi << " != " << pf << "  (" << cml.p2 << " * [-1]^" << cml.L << ")" << endl;
        exit(4);
    }

    if (!(abs(cml.L - cml.J2) <= cml.J1 && abs(cml.L + cml.J2) >= cml.J1)) {
        cerr << "Spin violation!" << endl;
        exit(5);
    }

    return cml;
}