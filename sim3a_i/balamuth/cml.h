//
// Created by munk on 08-11-17.
//

#ifndef JONAS_CML_H
#define JONAS_CML_H

#include <string>

struct CommandLineOptions {
    double ExC;
    double ExBe;
    size_t N;
    int L, J1, J2, L2;
    int p1, p2;
    double r0p, r0s;
    double gs;
    bool coloumbCorrection, interference;
    double phase, ratio;
    double coloumbRadius;
    std::string output;
};

static const char *const J1 = "JC";
static const char *const J2 = "JBe";
static const char *const P1 = "pC";
static const char *const P2 = "pBe";
static const char *const R0C = "r0C";
static const char *const R0Be = "r0Be";
static const char *const gBe = "gBe";
static const char *const ExC = "ExC";
static const char *const ExBe = "ExBe";
static const char *const SEED = "seed";
static const char *const COLOUMB_R = "coulomb-radius";
static const char *const OUTPUT = "output";
static const char *const L2 = "L2";
static const char *const ratio = "ratio";
static const char *const phase = "phase";
static const char *const HELP = "help";

CommandLineOptions parseCml(int argc, char* argv[]);


#endif //JONAS_CML_H
