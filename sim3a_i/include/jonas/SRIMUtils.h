#ifndef SRIMUTILS_H
#define SRIMUTILS_H
#include <string>
#include <TH1I.h>
#include <TF1.h>
#include <memory>

using namespace std;

/**
* Given the name and path to a standard SRIM range file this function returns
* a ROOT histogram with the ranges from the file.
*/
TH1I * readRangeFile(string);

/**
* Given the name and path to a standard SRIM range file this function returns
* a pointer to a ROOT TF1 describing the implantation depth profile. This 
* function can then be sampled in a Monte Carlo simulation.
*/
unique_ptr<TF1> getRangeProfile(string);
#endif
