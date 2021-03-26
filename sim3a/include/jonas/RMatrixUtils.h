#ifndef RMATRIX_UTILS_C
#define RMATRIX_UTILS_C
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <iostream>
#include <complex>
#include <TMath.h>

using namespace TMath;
using namespace std;

//Main functions.
double P(int, double, double);
double ShiftNeg(int, double, double);
double ShiftPos(int, double, double);

//Helper functions.
double CoulombF(int, double, double);
double CoulombG(int, double, double);
int CoulombFunctions(int, double, double, double&, double&);
int CoulombFunctions(int, double, double, double&, double&, double&, double&);
double WhittakerW(double, double, double);
struct Wparams;
double fWhittakerW(double, void *);
double WhittakerWDiff(double, double, double);
bool CheckFinity(complex<double>&, complex<double>&, complex<double>&,complex<double>&);

#endif
