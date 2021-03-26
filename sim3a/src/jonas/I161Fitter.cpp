#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include <functional>
#include <TMath.h>
#include "jonas/I161Fitter.h"

using namespace TMath;
using namespace std;
using namespace arma;


I161Fitter::I161Fitter()
{ 
  double alphaMass = 3727379.24;
  function<double(Mat<double>)> fPx = [=](Mat<double> x) -> double {
    double result = 0;
    for(int i=0; i<3; i++){
      double dx = x[4*i+1] - x[12];
      double dy = x[4*i+2] - x[13];
      double dz = x[4*i+3] - x[14];
      double d = Sqrt(Power(dx,2) + Power(dy,2) + Power(dz,2));
      result += Sqrt(2*alphaMass*x[4*i]) * dx / d;
    }
    return result;
  };
  function<double(Mat<double>)> fPy = [=](Mat<double> x) -> double {
    double result = 0;
    for(int i=0; i<3; i++){
      double dx = x[4*i+1] - x[12];
      double dy = x[4*i+2] - x[13];
      double dz = x[4*i+3] - x[14];
      double d = Sqrt(Power(dx,2) + Power(dy,2) + Power(dz,2));
      result += Sqrt(2*alphaMass*x[4*i]) * dy / d;
    }
    return result;
  };
  function<double(Mat<double>)> fPz = [=](Mat<double> x) -> double {
    double result = 0;
    for(int i=0; i<3; i++){
      double dx = x[4*i+1] - x[12];
      double dy = x[4*i+2] - x[13];
      double dz = x[4*i+3] - x[14];
      double d = Sqrt(Power(dx,2) + Power(dy,2) + Power(dz,2));
      result += Sqrt(2*alphaMass*x[4*i]) * dz / d;
    }
    return result;
  };


  AddConstraint(fPx);
  AddConstraint(fPy);
  AddConstraint(fPz);
}

I161Fitter::~I161Fitter(){ }
