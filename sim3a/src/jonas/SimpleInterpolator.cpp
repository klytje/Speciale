#include "jonas/SimpleInterpolator.h"
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>


using namespace std;

SimpleInterpolator::SimpleInterpolator(const char *dataFile)
{
  ifstream data(dataFile);

  vector<double> xv, yv;

  string line;
  while(getline(data,line)){
    double x, y;
    stringstream ss(line);
    ss >> x >> y;
    xv.push_back(x);
    yv.push_back(y);
  }

  int N = xv.size();
  interpolator = new ROOT::Math::Interpolator(N,ROOT::Math::Interpolation::kCSPLINE);
  interpolator->SetData(xv, yv);
}

SimpleInterpolator::SimpleInterpolator(TH1 *data)
{
  vector<double> xv, yv;

  int n = data->GetXaxis()->GetNbins();
  for(int i=1; i<n; i++){
    double x, y;
    x = data->GetBinCenter(i);
    y = data->GetBinContent(i);
    xv.push_back(x);
    yv.push_back(y);
  }

  int N = xv.size();
  interpolator = new ROOT::Math::Interpolator(N,ROOT::Math::Interpolation::kCSPLINE);
  interpolator->SetData(xv, yv);
}

SimpleInterpolator::~SimpleInterpolator()
{
  delete interpolator;
}

double SimpleInterpolator::Eval(double arg)
{
  return interpolator->Eval(arg);
}

double SimpleInterpolator::Deriv(double arg)
{
  return interpolator->Deriv(arg);
}

double SimpleInterpolator::Integ(double arg1, double arg2)
{
  return interpolator->Integ(arg1,arg2);
}
