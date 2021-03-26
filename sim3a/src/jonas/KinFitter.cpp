#include <iostream>
#include <TMath.h>
#include "jonas/KinFitter.h"

using namespace TMath;

KinFitter::KinFitter()
{
  initialised = false;
  h = 0.1;
  nMax = 50;
  eps = 0.01;
}

KinFitter::~KinFitter(){ }

void KinFitter::AddConstraint(function<double(Mat<double>)> c)
{
  /*
  * This function is used to add a constraint to the fitter. The fitter will try
  * to find a solution where all the constraints evaluate to zero.
  */
  constraints.push_back(c);
  int K = constraints.size();
  g.zeros(K,1);
}

void KinFitter::SetDiffStep(double size)
{
  /*
  * The fitter uses a crude, homemade, method for differentiation of the multidimensional
  * constraint functions. Basically
  *
  *  (dg/dxi)(x) = (g(...,xi+h/2,...) - g(...,xi-h/2,...))/h
  *
  * If the constraints are well behaved, a small step size should not be a problem.
  */
  h = size;
}

void KinFitter::SetMaxIterations(int n)
{
  /*
  * How many iterations you want the minimiser to spend on a single fit.
  */
  nMax = n;
}

void KinFitter::SetEpsilon(double num)
{
  /*
  * Define the success criterion for the minimisation. If the relative improvement
  * in the constructed chi-square is smaller than eps from one iteration to the next,
  * the fitting is regarded as successful and terminates.
  */
  eps = num;
}

void KinFitter::SetMeasurement(vector<double> m, Mat<double> M)
{
  /*
  * Use this method to give a measurement and its associated covariance-matrix to the
  * fitter. The covariance matrix must be a square matrix with side-length equal to
  * the length of the measurement. Remember to call the Fit()-method before retrieving
  * the result from the fitter.
  */
  if(m.size() == M.n_rows && m.size() == M.n_cols){
    y = m;
    x = y;
    V = M;
    initialised = true;
  }
  else{
    initialised = false;
  }
}

void KinFitter::CalculateDerivatives()
{
  int K = constraints.size();
  int N = x.n_rows;
  G.zeros(K,N);
  for(int i=0; i<K; i++){
    for(int j=0; j<N; j++){
      double xmid = x[j];
      double xlo = xmid - h/2.;
      double xhi = xmid + h/2.;
      x[j] = xlo;
      double flo = constraints[i](x);
      x[j] = xhi;
      double fhi = constraints[i](x);
      x[j] = xmid;
      double d = (fhi - flo) / h;
      G(i,j) = d;
    }
  }
}

void KinFitter::EvaluateConstraints()
{

  int K = constraints.size();
  for(int i=0; i<K; i++){
    g[i] = constraints[i](x);
  }
}

bool KinFitter::Fit()
{
  /*
  * When you have added a measurement to the fitter, you are ready to call the Fit()-method,
  * which is the work-horse of the KinFitter. After successful fitting you can retrieve the
  * result with GetResult().
  */
  if(!initialised) return false;
  EvaluateConstraints();
  mat Vi;
  if(!pinv(Vi,V)) return false;
  int nIterations = 0;
  double chiSquareOld = 1000000; //Just some ridiculous value.
  double chiSquare = 0.99 * (1.-eps) * chiSquareOld;
  while(TMath::Abs(chiSquare-chiSquareOld)/chiSquareOld > eps && nIterations < nMax){
    chiSquareOld = chiSquare;
    mat xOld = x;
    CalculateDerivatives();
    mat r = g + G * (y - x);
    mat S = G * V * G.t();
    mat Si;
    if(!pinv(Si,S)){
      x = y;
      return false;
    }
    mat l = Si * r;
    x = y - V * G.t() * l;
    EvaluateConstraints();
    mat chiSquareM = ((y.t()-x.t()) * Vi * (y-x)) + (2. * l.t() * g);
    chiSquare = chiSquareM(0,0);
    nIterations++;
  }
  return true;
}

vector<double> KinFitter::GetResult()
{
  /*
  * Retrieve the result from the kinematic fitting. Should only be called after a measurement
  * has been added with AddMeasurement() and Fit() has been succesfully run.
  */
  return conv_to<vector<double>>::from(x);
}
