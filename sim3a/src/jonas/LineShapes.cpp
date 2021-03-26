#include "jonas/LineShapes.h"
#include <TMath.h>

using namespace TMath;

double MonoAlpha(double *var, double *par)
{
  /*
  This function describes the lineshape of a monoenergetic alpha-line in a Si-detector.
  It is a Gaussian folded with two exponential tails, see for instance 'Collaers and Bortels (1987)'.
  Ideally the parameters sigma, tau1, tau2 and eta are the same for all peaks in a spectrum (or are they...?).
  */
  double x= var[0];
  double area = par[0];  //Peak area.
  double mu = par[1];    //centroid of Gaussian.
  double sigma = par[2]; //Width of Gaussian.
  double tau1 = par[3];  //Fall-off of first tail.
  double tau2 = par[4];  //Fall-off of second tail.
  double eta = par[5];   //Weighting between the two tail components.

  return area/2. * ((1.-eta)/tau1 * Exp((x-mu)/tau1 + 0.5 * Power(sigma/tau1,2)) * Erfc(1./Sqrt(2) * ((x-mu)/sigma + sigma/tau1))
                       + eta/tau2 * Exp((x-mu)/tau2 + 0.5 * Power(sigma/tau2,2)) * Erfc(1./Sqrt(2) * ((x-mu)/sigma + sigma/tau2)));
}


double SimpleAlpha(double *var, double *par)
{
  //Also a Gaussian folded with only one tail is possible
  double x = var[0];
  double area = par[0];  //Peak area.
  double mu = par[1];    //centroid of Gaussian.
  double sigma = par[2]; //Width of Gaussian.
  double tau = par[3];  //Fall-off of tail.

  return area/tau * Exp((x-mu)/tau + 0.5 * Power(sigma/tau,2)) * Erfc(1./Sqrt(2) * ((x-mu)/sigma + sigma/tau));
}
