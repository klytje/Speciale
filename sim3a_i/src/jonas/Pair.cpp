#include "jonas/Pair.h"
#include <TMath.h>

using namespace TMath;

Pair::Pair(Nucleus n1 , Nucleus n2)
{
  SetPair(n1,n2);
}

Pair::~Pair() {}

void Pair::SetPair(Nucleus p1, Nucleus p2)
{
  pair.clear();
  pair.push_back(p1);
  pair.push_back(p2);
}

const vector<Nucleus> & Pair::Particles() const
{
  return pair;
}

double Pair::RedMass() const
{
  double m1 = pair.at(0).M();
  double m2 = pair.at(1).M();
  double mu = m1 * m2 / (m1 + m2);

  return mu;
}

double Pair::QProduct() const
{
  double q1 = pair.at(0).Q();
  double q2 = pair.at(1).Q();

  return q1 * q2;
}

double Pair::Mtot() const
{
  return pair.at(0).M() + pair.at(1).M();
}

int Pair::Qtot() const
{
  return pair.at(0).Q() + pair.at(1).Q();
}

double Pair::WaveNumber(double E)
{
  return Sqrt(2 * Abs(E) * RedMass()) / hbarc;
}

double Pair::Eta(double E)
{
  double k = WaveNumber(E);
  return RedMass() / (k * hbarc) * alpha * QProduct();
}

double Pair::GamowFactor(double E)
{
  double eta = Eta(E);
  return Exp(-2.*Pi()*eta);
}

double Pair::Radius() const
{
  int A1 = pair.at(0).A();
  int A2 = pair.at(1).A();
  return Power(A1,1./3.) + Power(A2,1./3.);
}
