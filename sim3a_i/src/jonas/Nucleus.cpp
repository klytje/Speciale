#include "jonas/Nucleus.h"
#include <cmath>
#include <assert.h>
#include <ausa/constants/Mass.h>

Nucleus::Nucleus(int A, int Z, double J, int parity)
{
  a = A;
  q = Z;
  m = AUSA::Constants::isotopeMass(Z,A);
  Level l(0.0,J,parity);
  levels.push_back(l);
}

Nucleus::~Nucleus(){}

int Nucleus::AddLevel(Level& l)
{
  int id = levels.size();
  levels.push_back(l);
  return id;
}

const vector<Level> & Nucleus::Levels()
{
  return levels;
}

void Nucleus::SetMass(double mass)
{
  m = mass;
}

void Nucleus::SetCharge(int Q)
{
  q = Q;
}

void Nucleus::SetJ(double J)
{
  levels.at(0).SetJ(J);
}

void Nucleus::SetParity(int parity)
{
  assert(parity == -1 || parity == 1);
  levels.at(0).SetParity(parity);
}

double Nucleus::Mass() const
{
  return m;
}

double Nucleus::M() const
{
  return m;
}

int Nucleus::Charge() const
{
  return q;
}

int Nucleus::Q() const
{
  return q;
}

double Nucleus::J() const
{
  return levels.at(0).J();
}

int Nucleus::Parity() const
{
  return levels.at(0).Pi();
}

int Nucleus::Pi() const
{
  return levels.at(0).Pi();
}

int Nucleus::A() const
{
  return a;
}

int Nucleus::Z() const
{
  return q;
}

int Nucleus::IsotopeID()
{
  return lrint(0.5 * (a + q) * (a + q + 1) + q);
}
