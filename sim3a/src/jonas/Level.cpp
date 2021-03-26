#include <vector>
#include <assert.h>
#include "jonas/Level.h"

Level::Level(double e, double j, int p)
{
  Elevel = e;
  Jlevel = j;
  SetParity(p);
}

Level::~Level() {}

void Level::SetEnergy(double e)
{
  Elevel = e;
}

void Level::SetJ(double j)
{
  Jlevel = j;
}

void Level::SetParity(int p)
{
  assert(p == -1 || p == 1);
  Plevel = p;
}

double Level::Energy() const
{
  return Elevel;
}

double Level::E() const
{
  return Elevel;
}

double Level::J() const
{
  return Jlevel;
}

int Level::Parity() const
{
  return Plevel;
}

int Level::Pi() const
{
  return Plevel;
}
