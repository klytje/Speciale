#include <iostream>
#include "jonas/DecayWeight.h"

DecayWeight::DecayWeight()
{
  scale = 1.0;
}

void DecayWeight::SetScale(double input)
{
  scale = input;
}

double DecayWeight::GetScale()
{
  return scale;
}
 
double DecayWeight::Calculate(vector<TLorentzVector> &p)
{
  return 1.0;
}
