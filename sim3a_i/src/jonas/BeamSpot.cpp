#include "jonas/BeamSpot.h"

BeamSpot::BeamSpot(unique_ptr<TF1> d, unique_ptr<TF2> t, double thickness) : 
  depthProfile(move(d)), transverseProfile(move(t)), foilThickness(thickness * 1e7) {}

BeamSpot::~BeamSpot() {}

double BeamSpot::GetDepth()
{
  double z = -1.;
  while(z < 0 || z > foilThickness){
    z = depthProfile->GetRandom();
  }
  return z * 1e-7;  //Return answer in mm.
}

TVector3 BeamSpot::GetDecayPoint()
{
  double z = GetDepth();
  double x, y;
  transverseProfile->GetRandom2(x, y);

  return TVector3(x,y,z);
}

void BeamSpot::GetDecayPoint(TVector3 *v)
{
  double z = GetDepth();
  double x, y;
  transverseProfile->GetRandom2(x, y);

  v->SetXYZ(x,y,z);
  return;
}
