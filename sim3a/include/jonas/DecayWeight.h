#ifndef DECAY_WEIGHT_H
#define DECAY_WEIGHT_H
#include <vector>
#include <TLorentzVector.h>

using namespace std;

class DecayWeight {
  protected:
    double scale;

  public:
    DecayWeight();
    virtual ~DecayWeight() = default;

    void SetScale(double);
    double GetScale();
    virtual double Calculate(vector<TLorentzVector> &);
};
#endif
