//
// Created by kuhlwein on 1/9/20.
//

#ifndef JONAS_BALAMUTHINTERFERENCEWEIGHT_H
#define JONAS_BALAMUTHINTERFERENCEWEIGHT_H
#include <vector>
#include <memory>
#include "DecayWeight.h"
#include "Level.h"
#include "Channel.h"
#include "ClebschGordan.h"
#include "SphericalHarmonic.h"
#include "BalamuthWeight.h"

using namespace std;

class BalamuthInterferenceWeight : public DecayWeight {
private:
    unique_ptr<BalamuthWeight> balamuth1;
    unique_ptr<BalamuthWeight> balamuth2;
public:
    BalamuthInterferenceWeight(unique_ptr<BalamuthWeight>, unique_ptr<BalamuthWeight>);
    ~BalamuthInterferenceWeight();

    double Calculate(vector<TLorentzVector> &);
    std::vector<double> CalculateInterference(vector<TLorentzVector> &);
};


#endif //JONAS_BALAMUTHINTERFERENCEWEIGHT_H
