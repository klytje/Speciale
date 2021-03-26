//
// Created by kuhlwein on 1/15/20.
//

#include "jonas/TripleInterferenceDecay.h"

#include <utility>

double TripleInterferenceDecay::Generate() {
    double w = 1.;
    if(doGroundState){
        wGenerator = GenerateGroundState();
        wCalculator = 1.;
    }
    else{
        wGenerator = GenerateUniform();
        wCalculator = weight->Calculate(products);
        factors.clear();
        for (auto &weight : weights) {
            factors.push_back(weight->CalculateInterference(products));
        }
        //factors = {weight->CalculateInterference(products)};
    }

    w *= wGenerator;
    w *= wCalculator;

    rawProducts = products; //Store raw decay.
    if(doRecoil) RecoilBoost();

    return w;
}

void TripleInterferenceDecay::SetDecayWeights(vector<BalamuthInterferenceWeight*> weights) {
    this->weights = std::move(weights);
}

vector<vector<double>> TripleInterferenceDecay::GetFactors() {
    return factors;
}
