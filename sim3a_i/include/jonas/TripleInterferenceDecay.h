//
// Created by kuhlwein on 1/15/20.
//

#ifndef JONAS_TRIPLEINTERFERENCEDECAY_H
#define JONAS_TRIPLEINTERFERENCEDECAY_H


#include <jonas/TripleDecay.h>
#include <jonas/BalamuthInterferenceWeight.h>

class TripleInterferenceDecay : public TripleDecay {
private:
    vector<BalamuthInterferenceWeight*> weights;
    vector<vector<double>> factors;
public:
    double Generate();
    void SetDecayWeights(vector<BalamuthInterferenceWeight*>);
    vector<vector<double>> GetFactors();
};


#endif //JONAS_TRIPLEINTERFERENCEDECAY_H
