//
// Created by kuhlwein on 1/9/20.
//

#include "jonas/BalamuthInterferenceWeight.h"

BalamuthInterferenceWeight::~BalamuthInterferenceWeight() {}

BalamuthInterferenceWeight::BalamuthInterferenceWeight(unique_ptr<BalamuthWeight> b1, unique_ptr<BalamuthWeight> b2) {
    balamuth1 = std::move(b1);
    balamuth2 = std::move(b2);
}

double BalamuthInterferenceWeight::Calculate(vector<TLorentzVector> &p) {
    auto a = balamuth1->CalculateAmplitudes(p);
    auto b = balamuth2->CalculateAmplitudes(p);
    double sum = 0;

    return sum;
}

std::vector<double> BalamuthInterferenceWeight::CalculateInterference(vector<TLorentzVector> &p) {
    auto a = balamuth1->CalculateAmplitudes(p);
    auto b = balamuth2->CalculateAmplitudes(p);
    double f1=0, f2=0, re=0, im=0;
    for (int i=0; i<a.size(); i++) {
        f1 += norm(a[i]);
        f2 += norm(b[i]);
        re += real(a[i]*conj(b[i]));
        im += imag(a[i]*conj(b[i]));
    }

    return {f1,f2,re,im};
}