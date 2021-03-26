//
// Created by munk on 10-11-17.
//

#include <cassert>
#include <ausa/util/memory>
#include "jonas/SphericalHarmonic.h"
#include "Math/SpecFuncMathMore.h"
#include <iostream>
#include <TMath.h>

using namespace std;

SphericalHarmonic::SphericalHarmonic(unsigned int L) : L(L), interpolated(false) {
    assert(L >= 0);
}

std::complex<double> SphericalHarmonic::eval(int m, double theta, double phi) {
    double v;
    if (interpolated) {
        v = interpolations[abs(m)]->Eval(cos(theta));
    }
    else v = legendre(abs(m), theta);

    bool r_sign = false;
    bool i_sign = false;
    if(m < 0)
    {
        // Reflect and adjust sign if m < 0:
        r_sign = m&1;
        i_sign = !(m&1);
        m = abs(m);
    }
    if(m&1)
    {
        // Check phase if theta is outside [0, PI]:
        auto mod = fmod(theta, 2*TMath::Pi());
        if(mod < 0)
            mod += 2*TMath::Pi();
        if(mod > 2*TMath::Pi())
        {
            r_sign = !r_sign;
            i_sign = !i_sign;
        }
    }

    auto r = v * cos(m * phi);
    auto i = v * sin(m * phi);
    if(r_sign)
        r = -r;
    if(i_sign)
        i = -i;

    return {r, i};
}

bool SphericalHarmonic::isInterpolated() const {
    return interpolated;
}

void SphericalHarmonic::interpolate(bool value, size_t nPoints) {
    if (!value) {
        interpolations.clear();
        return;
    }

    double diff = 2. / (nPoints - 1);
    for (unsigned int m = 0; m <= L; ++m) {
        std::vector<double> x, y;

        for (size_t i = 0; i < nPoints; ++i) {
            auto cosTheta = -1 + i * diff;
            auto theta = acos(cosTheta);
            x.push_back(cosTheta);
            y.push_back(legendre(m, theta));
        }
        interpolations.push_back(std::make_unique<ROOT::Math::Interpolator>(x, y));
        interpolations.back()->Eval(0);
    }
    interpolated = true;
}


double SphericalHarmonic::legendre(unsigned int m, double theta) {
    return ROOT::Math::sph_legendre(L, m, theta);
}
