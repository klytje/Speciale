//
// Created by munk on 10-11-17.
//

#ifndef JONAS_SPHERICALHARMONIC_H
#define JONAS_SPHERICALHARMONIC_H

#include <complex>
#include <memory>
#include <vector>
#include <Math/Interpolator.h>

class SphericalHarmonic {

    public:
    explicit SphericalHarmonic(unsigned int L);
    bool isInterpolated() const;
    void interpolate(bool value, size_t nPoints);

    std::complex<double> eval(int m, double theta, double phi);

    private:
    unsigned int L;
    bool interpolated;
    double legendre(unsigned int m, double theta);

    std::vector<std::unique_ptr<ROOT::Math::Interpolator>> interpolations;
};

inline std::complex<double> spherical_harmonic(unsigned int L, int m, double theta, double phi) {
    return SphericalHarmonic{L}.eval(m, theta, phi);
}


#endif //JONAS_SPHERICALHARMONIC_H
