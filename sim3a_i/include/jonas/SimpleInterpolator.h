#ifndef SIMPLE_INTERPOLATOR_H
#define SIMPLE_INTERPOLATOR_H
#include "Math/Interpolator.h"
#include <TH1.h>

class SimpleInterpolator {

  private:

    ROOT::Math::Interpolator *interpolator;

  public:

    /**
    * A simple wrapper around the ROOT::Math::Interpolator. The constructor should be given the
    * path to a two-column .txt-file, which are then read in and interpolated with a cubic spline.
    */
    SimpleInterpolator(const char *);
    SimpleInterpolator(TH1 *);
    ~SimpleInterpolator();

    double Eval(double);
    double Deriv(double);
    double Integ(double, double);
};
#endif
