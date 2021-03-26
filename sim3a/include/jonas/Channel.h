#ifndef CHANNEL_H
#define CHANNEL_H
#include "Math/Interpolator.h"
#include "Pair.h"
#include <memory>

class Channel {
  /**
  * This class will calculate some properties often used in R-matrix analysis.
  * It needs to know the constituent particles of the breakup, the reduced width,
  * relative angular momentum and the desired channel radius.
  * Units: Mass and energy in keV, lengths in fm.
  */

  private:
    Pair pair;       //Which particle pair does the breakup proceed through.
    int l;           //Orbital angular momentum between the fragments.
    double r;        //Channel radius in fm.
    double g;

    static constexpr double hbarc = 197326.9788;       //Planck const x speed of light [keV * fm]
    static constexpr double alpha = 7.2973525664e-3;   //Fine structure constant ~1/137.

    //The following is for fast evaluation.
    bool interpolate;  //Do interpolation or not?
    double emax, emin; //Upper and lower energy limit of the existing interpolation. emin = -0.5*emax.
    double estep;      //Grid of 'exact' points on which to perform the interpolation.

    std::unique_ptr<ROOT::Math::Interpolator> Pl;      //Penetrability.
    std::unique_ptr<ROOT::Math::Interpolator> S;       //Shift function.
    std::unique_ptr<ROOT::Math::Interpolator> dS;      //Derivative of shift function.
    std::unique_ptr<ROOT::Math::Interpolator> Phi;     //Hard-sphere phase shift.

    void MakeInterpolation(double, double);  //Initialises the interpolation.
    void CheckAndExpand(double);

    double ExactPenetrability(double);    //Evaluation using high-precision Coulomb functions.
    double FastPenetrability(double);     //Evaluation using interpolation.

    double ExactShiftFunction(double);
    double FastShiftFunction(double);

    double ExactShiftDeriv(double);
    double FastShiftDeriv(double);

    double ExactHardSphere(double);
    double FastHardSphere(double);

  public:

    /**
     *
     * @param pc
     * @param lc
     * @param gc Reduced width of channel
     * @param r0 Reduced channel radius
     * @param interpc
     */
    Channel(Pair pc = Pair(), int lc = 0, double gc = 1., double r0 = 1.42, bool interpc = false);
    ~Channel();

    void SetPair(Pair&);

    /**
    * Change the relative angular momentum between the breakup-fragments.
    */
    void SetL(int);
 
    /**
    * The channel radius may be varied. The default value is r0 = 1.42fm.
    */
    void SetRadius(double);

    const Pair & GetPair();
    int L();
    double Radius();

    /**
    * The R-matrix functions may be evaluated faster using interpolation. In this
    * implementation we use cubic spline interpolation with user defined limits
    * and step size. If the maximum energy is later exceeded, the interpolation
    * is automatically expanded to incorporate a larger energy region. One should
    * be careful with the fast option, since it takes time to initially set up the
    * interpolation.
    */
    void UseInterpolation(bool interp = true, double Emax = 10000., double Estep = 100.);

    /**
    * The following are some common functions needed in R-matrix analysis. 
    * Arguments are all energy above the channel threshold in keV.
    */
    double Penetrability(double);
    double ShiftFunction(double);
    double ShiftFunctionDeriv(double);
    double CoulombShift(double);
    double HardSphereShift(double);
    double BarrierHeight();        //Height of Coulomb barrier in [keV].
    double Rho(double);            //'Distance' parameter, (wave number) * r. 


    //TODO:Following functions should be moved elsewhere.
    double PartialWidth(double);
    double ReducedWidth();
    void SetReducedWidth(double);
};
#endif
