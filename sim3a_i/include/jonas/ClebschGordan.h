#ifndef CLEBSCH_GORDAN_H
#define CLEBSCH_GORDAN_H
#include <vector>

using namespace std;

class ClebschGordan {
  /**
  * A class that can be used to calculate the Clebsch-Gordan coefficients of two
  * coupled angular momenta. The coefficients are calculated once at construction,
  * so if you need to change the angular momentum of one sub-system, you must construct
  * a new ClebschGordan object.
  *
  * Hint: It is also possible to use GSL, where the Wigner 3j-symbols are declared in
  * the header 'gsl_sf_coupling.h'. Check out the AZURE2 source code for an example.
  */

  private:
    vector<vector<vector<double>>> coefficients;
    double j1, j2;
    double eps;

    /**
    * This method makes sure that the CGformula() is called with the proper arguments.
    */
    double CalculateCG(double, double, double);

    /**
    * This is an implementation of a general formula for calculation of Clebsch-Gordan coefficients
    * found on Wikipedia. The formula only works if m1 + m2 > 0 and J1 > J2.
    * (https://en.wikipedia.org/wiki/Table_of_Clebsch%E2%80%93Gordan_coefficients).
    */
    double CGformula(double, double, double, double, double);

  public:
    /**
    * The constructor needs the angular momenta of the two sub-systems which
    * are coupled. It then calculates the Clebsch-Gordan coefficients and fills 
    * an array that can later be accessed by the user.
    */
    ClebschGordan(double J1 = 0.5, double J2 = 0.5);
    ~ClebschGordan();

    /**
    * If you want to couple some different angular momenta than was passed to the constructor.
    * A bit expensive, so shouldn't be called very often.
    */
    void Couple(double,double);

    /**
    * Call this function to extract the calculated Clebsch-Gordan coefficients. If you
    * call it with unphysical arguments it returns zero. The arguments are the magnetic
    * quantum numbers of the two sub-systems and the total angular momentum of the coupled system.
    */
    double Coefficient(double, double, double);
    double J1();
    double J2();
};
#endif
