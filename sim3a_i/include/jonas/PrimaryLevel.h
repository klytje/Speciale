#ifndef PRIMARY_LEVEL_H
#define PRIMARY_LEVEL_H
#include <armadillo>

using namespace std;
using namespace arma;

class PrimaryLevel {
  private:
    Mat<double> widths;  //Sqrt of formal reduced widths, including proper signs. nl x nmu.

  public:
    PrimaryLevel();
    ~PrimaryLevel();

    double El; //Level energy in keV above ground state.
    int J;     //Spin
    double Bg; //Beta decay feeding factor.
    
    void SetWidths(Mat<double>);
    Mat<double> & GetWidths();

    /**
    * Returns the width associated with l (index of l is first parameter) and mu (second parameter).
    */
    double GetWidth(int, int);
};
#endif
