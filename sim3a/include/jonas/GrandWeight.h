#ifndef GRAND_WEIGHT_H
#define GRAND_WEIGHT_H
#include <vector>
#include <armadillo>
#include <memory>
#include <complex>
#include "jonas/DecayWeight.h"
#include "jonas/DecayScheme.h"
#include "jonas/Level.h"
#include "jonas/Channel.h"
#include "jonas/PrimaryLevel.h"
#include "jonas/SecondaryLevel.h"

using namespace std;
using namespace arma;

class GrandWeight : public DecayWeight {

  private:
    unique_ptr<DecayScheme> scheme;

    double primThreshold;  //8Be + alpha threshold
    double secThreshold;   //alpha + alpha threshold
    double tripleThreshold;  //3alpha threshold

    /**
    * We calculate the Clebsch-Gordan coefficients from the Wigner 3j-symbols instead
    * of using the otherwise excellent class provided in the present library.
    */
    double ClebschGordan(double,double,double,double,double,double);

    Mat<complex<double>> ConstructPrimaryMatrix(int, double);
    Mat<complex<double>> ConstructSecondaryMatrix(int, double);

    //Variables related to the final state Coulomb interaction.
    bool doCorrection;
    double correctionRadius;
    vector<unique_ptr<Channel>> c1, c2;

    //A couple of simple utility functions. Should probably live somewhere else.
    int delta(int, int); //Kronecker delta
    bool IsAlmostZero(double);
    vector<int> GetPrimaryIndices(int);  //Indices of primary levels with specified spin.
    vector<int> GetSecondaryIndices(int);  //Do. for secondary levels.
    double CoulombCorrection(double, double, double, double, int, int);
    double CorrectionRadius(double, double, int, int);//(primChannels.at(k),secChannels.at(k),secLevels.at(k),E,E23);

  public:
    GrandWeight(unique_ptr<DecayScheme>);
    ~GrandWeight();

    /**
    * Activates correction for Coulomb interaction in the final state for near-direct decays.
    */
    void DoCorrection(double rt = -1.);

    double Calculate(vector<TLorentzVector> &);
};
#endif
