#ifndef BALAMUTH_WEIGHT_H
#define BALAMUTH_WEIGHT_H
#include <vector>
#include <memory>
#include "DecayWeight.h"
#include "Level.h"
#include "Channel.h"
#include "ClebschGordan.h"
#include "SphericalHarmonic.h"

using namespace std;

class BalamuthWeight : public DecayWeight {

  private:
    static constexpr double hbarc = 197326.9788;       //Planck const x speed of light [keV * fm]
    Level primLevel;
    unique_ptr<Channel> primChannel;
    Level secLevel;
    unique_ptr<Channel> secChannel;
    unique_ptr<ClebschGordan> CG;
    unique_ptr<SphericalHarmonic> primSH, secSH;

    //Following variables are for the final-state Coulomb interaction.
    bool doCorrection;
    unique_ptr<Channel> c1, c2;
    double rt;
 
   /**
    * We calculate the Clebsch-Gordan coefficients from the Wigner 3j-symbols.
    */
    double GetClebschGordan(double, double, double, double, double, double);

    /**
     * Compute Clebsch Gordan
     */
    void spinCouple();

    /**
    * Phase shift associated with the resonance in the secondary (intermediate) system.
    */
    //double SecondaryPS(double);
    //double ROOTPS(double *, double *);

    /**
    * Correction radius for final state Coulomb interactions. Dependent on the relative
    * a-Be energy and the a-a energy.
    */
    //double CorrectionRadius(double, double);
    //double CorrectionRadius(unique_ptr<Channel> &, Level &, double, double);

  public:
    BalamuthWeight();
    ~BalamuthWeight();

    void SetPrimary(Level&, std::unique_ptr<Channel>);
    void SetSecondary(Level&, std::unique_ptr<Channel>);
    void DoCorrection(bool input = true, double radius = 15.);

    double Calculate(vector<TLorentzVector> &);
};
#endif
