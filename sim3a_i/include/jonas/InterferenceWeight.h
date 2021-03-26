#ifndef INTERFERENCE_WEIGHT_H
#define INTERFERENCE_WEIGHT_H
#include <vector>
#include "DecayWeight.h"
#include "Level.h"
#include "Channel.h"

using namespace std;

class InterferenceWeight : public DecayWeight {

  private:
    static constexpr double hbarc = 197326.9788;       //Planck const x speed of light [keV * fm]
//    vector<vector<int>> contributions;  //Decays models {Ja,L,Jb}.
    vector<double> weights;             //Integrals of the decay models.
    vector<Level> primLevels;
    vector<unique_ptr<Channel>> primChannels;
    vector<Level> secLevels;
    vector<unique_ptr<Channel>> secChannels;

    //Following variables are for the final-state Coulomb interaction.
    bool doCorrection;
    vector<unique_ptr<Channel>> c1, c2;

    /**
    * We calculate the Clebsch-Gordan coefficients from the Wigner 3j-symbols instead
    * of using the otherwise excellent class provided in the present library.
    */
    double ClebschGordan(double,double,double,double,double,double);

    /**
    * Calculates the radius for the final state Coulomb correction (the 'smudge factor')
    * based on the lifetime of the intermediate resonance calculated from the resonant
    * phase shift (actually its derivative).
    */
    double CorrectionRadius(unique_ptr<Channel> &, unique_ptr<Channel> &, Level &, double, double);

  public:
    InterferenceWeight();
    ~InterferenceWeight();

    //Final state Coulomb correction, input in fm.
    void DoCorrection(bool input = true, double radius = 15.);

    /**
    * Set the various decay modes and the associated weights. Each mode is
    * specified with a vector containing integers {Ja, L, Jb}.
    */
    void SetContributions(vector<vector<int>>, vector<double>, double);

    /**
    * Specify whether to use the fast interpolation method when calculating the
    * R-matrix functions.
    */
    void UseInterpolation(bool input = true);

    double Calculate(vector<TLorentzVector> &);
};
#endif
