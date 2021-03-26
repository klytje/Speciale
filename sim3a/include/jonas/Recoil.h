#ifndef RECOIL_H
#define RECOIL_H
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TRandom3.h>

using namespace std;
using namespace TMath;

class Recoil {
  /**
  * Class to generate beta-neutrino recoil momentum for beta decays. In order to do
  * this efficiently, the class divides the Q-value range into bins and samples the
  * recoil distributions in the relevant bin.
  */

  private:
    double Qmax;    //Maximum possible Q-value in keV.
    double Wmax;
    vector<TF1> fTbeta;   //Distributions of electron kinetic energies for different Q-value bins.
    vector<TF1> fCos;     //Distributions of angle between electron and neutrino momenta.
    TRandom3 rGen;

    static constexpr double mBeta = 510.9989461;

    int FindQbin(double);
    int FindWbin(double);
    int NBins();  //Number of bins on the energy 'grid'. More bins give more accurate behaviour.

  public:
    /**
    * The constructor passes the arguments on to the Initialise()-function.
    */
    Recoil(double q = 1000., int n = 1000);
    ~Recoil();

    /**
    * Initialise() takes the maximum Q-value of the beta-decay and the desired number of bins
    * on the energy grid as arguments.
    */   
    void Initialise(double, int);
  
    /**
    * Generate() takes the actual Q-value of the decay as argument and returns recoil momentum
    * in units of keV/c. If provided with the two last arguments it will also tell the
    * four-momenta of the generated beta and neutrino.
    */    
    TVector3 Generate(double, TLorentzVector *pBeta = 0, TLorentzVector *pNeutrino = 0);
};
#endif
