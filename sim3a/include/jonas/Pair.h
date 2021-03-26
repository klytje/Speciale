#ifndef PAIR_H
#define PAIR_H
#include <vector>
#include "Nucleus.h"

using namespace std;

class Pair {
  /**
  * A class which holds two Nucleus-objects.
  */

  private:
    vector<Nucleus> pair;

    static constexpr double hbarc = 197326.9788;       //Planck const x speed of light [keV * fm]
    static constexpr double alpha = 7.2973525664e-3;   //Fine structure constant ~1/137.

  public:
    /**
    * The constructor takes the two Nucleus-objects that will make up the pair.
    */
    Pair(Nucleus n1 = Nucleus(), Nucleus n2 = Nucleus());
    ~Pair();

    /**
    * It is possible to change the consituents of the pair.
    */
    void SetPair(Nucleus, Nucleus);
   
    /**
    * Retrieve the constituents in order to learn more about them.
    */
    const vector<Nucleus> & Particles() const;

    /**
    * The following functions provide easy acces to often needed quantities.
    * Masses are in keV and charges in unit of the elementary charge.
    */
    double RedMass() const;
    double QProduct() const;
    double Mtot() const;
    int Qtot() const;
    
    /**
    * Returns the wave number in [1/fm] for a relative kinetic energy in [keV].
    */
    double WaveNumber(double);

    /**
    * Calculates the Sommerfeld parameter for a given relative kinetic energy in [keV].
    */
    double Eta(double);

    /**
    * Calculates the Gamow factor exp[-2*pi*eta] for a relative kinetic energy in [keV].
    */
    double GamowFactor(double);

    /**
    * Returns the number (A1^(1/3) + A2^(1/3)) which is a measure of the 'radius' of the
    * pair in units of the nucleon radius.
    */
    double Radius() const;    
};
#endif
