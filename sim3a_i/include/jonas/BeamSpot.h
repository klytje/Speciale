#ifndef BEAMSPOT_H
#define BEAMSPOT_H
#include <TF1.h>
#include <TF2.h>
#include <TVector3.h>
#include <memory>

using namespace std;

class BeamSpot {
  /**
  * A class that, given an implantation depth distribution (from SRIM) and a transverse
  * implantation point distribution (2D), is able to generate random decay points for
  * a Monte Carlo simulation.
  */

  private:
    unique_ptr<TF1> depthProfile;     //From SRIM, (i.e. units of Ångstrøm).
    unique_ptr<TF2> transverseProfile;
    double foilThickness;
    double GetDepth();

  public:
    /**
    * The depth and transverse distributions are specified at construction. The implantation depth
    * distribution is expected to come directly from SRIM, and should give the distribution in
    * units of Ångstrøm. The transverse profile should be in mm, though. The third argument
    * is the thickness of the foil (should also be given in mm)
    */
    BeamSpot(unique_ptr<TF1>, unique_ptr<TF2>, double);
    ~BeamSpot();
    
    /**
    * Generate a random decay point according to the specified distributions. All positions
    * are return in units of mm.
    */
    TVector3 GetDecayPoint();

    /**
    * Generate a random decay point according to the specified distributions. All positions
    * are return in units of mm.
    */
    void GetDecayPoint(TVector3 *);
};
#endif
