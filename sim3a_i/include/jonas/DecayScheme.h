#ifndef DECAY_SCHEME_H
#define DECAY_SCHEME_H
#include <vector>
#include <memory>
#include <armadillo>
#include "jonas/Channel.h"
#include "jonas/PrimaryLevel.h"
#include "jonas/SecondaryLevel.h"

using namespace std;
using namespace arma;

class DecayScheme {

  private:
    vector<PrimaryLevel> primLevels;         //Added one by one.
    vector<unique_ptr<Channel>> primChannels; //Constructed at instantiation.
    vector<SecondaryLevel> secLevels;         //Set at instantiation.
    vector<unique_ptr<Channel>> secChannels;  //Constructed at instantiation.

    //Boundary conditions for the decay channels.
    Mat<double> primBoundaries; //Calculated each time a primary level is added.
    vector<double> secBoundaries;  //Calculated at instantiation and if channel radius is changed.

    //unique_ptr<ThreeBodyFunctions> threeBodyFcn; 

    //Updates boundary conditions for both primary and secondary channels.
    void CalculateBoundaries();

    //Checks if a double is zero to within something like machine precision.
    bool IsAlmostZero(double);

    //Calculates the 'density of states' for secondary state specified by the index.
    double DensityFunction(int, double);

  public:
    /**
    * The constructor takes a vector of the l's contributing to the decay of
    * the primary system (first argument) and a vector of levels in the 
    * intermediate system (second argument);
    */
    DecayScheme(vector<int>, vector<SecondaryLevel>, double r1 = 2.47, double r2 = 1.9);
    ~DecayScheme();

    /**
    * Add a primary level to the scheme. The width matrix of the primary level
    * must have dimensions matching the number of l's and secondary levels
    * participating in the decay (they can, however, be set to zero).
    */
    void AddPrimaryLevel(PrimaryLevel);

    /**
    * We allow the user to modify the channel radius for the primary breakup.
    */
    void SetPrimaryRadius(double);

    /**
    * Do. for the secondary breakup.
    */
    void SetSecondaryRadius(double);

    //Simple getter methods...    
    vector<PrimaryLevel> & GetPrimaryLevels();
    vector<unique_ptr<Channel>> & GetPrimaryChannels();
    vector<SecondaryLevel> & GetSecondaryLevels();
    vector<unique_ptr<Channel>> & GetSecondaryChannels();

    /**
    * First parameter is the index of the l (typically 0,1,2 would correspond to l=0,2,4)
    * and the index for the level in the intermediate system.
    */
    double GetPrimaryBoundary(int, int);

    /**
    * Returns the boundary condition given index of l (typically 0,1 would correspond to l=0,2)
    */
    double GetSecondaryBoundary(int);

    //XXX: The three-body functions may be terribly slow in this implementation,
    //     but I think interpolation of the functions is unfeasible.
    /**
    * Calculates the renormalisation factor for the three-body R-matrix functions.
    */
    double Norm(int, double);

    double AvgPenetrability(int, int, double);

    double AvgShiftFunction(int, int, double);
};
#endif
