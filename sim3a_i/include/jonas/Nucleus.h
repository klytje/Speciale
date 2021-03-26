#ifndef NUCLEUS_H
#define NUCLEUS_H
#include <vector>
#include "Level.h"

using namespace std;

class Nucleus {
  /**
  * The Nucleus object is a container to hold basic information describing
  * a nucleus. It holds a vector of Level-objects, of which the first one
  * (index 0) is the ground state.
  */

  private:
    int a;
    double m;
    int q;
    vector<Level> levels;

  public:
    /**
    * The default constructor will produce a 4He-nucleus.
    */
    Nucleus(int A = 4, int Z = 2, double J = 0., int parity = 1);
    ~Nucleus();

    int AddLevel(Level&);
    const vector<Level> & Levels();

    void SetMass(double);
    void SetJ(double);
    void SetCharge(int);
    void SetParity(int);
    
    double Mass() const;
    double M() const;
    int Charge() const;
    int Q() const;
    double J() const;
    int Parity() const;
    int Pi() const;
    int A() const;
    int Z() const;

    /**
    * Generates a unique integer for each A,Z using the Cantor pairing function.
    * Can be used to see if two instances of the class describes the same nucleus.
    */
    int IsotopeID();
};
#endif
