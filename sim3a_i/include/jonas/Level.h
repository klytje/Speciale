#ifndef LEVEL_H
#define LEVEL_H

class Level {
  /*
  * The Level class stores basic information on a nuclear level/state/resonance.
  * The energy should be given as energy above the ground state.
  */

  private:
    double Elevel;
    double Jlevel;
    int Plevel;
 
  public:
    Level(double e = 0.0, double j = 0.0, int p = 1);
    ~Level();

    void SetEnergy(double);
    void SetJ(double);
    void SetParity(int);
  
    double Energy() const;
    double E() const;
    double J() const;
    int Parity() const;
    int Pi() const;
};
#endif
