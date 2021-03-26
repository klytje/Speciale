#ifndef MULTI_LEVEL_WEIGHT_H
#define MULTI_LEVEL_WEIGHT_H
#include <vector>
#include <memory>
#include "DecayWeight.h"
#include "Level.h"
#include "Channel.h"

using namespace std;

class MultiLevelWeight : public DecayWeight {

  private:
    Level primLevel;
    unique_ptr<Channel> primChannel;
    Level secLevel;
    unique_ptr<Channel> secChannel;

    /**
    * We calculate the Clebsch-Gordan coefficients from the Wigner 3j-symbols instead
    * of using the otherwise excellent class provided in the present library.
    */
    double ClebschGordan(double,double,double,double,double,double);

  public:
    MultiLevelWeight();
    ~MultiLevelWeight();

    void SetPrimary(Level&, unique_ptr<Channel>);
    void SetSecondary(Level&, unique_ptr<Channel>);

    double Calculate(vector<TLorentzVector> &);
};
#endif
