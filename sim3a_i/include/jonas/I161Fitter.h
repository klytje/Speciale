#ifndef I161FITTER_C
#define I161FITTER_C
#include "KinFitter.h"

class I161Fitter : public KinFitter {
  /**
  * Specific version of the kinematic fitter which minimises the three components
  * of the total momentum in a triple-alpha decay.
  */
  private:

  public:
    I161Fitter();
    ~I161Fitter();
};
#endif
