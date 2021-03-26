#ifndef SECONDARY_LEVEL_H
#define SECONDARY_LEVEL_H
#include <vector>

using namespace std;

struct SecondaryLevel {
  double El;
  int J;
  vector<double> widths;  //Formal reduced widths.
};

#endif
