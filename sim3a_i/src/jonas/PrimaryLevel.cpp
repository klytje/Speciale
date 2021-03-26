#include "jonas/PrimaryLevel.h"

PrimaryLevel::PrimaryLevel() {}

PrimaryLevel::~PrimaryLevel() {}

void PrimaryLevel::SetWidths(Mat<double> input)
{
  widths = input;
}

Mat<double> & PrimaryLevel::GetWidths()
{
  return widths;
}

double PrimaryLevel::GetWidth(int il, int mu)
{
  return widths(il,mu);
}
