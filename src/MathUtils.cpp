/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// namespace MathUtils: mathematical utility functions
//
/******************************************************************************/

#include "MathUtils.h"
#include <cmath>

namespace MathUtils
{
int sign(float value)
{
  return (value > 0.f) - (value < 0.f);
}

double solveQuadratic(double a,
                      double b,
                      double c,
                      int solution)
{
  if (a == 0.)
    return -c / b;

  double d = b*b - 4*a*c;
  if (d < 0.)
    return -9.999999999e9;

  d = std::sqrt(d);
  if (solution > 0)
    return (-b + d)/(2*a);
  return (-b - d)/(2*a);
}
}
