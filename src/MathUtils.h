/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// namespace MathUtils: mathematical utility functions
//
/******************************************************************************/

#ifndef __MATHUTILS__
#define __MATHUTILS__

namespace MathUtils
{
  int sign(float value);
  double solveQuadratic(double a,
                        double b,
                        double c,
                        int solution);
}

#endif
