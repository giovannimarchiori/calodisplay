#ifndef __UNITS_H__
#define __UNITS_H__

#include <string>
namespace Units
{
  // units
  // G4 uses mm but Root uses cm..
  inline constexpr double m  = 100.0;
  inline constexpr double cm = 1.0;
  inline constexpr double mm = 0.1;
  inline const std::string units = "cm";
}
#endif
