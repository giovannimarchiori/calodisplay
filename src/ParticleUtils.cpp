/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// namespace ParticleUtils: names of the particles from their PDG ID
//
/******************************************************************************/

#include "ParticleUtils.h"

#include <unordered_map>

namespace
{
const std::unordered_map<int, std::wstring> names =
{
  {111, L"\u03c0\u2070"},  // pi0
  {211, L"\u03c0\u207a"},  // pi+
  {-211, L"\u03c0\u207b"}, // pi-
  {22, L"\u03b3"},          // gamma
  {11, L"e\u207b"},         // e-
  {-11, L"e\u207a"},        // e+
  {13, L"\u03bc\u207b"},   // mu-
  {-13, L"\u03bc\u207a"},  // mu+
  {130, L"K\u2070L"}        // K0L (no subscript L in the font of the labels)
};
}

std::wstring ParticleUtils::name(int pdg)
{
  auto it = names.find(pdg);

  if (it != names.end())
    return it->second;

  return L"PDG " + std::to_wstring(pdg);
}
