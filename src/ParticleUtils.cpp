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
  // (no subscript L, S, mu, tau in the font of the labels; the antiparticles
  // with a combining overline)
  {130, L"K\u2070L"},       // K0L
  {310, L"K\u2070S"},       // K0S
  {311, L"K\u2070"},        // K0
  {-311, L"K\u0305\u2070"}, // anti-K0
  {321, L"K\u207a"},        // K+
  {-321, L"K\u207b"},       // K-
  {2212, L"p"},             // proton
  {-2212, L"p\u0305"},      // antiproton
  {2112, L"n"},             // neutron
  {-2112, L"n\u0305"},      // antineutron
  {3122, L"\u039b"},        // Lambda
  {-3122, L"\u039b\u0305"}, // anti-Lambda
  {12, L"\u03bde"},         // nu_e
  {-12, L"\u03bd\u0305e"},  // anti-nu_e
  {14, L"\u03bd\u03bc"},    // nu_mu
  {-14, L"\u03bd\u0305\u03bc"}, // anti-nu_mu
  {16, L"\u03bd\u03c4"},    // nu_tau
  {-16, L"\u03bd\u0305\u03c4"}  // anti-nu_tau
};
}

std::wstring ParticleUtils::name(int pdg)
{
  auto it = names.find(pdg);

  if (it != names.end())
    return it->second;

  return L"PDG " + std::to_wstring(pdg);
}
