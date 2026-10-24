/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// namespace ParticleUtils: names of the particles from their PDG ID
//
/******************************************************************************/

#ifndef __PARTICLEUTILS__
#define __PARTICLEUTILS__

#include <string>

namespace ParticleUtils
{
  // name of a particle, with greek letters and superscripts (unicode),
  // or its PDG code if it is not known
  std::wstring name(int pdg);
}

#endif
