#include "ParticleUtils.h"

#include <unordered_map>

namespace
{
const std::unordered_map<int, const char*> names =
{
  {111,"pi0"},
  {211,"pi+"},
  {-211,"pi-"},
  {22,"γ"},
  {11,"e-"},
  {-11,"e+"},
  {13,"mu-"},
  {-13,"mu+"}
};
}

const char* ParticleUtils::name(int pdg)
{
  auto it = names.find(pdg);

  if (it != names.end())
    return it->second;

  return "unknown";
}
