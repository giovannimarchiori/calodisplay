/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// namespace PodioUtils: utilities for podio files (collection IDs)
//
/******************************************************************************/

#include "PodioUtils.h"

#include <cstring>

namespace
{
uint32_t rotl32(uint32_t x, int8_t r)
{
  return (x << r) | (x >> (32 - r));
}

// MurmurHash3_x86_32 by Austin Appleby (public domain)
uint32_t murmurHash3(const void* key, int len, uint32_t seed)
{
  const uint8_t* data = (const uint8_t*) key;
  const int nblocks = len / 4;
  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;
  uint32_t h1 = seed;

  for (int i = 0; i < nblocks; i++)
  {
    uint32_t k1;
    std::memcpy(&k1, data + 4*i, 4);
    k1 *= c1;
    k1 = rotl32(k1, 15);
    k1 *= c2;
    h1 ^= k1;
    h1 = rotl32(h1, 13);
    h1 = h1*5 + 0xe6546b64;
  }

  const uint8_t* tail = data + nblocks*4;
  uint32_t k1 = 0;
  switch (len & 3)
  {
  case 3: k1 ^= tail[2] << 16; [[fallthrough]];
  case 2: k1 ^= tail[1] << 8;  [[fallthrough]];
  case 1: k1 ^= tail[0];
          k1 *= c1;
          k1 = rotl32(k1, 15);
          k1 *= c2;
          h1 ^= k1;
  }

  h1 ^= len;
  h1 ^= h1 >> 16;
  h1 *= 0x85ebca6b;
  h1 ^= h1 >> 13;
  h1 *= 0xc2b2ae35;
  h1 ^= h1 >> 16;
  return h1;
}
}

uint32_t PodioUtils::collectionID(const std::string& name)
{
  return murmurHash3(name.data(), name.size(), 0);
}
