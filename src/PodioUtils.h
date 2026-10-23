#ifndef __PODIOUTILS__
#define __PODIOUTILS__

#include <cstdint>
#include <string>

namespace PodioUtils
{
  // podio collection ID of a collection: 32-bit MurmurHash3 of its name
  // (as done by podio::CollectionIDTable), used to resolve the
  // collectionID of the podio links (ObjectIDs) without the podio libraries
  uint32_t collectionID(const std::string& name);
}

#endif
