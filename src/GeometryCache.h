/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// GeometryCache: cache of the detector geometry, pruned of the volumes that
// are never drawn.
//
// The geometry files contain all the volumes of the detector (~850k for
// ALLEGRO o1_v03), but only a small fraction of them is drawn: reading,
// building and projecting all of them dominates the startup time.
// The first time a geometry is used with a set of display options, the full
// geometry is read; once the display is set up, the volumes that can be drawn
// in one of the views are determined, and the geometry pruned of the others
// is written to a cache file, which is used the following times.
//
/******************************************************************************/

#ifndef GEOMETRYCACHE_H
#define GEOMETRYCACHE_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class TEveElement;
class TEveGeoShapeExtract;

namespace GeometryCache
{
// the element created by TEveGeoShape::ImportShapeExtract for each extract
using ExtractMap = std::unordered_map<TEveElement*, TEveGeoShapeExtract*>;

// the cached geometry of geomFile for the display options in tag, or
// nullptr if there is none or it is not up to date (the geometry file or
// the code that decides what is drawn, identified by codeVersion, changed)
TEveGeoShapeExtract* read(const std::string& geomFile, const std::string& tag,
                          const std::string& codeVersion);

// map the elements created by TEveGeoShape::ImportShapeExtract(extract) to
// their extract: to be called right after the import, before the elements
// are moved around. The extracts also take a reference to their shapes,
// shared with the elements, so that they stay valid if the elements are
// deleted (they are released by deleteExtract with releaseShapes)
void mapExtracts(TEveGeoShapeExtract* extract, TEveElement* element, ExtractMap& extracts);

// prune the geometry (the extract of the world volume, mapped with
// mapExtracts) of the volumes that cannot be drawn in the views of the display
// (their roots: the geometry in 3D and the scenes of the projections), and
// write it to the cache.
// A volume can be drawn if all its ancestors in a view draw their children
// (then it is drawn, or it can be turned on in the browser).
// The extracts that are pruned are deleted.
void write(const std::string& geomFile, const std::string& tag, const std::string& codeVersion,
           TEveGeoShapeExtract* world, const ExtractMap& extracts,
           const std::vector<TEveElement*>& views);

// delete an extract and its children, but not their shapes, which are
// shared with the elements created from them: with releaseShapes, release the
// references taken by mapExtracts (deleting the shapes no longer used)
void deleteExtract(TEveGeoShapeExtract* extract, bool releaseShapes = false);
}

#endif
