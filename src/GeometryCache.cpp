/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// GeometryCache: cache of the detector geometry, pruned of the volumes that
// are never drawn
//
/******************************************************************************/

#include "GeometryCache.h"

#include <TEveElement.h>
#include <TEveGeoShapeExtract.h>
#include <TEveProjectionBases.h>
#include <TFile.h>
#include <TGeoShape.h>
#include <TList.h>
#include <TNamed.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <unistd.h>

namespace fs = std::filesystem;

namespace
{
// the cache file of a geometry file: next to it if possible, otherwise in
// the user's cache directory
std::string cachePath(const std::string& geomFile)
{
  fs::path geom(geomFile);
  std::string name = geom.stem().string() + ".displaycache.root";
  fs::path dir = geom.has_parent_path() ? geom.parent_path() : fs::path(".");
  if (access(dir.c_str(), W_OK) == 0)
    return (dir / name).string();
  const char* home = std::getenv("HOME");
  fs::path cacheDir = fs::path(home ? home : ".") / ".cache" / "calodisplay";
  std::error_code ec;
  fs::create_directories(cacheDir, ec);
  return (cacheDir / name).string();
}

// identifies the version of the geometry file and of the code
std::string versionInfo(const std::string& geomFile, const std::string& codeVersion)
{
  std::error_code ec;
  std::ostringstream info;
  info << "size=" << fs::file_size(geomFile, ec)
       << " mtime=" << fs::last_write_time(geomFile, ec).time_since_epoch().count()
       << " code=" << codeVersion;
  return info.str();
}

// the elements that can be drawn: all their ancestors in the view that are
// volumes draw their children (for the projections, the elements that are
// projected). The other elements (e.g. the lists of the subdetectors) can be
// turned on in the browser, so they are not required to draw their children
void collectDrawable(TEveElement* element, const GeometryCache::ExtractMap& extracts,
                     std::unordered_set<TEveElement*>& drawable)
{
  TEveProjected* projected = dynamic_cast<TEveProjected*>(element);
  TEveElement* original = projected ? dynamic_cast<TEveElement*>(projected->GetProjectable()) : element;
  drawable.insert(original);
  if (!element->GetRnrChildren() && extracts.count(original) > 0)
    return;
  for (TEveElement::List_i itr = element->BeginChildren(); itr != element->EndChildren(); ++itr)
    collectDrawable(*itr, extracts, drawable);
}

// remove from the extract the children that are not kept and have no kept
// descendants (and delete them); returns true if the extract is kept
bool prune(TEveGeoShapeExtract* extract, const std::unordered_set<TEveGeoShapeExtract*>& keep, long& nKept)
{
  bool keepChildren = false;
  if (extract->HasElements())
  {
    std::vector<TEveGeoShapeExtract*> kept;
    for (TObject* obj : *extract->GetElements())
    {
      TEveGeoShapeExtract* child = static_cast<TEveGeoShapeExtract*>(obj);
      if (prune(child, keep, nKept))
        kept.push_back(child);
      else
        GeometryCache::deleteExtract(child, true);
    }
    extract->GetElements()->Clear("nodelete");
    for (TEveGeoShapeExtract* child : kept)
      extract->GetElements()->Add(child);
    keepChildren = !kept.empty();
  }
  bool keepThis = keepChildren || keep.count(extract) > 0;
  if (keepThis)
    nKept++;
  return keepThis;
}

long countExtracts(TEveGeoShapeExtract* extract)
{
  long n = 1;
  if (extract->HasElements())
    for (TObject* obj : *extract->GetElements())
      n += countExtracts(static_cast<TEveGeoShapeExtract*>(obj));
  return n;
}
}

TEveGeoShapeExtract* GeometryCache::read(const std::string& geomFile, const std::string& tag,
                                         const std::string& codeVersion)
{
  std::string path = cachePath(geomFile);
  if (!fs::exists(path))
    return nullptr;
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie())
    return nullptr;
  TNamed* info = dynamic_cast<TNamed*>(file->Get(("info_" + tag).c_str()));
  if (!info || versionInfo(geomFile, codeVersion) != info->GetTitle())
    return nullptr;
  TEveGeoShapeExtract* world = dynamic_cast<TEveGeoShapeExtract*>(file->Get(("world_" + tag).c_str()));
  if (world)
    std::cout << "Reading the geometry (without the volumes that are never drawn) from the cache " << path << std::endl;
  return world;
}

void GeometryCache::mapExtracts(TEveGeoShapeExtract* extract, TEveElement* element, ExtractMap& extracts)
{
  // reference counting of the shapes as in TEveGeoShape::SetShape
  if (TGeoShape* shape = extract->GetShape())
    shape->SetUniqueID(shape->GetUniqueID() + 1);
  if (element)
    extracts[element] = extract;
  if (!extract->HasElements())
    return;
  // ImportShapeExtract creates the children in the order of the extracts
  // (the shapes of all the extracts are referenced, even without element)
  TEveElement::List_i itr;
  if (element)
    itr = element->BeginChildren();
  for (TObject* obj : *extract->GetElements())
  {
    TEveElement* child = nullptr;
    if (element && itr != element->EndChildren())
      child = *(itr++);
    mapExtracts(static_cast<TEveGeoShapeExtract*>(obj), child, extracts);
  }
}

void GeometryCache::write(const std::string& geomFile, const std::string& tag, const std::string& codeVersion,
                          TEveGeoShapeExtract* world, const ExtractMap& extracts,
                          const std::vector<TEveElement*>& views)
{
  // the volumes that can be drawn
  std::unordered_set<TEveElement*> drawable;
  for (TEveElement* view : views)
    if (view)
      collectDrawable(view, extracts, drawable);
  std::unordered_set<TEveGeoShapeExtract*> keep;
  for (TEveElement* element : drawable)
  {
    auto it = extracts.find(element);
    if (it != extracts.end())
      keep.insert(it->second);
  }

  long nTotal = countExtracts(world);
  long nKept = 0;
  prune(world, keep, nKept);

  std::string path = cachePath(geomFile);
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "UPDATE"));
  if (!file || file->IsZombie())
  {
    std::cout << "WARNING: cannot write the geometry cache " << path << std::endl;
    return;
  }
  world->Write(("world_" + tag).c_str(), TObject::kOverwrite);
  TNamed info(("info_" + tag).c_str(), versionInfo(geomFile, codeVersion).c_str());
  info.Write(nullptr, TObject::kOverwrite);
  file->Close();
  std::cout << "Geometry cache written to " << path << " (" << nKept << " of " << nTotal
            << " volumes can be drawn): the display will start faster next time" << std::endl;
}

void GeometryCache::deleteExtract(TEveGeoShapeExtract* extract, bool releaseShapes)
{
  if (extract->HasElements())
  {
    for (TObject* obj : *extract->GetElements())
      deleteExtract(static_cast<TEveGeoShapeExtract*>(obj), releaseShapes);
    extract->GetElements()->Clear("nodelete");
  }
  TGeoShape* shape = extract->GetShape();
  if (releaseShapes && shape)
  {
    shape->SetUniqueID(shape->GetUniqueID() - 1);
    if (shape->GetUniqueID() == 0)
      delete shape;
  }
  extract->SetShape(nullptr);
  delete extract;
}
