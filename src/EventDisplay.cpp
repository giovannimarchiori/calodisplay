/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventDisplay: creates and populates the event display
//
/******************************************************************************/

#include "EventDisplay.h"
#include "MathUtils.h"
#include "ParticleUtils.h"
#include "DetectorGeometry.h"
#include "MagField.h"
#include "CylField.h"
#include "Units.h"
#include "PathUtils.h"
// (generated in build/ by scripts/geometry_code_version.sh)
#include "GeometryCodeVersion.h"

// ROOT core
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTimer.h>
#include <TApplication.h>
#include <TTimeStamp.h>

// TG
#include <TGButton.h>
#include <TGTab.h>
#include <TGMsgBox.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGStatusBar.h>
#include <TVirtualX.h>

// TGeo
#include <TGeoTube.h>
#include <TGeoMatrix.h>

// TEve
#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEveProjectionAxes.h>
#include <TEveRGBAPaletteOverlay.h>
#ifdef __APPLE__
// OpenGL is deprecated on macOS (but still available, and used by ROOT)
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <TEveGeoNode.h>
#include <TEveGeoShapeExtract.h>
#include <TEvePolygonSetProjected.h>
#include <TEveTrack.h>
#include <TGLOrthoCamera.h>
#include <TGLRnrCtx.h>
#include <TGLPerspectiveCamera.h>


#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <iostream>
#include <sstream>

using Units::units;
using Units::cm;
using Units::mm;

namespace
{
// name of the main tracker, from the name of its collection of hits or digis
// (e.g. DCHCollection, STTDigis)
std::string mainTrackerName(const std::string &collection)
{
  if (collection.find("STT") != std::string::npos)
    return "STT";
  if (collection.find("DCH") != std::string::npos)
    return "DCH";
  return "Main tracker";
}

// bin of an ECAL endcap cell in the rho-z projection, from its wheel, rho bin,
// side and sign of y (the cells on the two sides and above and below the beam
// axis are drawn separately)
int ecalEndcapRhoZBin(ULong_t cellID, float y)
{
  int bin = (int)DetectorGeometry::ECalEndCapWheel(cellID) * 256 + (int)DetectorGeometry::ECalEndCapRhoBin(cellID);
  return (bin * 2 + (DetectorGeometry::ECalEndCapSide(cellID) > 0)) * 2 + (y > 0);
}

// version of the code that decides which volumes of the G4 geometry are drawn
// (a hash of the functions listed in scripts/geometry_code_version.sh): the
// geometry caches written by the other versions are rebuilt
const std::string geometryCodeVersion = GEOMETRY_CODE_VERSION;

// vertical legend of a palette, with the axis on the right of the colour bar,
// on a semi-transparent dark background so that it can be read on top of the
// detector (the background also covers the title, drawn separately above it).
// The palette values are integers, shown on the axis multiplied by scale.
class PaletteLegend : public TEveRGBAPaletteOverlay
{
public:
  PaletteLegend(TEveRGBAPalette *p, Float_t posx, Float_t posy, Float_t width, Float_t height, Double_t scale)
      : TEveRGBAPaletteOverlay(p, posx, posy, width, height), fScale(scale)
  {
    fAxis.SetNdivisions(505);
    fAxisPainter.SetLabelAlign(TGLFont::kLeft, TGLFont::kCenterV);
  }

  void Render(TGLRnrCtx &rnrCtx) override
  {
    const Int_t vMin = fPalette->GetMinVal();
    const Int_t vMax = fPalette->GetMaxVal();
    if (vMin >= vMax)
      return;

    glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT | GL_CURRENT_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // [0,1] window units
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, 0, 1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // background, including the axis labels and the title
    glColor4f(0.f, 0.f, 0.f, 0.7f);
    glBegin(GL_QUADS);
    glVertex2f(fPosX - 0.015f, fPosY - 0.03f);
    glVertex2f(fPosX + fWidth + 0.07f, fPosY - 0.03f);
    glVertex2f(fPosX + fWidth + 0.07f, fPosY + fHeight + 0.06f);
    glVertex2f(fPosX - 0.015f, fPosY + fHeight + 0.06f);
    glEnd();

    // colour bar, from the minimum (bottom) to the maximum (top)
    const Int_t step = std::max(1, (vMax - vMin) / 256);
    glBegin(GL_QUAD_STRIP);
    for (Int_t v = vMin;; v += step)
    {
      v = std::min(v, vMax);
      const Float_t y = fPosY + fHeight * (v - vMin) / (vMax - vMin);
      glColor4ubv(fPalette->ColorFromValue(v));
      glVertex2f(fPosX, y);
      glVertex2f(fPosX + fWidth, y);
      if (v == vMax)
        break;
    }
    glEnd();

    // frame and axis
    glColor4f(1.f, 1.f, 1.f, 1.f);
    glBegin(GL_LINE_LOOP);
    glVertex2f(fPosX, fPosY);
    glVertex2f(fPosX + fWidth, fPosY);
    glVertex2f(fPosX + fWidth, fPosY + fHeight);
    glVertex2f(fPosX, fPosY + fHeight);
    glEnd();
    // the axis painter draws along x, with the ticks towards -y: rotate it by
    // 90 degrees, so that it is vertical with the ticks and labels on the right
    const Double_t low = vMin * fScale;
    const Double_t high = vMax * fScale;
    fAxis.SetLimits(low, high);
    fAxis.SetTickLength(0.3 * fWidth);
    glTranslatef(fPosX + fWidth, fPosY, 0);
    glRotatef(90, 0, 0, 1);
    glScalef(fHeight / (high - low), 1, 1);
    glTranslatef(-low, 0, 0);
    fAxisPainter.RefTMOff(0).Set(0, -1, 0);
    fAxisPainter.PaintAxis(rnrCtx, &fAxis);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glPopAttrib();
  }

private:
  Double_t fScale;
};

// palette for the energies of cells, with values in MeV
TEveRGBAPalette *makeCellPalette()
{
  TEveRGBAPalette *palette = new TEveRGBAPalette(0, 1000);
  // ref-counted by the digit sets that use it: keep a reference so that it
  // is not deleted when the clusters of an event are destroyed
  palette->IncRefCount();
  return palette;
}

// set the range of a cell palette from 0 to the maximum value, rounded up to
// 1, 2 or 5 times a power of 10 (default range if there are no values)
void setCellPaletteRange(TEveRGBAPalette *palette, int maxValue)
{
  int high = 1000;
  if (maxValue > 0)
  {
    int power = 1;
    while (power * 10 <= maxValue)
      power *= 10;
    high = (maxValue <= power) ? power : (maxValue <= 2 * power) ? 2 * power : (maxValue <= 5 * power) ? 5 * power : 10 * power;
  }
  palette->SetLimits(0, high);
  palette->SetMinMax(0, high);
  // the legend is drawn only if the colour array is set up
  palette->SetupColorArray();
}

// largest value of the digits of a digit set
int maxDigitValue(TEveDigitSet *digits)
{
  int maxValue = 0;
  if (digits)
    for (int i = 0; i < digits->GetPlex()->Size(); i++)
      maxValue = std::max(maxValue, digits->GetDigit(i)->fValue);
  return maxValue;
}

// identifies a cell of a cluster in a 2D projection: the system, the index of
// the cluster, the layer, and the bin in the projected view (theta bin for
// rho-z; module for the ECAL and phi bin for the HCAL for rho-phi)
struct ProjectedCellKey
{
  int system;
  int cluster;
  int layer;
  int bin;
  bool operator<(const ProjectedCellKey &other) const
  {
    return std::tie(system, cluster, layer, bin) < std::tie(other.system, other.cluster, other.layer, other.bin);
  }
};

// find the first element (depth-first) whose name begins with prefix
TEveElement *findElementBeginningWith(TEveElement *el, const char *prefix)
{
  for (TEveElement::List_i itr = el->BeginChildren(); itr != el->EndChildren(); ++itr)
  {
    if (TString((*itr)->GetElementName()).BeginsWith(prefix))
      return *itr;
    if (TEveElement *found = findElementBeginningWith(*itr, prefix))
      return found;
  }
  return nullptr;
}

// count the ECAL barrel readout layers in the G4 geometry
// (children named layer* of the active volume in the LAr bath)
// returns -1 if the active volume is not found
int countECalBarrelLayers(TEveElement *world)
{
  TEveElement *ecalBarrel = findElementBeginningWith(world, "ECalBarrel");
  if (!ecalBarrel)
    return -1;
  TEveElement *active = findElementBeginningWith(ecalBarrel, "active");
  if (!active)
    return -1;
  int nLayers = 0;
  for (TEveElement::List_i itr = active->BeginChildren(); itr != active->EndChildren(); ++itr)
    if (TString((*itr)->GetElementName()).BeginsWith("layer"))
      nLayers++;
  return nLayers;
}
}


void EventDisplay::fillClusters(ClusterSet &set)
{
  // the clusters and their cells in the current event
  const EventReader::Clusters &clusters = *set.clusters;
  const EventReader::ClusterCells &cells = clusters.cells;
  // read once, not in the loop over the clusters
  const bool calculateBaryCenters = displayConfig.getBoolConfig("calculateClusterBaryCenters");

  std::cout << "Creating " << set.label << " clusters" << std::endl;
  std::vector<std::unique_ptr<CaloCluster>> *clusterData = &set.data;
  unsigned int nClusters = clusters.size();

  if (debug)
  {
    std::cout << "  n(clusters) = " << nClusters << std::endl;
    std::cout << "  Looping over clusters to fill shape variables" << std::endl;
  }

  // delete clusters from previous event
  if (clusterData->size() > 0)
  {
    if (debug)
    {
      std::cout << "  Deleting old cluster data" << std::endl;
    }
    clusterData->clear();
  }

  // the number of layers of each part of the calorimeters
  std::array<int, CaloCluster::kNParts> nLayers;
  nLayers[CaloCluster::kECalBarrel] = geomReader->nLayers;
  nLayers[CaloCluster::kECalEndcap] = geomReader->nLayersECalEndCap;
  nLayers[CaloCluster::kHCalBarrel] = geomReader->nLayersHCal;
  nLayers[CaloCluster::kHCalEndcap] = geomReader->nLayersHCalEndCap;
  nLayers[CaloCluster::kMuonBarrel] = geomReader->nLayersMuonBarrel;
  nLayers[CaloCluster::kMuonEndcap] = geomReader->nLayersMuonEndcap;
  // the part of the calorimeters and the layer of a cell (false if the system is unknown)
  auto cellPartAndLayer = [this](ULong_t cellID, CaloCluster::Part &part, ULong_t &layer) -> bool
  {
    ULong_t systemID = geomReader->SystemID(cellID);
    if (systemID == DetectorGeometry::kECalBarrel)
    {
      part = CaloCluster::kECalBarrel;
      layer = geomReader->ECalBarrelLayer(cellID);
    }
    else if (systemID == DetectorGeometry::kHCalBarrel)
    {
      part = CaloCluster::kHCalBarrel;
      layer = geomReader->HCalBarrelLayer(cellID);
    }
    else if (systemID == DetectorGeometry::kECalEndcap)
    {
      part = CaloCluster::kECalEndcap;
      layer = geomReader->ECalEndCapZBin(cellID);
    }
    else if (systemID == DetectorGeometry::kHCalEndcap)
    {
      part = CaloCluster::kHCalEndcap;
      layer = geomReader->HCalEndCapLayer(cellID);
    }
    else if (systemID == DetectorGeometry::kMuonBarrel)
    {
      part = CaloCluster::kMuonBarrel;
      layer = geomReader->MuonBarrelLayer(cellID);
    }
    else if (systemID == DetectorGeometry::kMuonEndcap)
    {
      part = CaloCluster::kMuonEndcap;
      layer = geomReader->MuonEndcapLayer(cellID);
    }
    else
      return false;
    return true;
  };

  // loop over clusters and fill the relevant info
  for (unsigned int i = 0; i < nClusters; i++)
  {
    // keep only those above threshold
    float E = (*clusters.energy)[i];
    if (E < ClusterEnergyThreshold)
    {
      if (debug)
      {
        std::cout << "  Cluster energy " << E << " is below threshold, skipping .." << std::endl;
      }
      continue;
    }
    // create empty cluster object
    clusterData->emplace_back(std::make_unique<CaloCluster>());
    CaloCluster *cluster = clusterData->back().get();

    // set index and energy
    cluster->setIndex(i);
    cluster->setEnergy(E);

    // set barycenter
    if (debug)
      std::cout << "  Setting cluster barycenter position .." << std::endl;
    cluster->setBarycenterXYZ((*clusters.x)[i] * mm,
                              (*clusters.y)[i] * mm,
                              (*clusters.z)[i] * mm);

    // set cluster barycenters layer-by-layer

    // - recalculation from cells
    if (calculateBaryCenters) {
      // retrieve index of first and last+1 cell in cell vectors
      if (debug)
      {
        std::cout << "  Retrieving cluster cells and calculating energy per layer .." << std::endl;
      }
      // (no cells if the cells of the clusters are not available)
      const unsigned int nCells = cells.size();
      unsigned int first_hit = (*clusters.hits_begin)[i];
      unsigned int last_hit = (*clusters.hits_end)[i];

      // loop over cells to calculate energy per layer
      std::array<std::vector<float>, CaloCluster::kNParts> energyVsLayer;
      for (int part = 0; part < CaloCluster::kNParts; part++)
        energyVsLayer[part].assign(nLayers[part], 0.);
      for (unsigned int iCell = first_hit; iCell < std::min(last_hit, nCells); iCell++)
      {
        ULong_t cellID = cells.cellID[iCell];
        float energy = cells.energy[iCell];
        if (cellID == 0)
          continue; // cell not available (see EventReader::fillClusterCells)
        CaloCluster::Part part;
        ULong_t layer = 0;
        if (!cellPartAndLayer(cellID, part, layer))
        {
          std::cout << "WARNING: unknown system " << geomReader->SystemID(cellID) << " for cell " << cellID << ", skipping cell" << std::endl;
          continue;
        }
        if (layer >= energyVsLayer[part].size())
        {
          std::cout << "WARNING: layer " << layer << " of system " << geomReader->SystemID(cellID) << " out of range [0, " << energyVsLayer[part].size() << "), skipping cell" << std::endl;
          continue;
        }
        energyVsLayer[part][layer] += energy;
      }
      for (int part = 0; part < CaloCluster::kNParts; part++)
        cluster->setEnergyVsLayers(CaloCluster::Part(part), energyVsLayer[part]);

      // loop again over cells to calculate barycenter per layer
      if (debug)
      {
        std::cout << "  Calculating cell barycenter per layer .." << std::endl;
      }
      std::array<std::vector<TVector3>, CaloCluster::kNParts> barycenterVsLayer;
      std::array<std::vector<float>, CaloCluster::kNParts> sumWeightsVsLayer;
      for (int part = 0; part < CaloCluster::kNParts; part++)
      {
        barycenterVsLayer[part].resize(nLayers[part]);
        sumWeightsVsLayer[part].assign(nLayers[part], 0.0);
      }
      for (unsigned int iCell = first_hit; iCell < std::min(last_hit, nCells); iCell++)
      {
        ULong_t cellID = cells.cellID[iCell];
        float energy = cells.energy[iCell];
        float x = cells.x[iCell] * mm;
        float y = cells.y[iCell] * mm;
        float z = cells.z[iCell] * mm;
        CaloCluster::Part part;
        ULong_t layer = 0;
        // the muon tagger and unknown systems are not used for the barycenter;
        // out-of-range layers were already reported in the first loop
        if (!cellPartAndLayer(cellID, part, layer) || part == CaloCluster::kMuonBarrel || part == CaloCluster::kMuonEndcap ||
            layer >= energyVsLayer[part].size())
          continue;
        if (energyVsLayer[part][layer] > 0)
        {
          float weight = energy / energyVsLayer[part][layer];
          TVector3 &barycenter = barycenterVsLayer[part][layer];
          barycenter.SetXYZ(barycenter.X() + x * weight,
                            barycenter.Y() + y * weight,
                            barycenter.Z() + z * weight);
          sumWeightsVsLayer[part][layer] += weight;
        }
      }

      for (int part = 0; part < CaloCluster::kNParts; part++)
      {
        if (part != CaloCluster::kMuonBarrel && part != CaloCluster::kMuonEndcap)
        {
          for (int iLayer = 0; iLayer < nLayers[part]; iLayer++)
          {
            if (energyVsLayer[part][iLayer] > 0)
            {
              TVector3 &barycenter = barycenterVsLayer[part][iLayer];
              barycenter.SetXYZ(barycenter.X() / sumWeightsVsLayer[part][iLayer],
                                barycenter.Y() / sumWeightsVsLayer[part][iLayer],
                                barycenter.Z() / sumWeightsVsLayer[part][iLayer]);
            }
          }
        }
        cluster->setBarycenterVsLayers(CaloCluster::Part(part), barycenterVsLayer[part]);
      }
    }

    if (debug) {
      cluster->print();
    }
  }
}

void EventDisplay::drawParticles()
{
  if (displayConfig.getBoolConfig("drawGenParticles"))
  {
    std::cout << "Creating particles" << std::endl;
    if (particles == nullptr)
    {
      if (ParticleEnergyThreshold>0.)
        particles = new TEveTrackList(Form("particles (p>%.2f GeV)", ParticleEnergyThreshold));
      else
        particles = new TEveTrackList("particles");
      TEveTrackPropagator *trkProp = particles->GetPropagator();
      trkProp->SetMagFieldObj(magField.get(), false);
      trkProp->SetMaxR(rMax);
      // trkProp->SetMaxZ(geomReader->zMax);
      trkProp->SetMaxZ(geomReader->zMaxEndCap);
      trkProp->SetRnrReferences(false);
      trkProp->SetFitReferences(false);
      particles->SetMainColor(kWhite);
      particles->SetLineWidth(1.);
      particles->SetLineStyle(5);
      gEve->AddElement(particles);
    }
    else
      particles->DestroyElements();

    // loop over the MC particles to create the particle tracks
    for (unsigned int ip = 0; ip < eventReader->genParticles.generatorStatus->GetSize(); ip++)
    {
      // keep only the final-state particles of the generator (status 1) and
      // the particles created in the simulation (status 0), not the
      // intermediate particles of the generator (status > 1)
      if ((*eventReader->genParticles.generatorStatus)[ip] > 1)
        continue;

      float m = (*eventReader->genParticles.mass)[ip];
      float px = (*eventReader->genParticles.momentum_x)[ip];
      float py = (*eventReader->genParticles.momentum_y)[ip];
      float pz = (*eventReader->genParticles.momentum_z)[ip];
      float p = sqrt(px * px + py * py + pz * pz);
      if (p < ParticleEnergyThreshold)
        continue;

      double t = (*eventReader->genParticles.time)[ip];
      double x1 = (*eventReader->genParticles.vertex_x)[ip] * mm;
      double y1 = (*eventReader->genParticles.vertex_y)[ip] * mm;
      double z1 = (*eventReader->genParticles.vertex_z)[ip] * mm;

      TEveMCTrack mct;
      mct.SetPdgCode((*eventReader->genParticles.PDG)[ip]);
      mct.SetMomentum(px, py, pz, sqrt(p * p + m * m));
      mct.SetProductionVertex(x1, y1, z1, t);

      TEveVectorF v;
      v[0] = (*eventReader->genParticles.endpoint_x)[ip] * mm;
      v[1] = (*eventReader->genParticles.endpoint_y)[ip] * mm;
      v[2] = (*eventReader->genParticles.endpoint_z)[ip] * mm;

      // for electrons, for which there could be some hard bremsstrahlung, we
      // split the trajectory into pieces between photon emissions
      if (std::abs((*eventReader->genParticles.PDG)[ip]) == 11)
      {
        // retrieve the various splittings
        auto segments = getElectronBremsstrahlungSegments(
                                                          ip,
                                                          *eventReader->genParticles.daughters_index,
                                                          *eventReader->genParticles.PDG,
                                                          *eventReader->genParticles.vertex_x,
                                                          *eventReader->genParticles.vertex_y,
                                                          *eventReader->genParticles.vertex_z,
                                                          *eventReader->genParticles.momentum_x,
                                                          *eventReader->genParticles.momentum_y,
                                                          *eventReader->genParticles.momentum_z,
                                                          *eventReader->genParticles.daughters_begin,
                                                          *eventReader->genParticles.daughters_end);
        // create a track for each segment
        for(const auto& seg : segments)
        {
          TEveMCTrack mctSegment;
          mctSegment.SetPdgCode((*eventReader->genParticles.PDG)[ip]);
          double pSeg = seg.momentum.Mag();
          double pxSeg = seg.momentum.X();
          double pySeg = seg.momentum.Y();
          double pzSeg = seg.momentum.Z();
          double xSeg = seg.start.X();
          double ySeg = seg.start.Y();
          double zSeg = seg.start.Z();
          double mSeg = (*eventReader->genParticles.mass)[ip];
          mctSegment.SetMomentum(pxSeg, pySeg, pzSeg, sqrt(pSeg*pSeg+mSeg*mSeg));
          mctSegment.SetProductionVertex(xSeg, ySeg, zSeg, 0);
          TEveTrackPropagator* prop = new TEveTrackPropagator("electron_brem_segment",
                                                              "electron brem segment",
                                                              magField.get(),
                                                              false);
          double propMaxR = rMax;
          double vSeg[3];
          vSeg[0] = seg.end.X();
          vSeg[1] = seg.end.Y();
          vSeg[2] = seg.end.Z();
          if ((std::sqrt(vSeg[0]*vSeg[0]+vSeg[1]*vSeg[1]) < rMax) && std::sqrt(xSeg*xSeg+ySeg*ySeg) < rMax) {
            if (xSeg*xSeg+ySeg*ySeg > vSeg[0]*vSeg[0]+vSeg[1]*vSeg[1]) {
              propMaxR = std::sqrt(xSeg*xSeg+ySeg*ySeg);
            }
            else {
              propMaxR = std::sqrt(vSeg[0]*vSeg[0]+vSeg[1]*vSeg[1]);
            }
          }
          prop->SetMaxR(propMaxR);
          prop->SetRnrReferences(false);
          prop->SetFitReferences(false);
          TEveTrack* segTrack = new TEveTrack(&mctSegment, prop);
          segTrack->SetAttLineAttMarker(particles);
          segTrack->SetElementTitle(Form("pSeg = %.3f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
                                      pSeg, acos(pzSeg / pSeg), atan2(pySeg, pxSeg),
                                      xSeg / cm, ySeg / cm, zSeg / cm));
          particles->AddElement(segTrack);
        }
      }
      else {
        // particles other than electrons
        // use per-track propagator to set maxR
        TEveTrackPropagator* prop = new TEveTrackPropagator(Form("propagator_%d", ip), "per-track propagator", magField.get(), false);
        prop->SetRnrReferences(false);
        // use this if we add track marks to electrons
        // prop->SetFitReferences(true);
        prop->SetFitReferences(false);
        double propMaxR = rMax;
        if ((std::sqrt(v[0]*v[0]+v[1]*v[1]) < rMax) && std::sqrt(x1*x1+y1*y1) < rMax) {
          if (x1*x1+y1*y1 > v[0]*v[0]+v[1]*v[1]) {
            propMaxR = std::sqrt(x1*x1+y1*y1);
          }
          else {
            propMaxR = std::sqrt(v[0]*v[0]+v[1]*v[1]);
          }
        }
        prop->SetMaxR(propMaxR);

        TEveTrack *track = new TEveTrack(&mct, prop);
        track->SetAttLineAttMarker(particles);
        track->SetElementTitle(Form("p = %.3f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
                                    p, acos(pz / p), atan2(py, px),
                                    x1 / cm, y1 / cm, z1 / cm));


        particles->AddElement(track);
      }
    /*
    // if the particle is a pi0, also draw the two photons, and set the endpoint
    // of the pi0 track
    if ((pdgID == 111) && (displayConfig.getBoolConfig("drawSimParticles")))
    {
      bool decayVtxSet = false;
      for (unsigned int i = 0; i < eventReader->simParticles.PDG->GetSize(); i++)
      {
        pdgID = (*eventReader->simParticles.PDG)[i];
        // keep only photons
        if (pdgID != 22)
          continue;
        // if ( (* eventReader->simParticles.generatorStatus)[i] != 1 ) continue;
        px = (*eventReader->simParticles.momentum_x)[i];
        py = (*eventReader->simParticles.momentum_y)[i];
        pz = (*eventReader->simParticles.momentum_z)[i];
        m = (*eventReader->simParticles.mass)[i];
        p = sqrt(px * px + py * py + pz * pz);
        float e = sqrt(p * p + m * m);
        // std::cout << "p = "<< p << std::endl;
        if (p < ParticleEnergyThreshold)
          continue;

        // std::cout << "PDG = "<< pdgID << std::endl;
        t = (*eventReader->simParticles.time)[i];
        x1 = (*eventReader->simParticles.vertex_x)[i] * mm;
        y1 = (*eventReader->simParticles.vertex_y)[i] * mm;
        z1 = (*eventReader->simParticles.vertex_z)[i] * mm;
        double r1 = sqrt(x1 * x1 + y1 * y1);
        // the two photons from a pi0 in the origin must come from small R
        // TODO IMPROVE FOR DDSIM FILES - WHERE PARTICLES PARENTS ARE STORED!
        if (r1 > 1.)
          continue;
        double sintheta = sqrt(px * px + py * py) / p;
        x2 = x1 + px / p * rMax / sintheta;
        y2 = y1 + py / p * rMax / sintheta;
        z2 = z1 + pz / p * rMax / sintheta;
        // double x2 = (* eventReader->simParticles.endpoint_x)[i] * mm;
        // double y2 = (* eventReader->simParticles.endpoint_y)[i] * mm;
        // double z2 = (* eventReader->simParticles.endpoint_z)[i] * mm;
        // double r2 = sqrt(x2*x2+y2*y2);
        // std::cout << "x1 y1 z1 x2 y2 z2 = "
        //      << x1 << " " << y1 << " " << z1 << " "
        //      << x2 << " " << y2 << " " << z2 << std::endl;
        // set pi0 decay point
        if (!decayVtxSet)
        {
          TEveVectorF v;
          v[0] = x1;
          v[1] = y1;
          v[2] = z1;
          TEvePathMark mark(TEvePathMark::kDecay, v);
          track->AddPathMark(mark);
          decayVtxSet = true;
        }
        TEveMCTrack mct;
        mct.SetPdgCode(pdgID);
        mct.SetMomentum(px, py, pz, sqrt(p * p + m * m));
        mct.SetProductionVertex(x1, y1, z1, t);
        TEveTrack *track = new TEveTrack(&mct, particles->GetPropagator());
        track->SetAttLineAttMarker(particles);
        track->SetElementTitle(Form("p = .3%f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
                                    p, acos(pz / p), atan2(py, px),
                                    x1 / cm, y1 / cm, z1 / cm));

        particles->AddElement(track);
      }
    }
    */
    }
    particles->MakeTracks();
  }
}

void EventDisplay::drawHits()
{
  createElementList(hits, "hits");

  // clear the hits of the previous event: the endcap collections are appended
  // to the barrel ones, so they would accumulate if the barrel is not drawn
  for (MyPointSet *hitSet : {vtxHits, mainTrkHits, siwrHits, ecalHits, hcalHits, muonHits})
    if (hitSet)
      hitSet->Reset();

  // hits in vertex
  if (displayConfig.getBoolConfig("drawVertexHits"))
  {
    std::cout << "Creating vertex barrel hits" << std::endl;
    drawHitCollection(vtxHits, "VTX",
                      eventReader->vertexBarrelHits,
                      -9999.);
    std::cout << "Creating vertex endcap hits" << std::endl;
    drawHitCollection(vtxHits, "VTX",
                      eventReader->vertexEndcapHits,
                      -9999., false, false, true);
  }

  // hits in the main tracker (DCH or STT)
  if (displayConfig.getBoolConfig("drawMainTrackerHits"))
  {
    std::cout << "Creating main tracker hits" << std::endl;
    drawHitCollection(mainTrkHits, mainTrackerName(displayConfig.getStringConfig("mainTrackerHits")),
                      eventReader->mainTrackerHits,
                      -9999.);
  }

  // hits in Si wrapper
  if (displayConfig.getBoolConfig("drawSiWrapperHits"))
  {
    std::cout << "Creating Si wrapper hits" << std::endl;
    drawHitCollection(siwrHits, "Si wrapper",
                      eventReader->siWrapperBarrelHits,
                      -9999.);
    drawHitCollection(siwrHits, "Si wrapper",
                      eventReader->siWrapperEndcapHits,
                      -9999., false, false, true);
  }

  // hits in ECAL
  if (displayConfig.getBoolConfig("drawECalBarrelHits"))
  {
    std::cout << "Creating ecal barrel hits" << std::endl;
    drawHitCollection(ecalHits, "ECAL",
                      eventReader->ecalBarrelHits,
                      HitEnergyThreshold);
  }

  if (displayConfig.getBoolConfig("drawECalEndcapHits"))
  {
    std::cout << "Creating ecal endcap hits" << std::endl;
    drawHitCollection(ecalHits, "ECAL",
                      eventReader->ecalEndcapHits,
                      HitEnergyThreshold, false, false, true);
  }

  // hits in HCAL
  if (displayConfig.getBoolConfig("drawHCalBarrelHits"))
  {
    std::cout << "Creating hcal barrel hits" << std::endl;
    drawHitCollection(hcalHits, "HCAL",
                      eventReader->hcalBarrelHits,
                      HitEnergyThreshold);
  }
  if (displayConfig.getBoolConfig("drawHCalEndcapHits"))
  {
    std::cout << "Creating hcal endcap hits" << std::endl;
    drawHitCollection(hcalHits, "HCAL",
                      eventReader->hcalEndcapHits,
                      HitEnergyThreshold, false, false, true);
  }

  // hits in muon detector
  if (displayConfig.getBoolConfig("drawMuonHits"))
  {
    std::cout << "Creating muon barrel hits" << std::endl;
    drawHitCollection(muonHits, "MUON",
                      eventReader->muonBarrelHits,
                      -9999.);
    std::cout << "Creating muon endcap hits" << std::endl;
    drawHitCollection(muonHits, "MUON",
                      eventReader->muonEndcapHits,
                      -9999., false, false, true);
  }
}

void EventDisplay::drawDigis()
{
  createElementList(digis, "digis");

  // clear the digis of the previous event: the endcap collections are appended
  // to the barrel ones, so they would accumulate if the barrel is not drawn
  for (MyPointSet *digiSet : {vtxDigis, mainTrkDigis, siwrDigis, ecalCells, hcalCells, muonCells, cells_merged})
    if (digiSet)
      digiSet->Reset();

  // digis in VTX
  if (displayConfig.getBoolConfig("drawVertexDigis"))
  {
    std::cout << "Creating vertex barrel digis" << std::endl;
    drawHitCollection(vtxDigis, "VTX",
                      eventReader->vertexBarrelDigis,
                      -9999., true);
    std::cout << "Creating vertex endcap digis" << std::endl;
    drawHitCollection(vtxDigis, "VTX",
                      eventReader->vertexEndcapDigis,
                      -9999., true, false, true);
  }

  // digis in the main tracker (DCH or STT)
  if (displayConfig.getBoolConfig("drawMainTrackerDigis"))
  {
    std::cout << "Creating main tracker digis" << std::endl;
    drawHitCollection(mainTrkDigis, mainTrackerName(displayConfig.getStringConfig("mainTrackerDigis")),
                      eventReader->mainTrackerDigis,
                      -9999., true);
  }

  // digis in Si wrapper
  if (displayConfig.getBoolConfig("drawSiWrapperDigis"))
  {
    std::cout << "Creating Si wrapper digis" << std::endl;
    drawHitCollection(siwrDigis, "Si wrapper",
                      eventReader->siWrapperBarrelDigis,
                      -9999., true);
    drawHitCollection(siwrDigis, "Si wrapper",
                      eventReader->siWrapperEndcapDigis,
                      -9999., true, false, true);
  }

  // digis (cells) in ECAL
  if (displayConfig.getBoolConfig("drawECalBarrelCells"))
  {
    std::cout << "Creating ecal barrel cells" << std::endl;
    drawHitCollection(ecalCells, "ECAL",
                      eventReader->ecalBarrelCells,
                      CellEnergyThreshold,
                      true, debug, false);
    // add tool tips
    if (debug) {
      for (unsigned int i = 0; i < eventReader->ecalBarrelCells.size(); i++) {
        float E = (*eventReader->ecalBarrelCells.energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        ULong_t cellID = (*eventReader->ecalBarrelCells.cellID)[i];
        ULong_t system = geomReader->SystemID(cellID);
        ULong_t layer = geomReader->ECalBarrelLayer(cellID);
        ULong_t module = geomReader->ECalBarrelModule(cellID);
        ULong_t theta = geomReader->ECalBarrelThetaBin(cellID);
        ecalCells->fTooltips.push_back(Form("Energy = %.3f GeV, System = %lu, Layer = %lu, Module = %lu, Theta = %lu", E, system, layer, module, theta));
      }
    }
  }

  if (displayConfig.getBoolConfig("drawECalEndcapCells"))
  {
    std::cout << "Creating ecal endcap cells" << std::endl;
    drawHitCollection(ecalCells, "ECAL",
                      eventReader->ecalEndcapCells,
                      CellEnergyThreshold,
                      true, debug, true);
    // add tool tips
    if (debug) {
      for (unsigned int i = 0; i < eventReader->ecalEndcapCells.size(); i++)
      {
        float E = (*eventReader->ecalEndcapCells.energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        ULong_t cellID = (*eventReader->ecalEndcapCells.cellID)[i];
        ULong_t system = geomReader->SystemID(cellID);
        Long_t side = geomReader->ECalEndCapSide(cellID);
        ULong_t wheel = geomReader->ECalEndCapWheel(cellID);
        ULong_t module = geomReader->ECalEndCapModule(cellID);
        ULong_t rho = geomReader->ECalEndCapRhoBin(cellID);
        ULong_t z = geomReader->ECalEndCapZBin(cellID);
        ULong_t layer = geomReader->ECalBarrelLayer(cellID);
        ecalCells->fTooltips.push_back(Form("Energy = %.3f GeV, System = %lu, Side = %ld, Wheel = %lu, Layer = %lu, Module = %lu, rho = %lu, z = %lu", E, system, side, wheel, layer, module, rho, z));
      }
    }
  }

  // digis (cells) in HCAL
  if (displayConfig.getBoolConfig("drawHCalBarrelCells"))
  {
    std::cout << "Creating hcal barrel cells" << std::endl;
    drawHitCollection(hcalCells, "HCAL",
                      eventReader->hcalBarrelCells,
                      CellEnergyThreshold,
                      true, debug, false);
    // add tool tips
    if (debug) {
      for (unsigned int i = 0; i < eventReader->hcalBarrelCells.size(); i++)
      {
        float E = (*eventReader->hcalBarrelCells.energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        ULong_t cellID = (*eventReader->hcalBarrelCells.cellID)[i];
        ULong_t system = geomReader->SystemID(cellID);
        ULong_t layer = geomReader->HCalBarrelLayer(cellID);
        ULong_t phi = geomReader->HCalBarrelPhiBin(cellID);
        ULong_t theta = geomReader->HCalBarrelThetaBin(cellID);
        hcalCells->fTooltips.push_back(Form("Energy = %.3f GeV, System = %lu, Layer = %lu, Phi = %lu, Theta = %lu", E, system, layer, phi, theta));
      }
    }
  }

  if (displayConfig.getBoolConfig("drawHCalEndcapCells"))
  {
    std::cout << "Creating hcal endcap cells" << std::endl;
    drawHitCollection(hcalCells, "HCAL",
                      eventReader->hcalEndcapCells,
                      CellEnergyThreshold,
                      true, debug, true);
    // add tool tips
    if (debug) {
      for (unsigned int i = 0; i < eventReader->hcalEndcapCells.size(); i++)
      {
        float E = (*eventReader->hcalEndcapCells.energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        ULong_t cellID = (*eventReader->hcalEndcapCells.cellID)[i];
        ULong_t system = geomReader->SystemID(cellID);
        ULong_t layer = geomReader->HCalEndCapLayer(cellID);
        ULong_t phi = geomReader->HCalEndCapPhiBin(cellID);
        ULong_t theta = geomReader->HCalEndCapThetaBin(cellID);
        hcalCells->fTooltips.push_back(Form("Energy = %.3f GeV, System = %lu, Layer = %lu, Phi = %lu, Theta = %lu", E, system, layer, phi, theta));
      }
    }
  }

  // digis in muon detector
  if (displayConfig.getBoolConfig("drawMuonCells"))
  {
    std::cout << "Creating muon barrel cells" << std::endl;
    drawHitCollection(muonCells, "MUON",
                      eventReader->muonBarrelCells,
                      CellEnergyThreshold,
                      true, false);

    std::cout << "Creating muon endcap cells" << std::endl;
    drawHitCollection(muonCells, "MUON",
                      eventReader->muonEndcapCells,
                      CellEnergyThreshold,
                      true, false, true);
  }

  // cells merged in ECAL barrel
  if (displayConfig.getBoolConfig("drawECalBarrelMergedCells"))
  {
    std::cout << "Creating merged ecal barrel cells" << std::endl;
    drawHitCollection(cells_merged, "ECAL merged",
                      eventReader->ecalBarrelMergedCells,
                      -9999.);
    // override default size and color
    cells_merged->SetMarkerSize(3);
    cells_merged->SetMarkerColor(kBlue);
  }
}

void EventDisplay::drawTracks()
{
  if (displayConfig.getBoolConfig("drawTracks"))
  {
    std::cout << "Creating tracks" << std::endl;
    if (tracks == nullptr)
    {
      tracks = new TEveTrackList("tracks");
      TEveTrackPropagator *trkProp = tracks->GetPropagator();
      trkProp->SetMagFieldObj(magField.get(), false);
      // trkProp->SetMaxR(rMax);
      // trkProp->SetMaxZ(geomReader->zMaxEndCap);
      // trkProp->SetMaxZ(10000.*mm);
      trkProp->SetRnrReferences(true);
      // the path marks of the track states in white (not to be confused with
      // the yellow digis of the tracker around them)
      trkProp->RefPMAtt().SetMarkerColor(kWhite);
      // trkProp->SetFitReferences(true);
      tracks->SetMainColor(kCyan);
      tracks->SetLineWidth(2.);
      tracks->SetLineStyle(5);
      gEve->AddElement(tracks);
    }
    else
      tracks->DestroyElements();

    // the transverse momentum of the tracks is obtained from their curvature
    // omega with the field used in the reconstruction: assume that it is the
    // one used to propagate the tracks (bField in the configuration), at the
    // origin (the reconstruction field is not stored in the event file)
    const float FCT = 2.99792458E-4f;
    float bField = std::fabs(magField->GetFieldD(0., 0., 0.).fZ);
    if (bField < 1e-6)
    {
      bField = geomReader->Bin;
      std::cout << "WARNING: the magnetic field is zero, the momentum of the tracks is computed from their curvature with B = "
                << bField << " T" << std::endl;
    }

    for (unsigned int ip = 0; ip < eventReader->tracks.subdetectorHitNumbers_begin->GetSize(); ip++)
    {
      if (debug) std::cout << "Track " << ip << std::endl;
      unsigned int trackStates_begin = (*eventReader->tracks.trackStates_begin)[ip];
      unsigned int trackStates_end = (*eventReader->tracks.trackStates_end)[ip];
      if (debug) std::cout << "  trackStates begin , end = " << trackStates_begin << " , " << trackStates_end << std::endl;
      // store in vectors the track state positions and momenta
      float x[5], y[5], z[5], omega[5], tanLambda[5], phi[5], d0[5], z0[5];
      for (unsigned int i=0; i<5; i++) x[i]=-1e9; // set to some large value to check later if track state has been filled
      for (unsigned int its = trackStates_begin; its < trackStates_end; its++) {
        int location = (*eventReader->tracks.trackState_location)[its];
        if (debug) {
          std::cout << "Trackstate " << its << std::endl;
          std::cout << "  location = " << location << std::endl;
        }
        if (location < 0 || location > 4) continue;
        x[location] = (*eventReader->tracks.trackState_referencePoint_x)[its];
        y[location] = (*eventReader->tracks.trackState_referencePoint_y)[its];
        z[location] = (*eventReader->tracks.trackState_referencePoint_z)[its];
        omega[location] = (*eventReader->tracks.trackState_omega)[its];
        tanLambda[location] = (*eventReader->tracks.trackState_tanLambda)[its];
        phi[location] = (*eventReader->tracks.trackState_phi)[its];
        d0[location] = (*eventReader->tracks.trackState_d0)[its];
        z0[location] = (*eventReader->tracks.trackState_z0)[its];
        // the position of the track state: the point of closest approach to
        // the reference point, displaced from it by d0 and z0
        x[location] -= d0[location] * std::sin(phi[location]);
        y[location] += d0[location] * std::cos(phi[location]);
        z[location] += z0[location];
        if (debug) {
          std::cout << "  x = " << x[location] << std::endl;
          std::cout << "  y = " << y[location] << std::endl;
          std::cout << "  z = " << z[location] << std::endl;
          std::cout << "  r = " << sqrt(x[location]*x[location] + y[location]*y[location]) << std::endl;
          std::cout << "  omega = " << omega[location] << std::endl;
          std::cout << "  phi = " << phi[location] << std::endl;
          std::cout << "  tanLambda = " << tanLambda[location] << std::endl;
          std::cout << "  d0 = " << d0[location] << std::endl;
          std::cout << "  z0 = " << z0[location] << std::endl;
          float radius = 1.f / std::fabs(omega[location]);
          float pxy = FCT * bField * radius;
          float pz = tanLambda[location] * pxy;
          float p = std::sqrt(pxy*pxy+pz*pz);
          std::cout << "  p = " << p << std::endl;
        }
      }
      const int nSubDetectors(5);
      int nhits[nSubDetectors];
      for (int ihit = 0 ; ihit < nSubDetectors; ++ihit) nhits[ihit]=0;
      unsigned int hits_begin = (*eventReader->tracks.subdetectorHitNumbers_begin)[ip];
      unsigned int hits_end = (*eventReader->tracks.subdetectorHitNumbers_end)[ip];
      if (debug) std::cout << "  subdetectorHitNumbers_begin , end = " << hits_begin << " , " << hits_end << std::endl;
      // (only the first nSubDetectors numbers are shown)
      for (unsigned int ih = hits_begin; ih < hits_end && ih - hits_begin < nSubDetectors; ih++) {
        nhits[ih-hits_begin] = (*eventReader->tracks.subdetectorHitNumbers)[ih];
      }

      // const int tsOrig=2; // which track state to use for track origin: 0=at other, 1=at IP, 2=at first hit, 3=at last hit, 4=at ECAL
      int tsOrig=1; // which track state to use for track origin: 0=at other, 1=at IP, 2=at first hit, 3=at last hit, 4=at ECAL
      // if the track has no state there, use the first one found along the
      // track (at first hit, at last hit, at ECAL, at other)
      if (x[tsOrig] < -1e8) {
        tsOrig = -1;
        for (int i : {2, 3, 4, 0}) {
          if (x[i] > -1e8) {
            tsOrig = i;
            break;
          }
        }
        if (tsOrig < 0)
          std::cout << "WARNING: track " << ip << " has no valid track state, it is not drawn" << std::endl;
        else if (debug)
          std::cout << "  no track state at IP, origin at location " << tsOrig << std::endl;
      }
      if (tsOrig >= 0) {
        TEveRecTrack t;
        float x1 = x[tsOrig]*mm;
        float y1 = y[tsOrig]*mm;
        float z1 = z[tsOrig]*mm;
        int charge = MathUtils::sign(omega[tsOrig]);
        float radius = 1.f / std::fabs(omega[tsOrig]);
        float pxy = FCT * bField * radius;
        float px = pxy * std::cos(phi[tsOrig]);
        float py = pxy * std::sin(phi[tsOrig]);
        float pz = tanLambda[tsOrig] * pxy;
        float p = std::sqrt(px*px+py*py+pz*pz);
        t.fV.Set(x1, y1, z1);
        t.fP.Set(px, py, pz);
        t.fSign = charge;
        TEveTrack* track = new TEveTrack(&t, tracks->GetPropagator());
        track->SetElementName(Form("Track %d", ip));
        track->SetAttLineAttMarker(tracks);
        track->SetElementTitle(Form("p = %.3f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm\nhits = %d/%d/%d/%d/%d",
          p, std::acos(pz / p), std::atan2(py, px),
          x1 / cm, y1 / cm, z1 / cm,
          nhits[0], nhits[1], nhits[2], nhits[3], nhits[4]));

        // add other track states (at last hit, at ECAL, and at other=2nd ECAL projection) references
        track->SetRnrPoints(true);
        track->SetMarkerStyle(4);
        for (int i=0; i<5; i++) {
          if (i==tsOrig) {
              if (debug) std::cout << "Skipping track state at location " << i << std::endl;
              // continue;
          }
          if (x[i]<-1e8) {
            if (debug) std::cout << "Skipping track state at location " << i << std::endl;
            continue;
          }
          float radiusState = 1.f / std::fabs(omega[i]);
          float pxyState = FCT * bField * radiusState;
          float pxState = pxyState * std::cos(phi[i]);
          float pyState = pxyState * std::sin(phi[i]);
          float pzState = tanLambda[i] * pxyState;
          track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kReference,
                                           TEveVectorD(x[i]*mm, y[i]*mm, z[i]*mm),
                                           TEveVectorD(pxState, pyState, pzState)));
        }
        if (debug) track->PrintPathMarks();
        // could also save decay ...
        tracks->AddElement(track);
      }
    }
    tracks->MakeTracks();
  }
}

void EventDisplay::drawClusters(ClusterSet &set)
{
  const std::string &clusterType = set.type;
  // the clusters and their cells in the current event
  const EventReader::Clusters &clusters = *set.clusters;
  const EventReader::ClusterCells &cells = clusters.cells;
  // read once, not in the loop over the clusters
  const bool drawDirection = displayConfig.getBoolConfig("drawClusterDirection");

  std::cout << "Drawing " << set.label << " clusters" << std::endl;
  std::vector<std::unique_ptr<CaloCluster>> *clusterData = &set.data;

  const double rMin = geomReader->rMin;
  const double zMinEC = geomReader->zMinEndCap;
  const double zMaxEC = geomReader->zMaxEndCap;
  const double alpha = geomReader->alpha;
  const double thetaGrid = geomReader->thetaGrid;
  const double gridPhi = geomReader->gridPhi;
  const double phiMin = geomReader->phiMin;
  const double thetaGridHCal = geomReader->thetaGridHCal;
  const double gridPhiHCal = geomReader->gridPhiHCal;
  const double phiMinHCal = geomReader->phiMinHCal;

  // create the graphic containers, or empty them
  // clusters in 3D
  if (debug)
  {
    std::cout << "  Creating the containers to hold the 3D clusters" << std::endl;
  }
  createElementList(clustersFolder, "clusters");
  if (set.list3D == nullptr)
  {
    set.list3D = new TEveElementList(Form("%s clusters (E>%.2f GeV)",
                                          clusterType.c_str(),
                                          ClusterEnergyThreshold));
    clustersFolder->AddElement(set.list3D);
  }
  else
  {
    set.list3D->DestroyElements();
    set.center = nullptr;
    set.direction = nullptr;
  }

  // clusters in 2D
  if (debug)
  {
    std::cout << "  Creating the containers to hold the 2D projections of the clusters" << std::endl;
  }
  if (set.listRhoZ == nullptr)
  {
    set.listRhoZ = new TEveElementList(Form("%s clusters in rho-z (E>%.2f GeV)",
                                            clusterType.c_str(),
                                            ClusterEnergyThreshold));
    // kept between the events, when the projected replica it is added to is deleted
    set.listRhoZ->IncDenyDestroy();
  }
  else
    set.listRhoZ->DestroyElements();

  if (set.listRhoPhi == nullptr)
  {
    set.listRhoPhi = new TEveElementList(Form("%s clusters in rho-phi (E>%.2f GeV)",
                                              clusterType.c_str(),
                                              ClusterEnergyThreshold));
    set.listRhoPhi->IncDenyDestroy();
  }
  else
    set.listRhoPhi->DestroyElements();

  TEveElementList *clusters_3D = set.list3D;
  TEveElementList *clusters_rhoz = set.listRhoZ;
  TEveElementList *clusters_rhophi = set.listRhoPhi;

  if (debug)
  {
    std::cout << "  Filling the containers" << std::endl;
  }
  unsigned int nClusters = clusters.size();
  if (debug)
  {
    std::cout << "  n(clusters) = " << nClusters << std::endl;
    std::cout << "  Looping over clusters to retrieve barycenter" << std::endl;
  }

  // the barycenters and directions of the clusters
  if (set.center == nullptr)
  {
    set.center = new TEvePointSet();
    set.center->SetName("cluster barycenters");
    set.center->SetMarkerColor(set.centerColor);
    set.list3D->AddElement(set.center);
  }
  else
    set.center->Reset();

  if (drawDirection) {
    if (set.direction == nullptr) {
      set.direction = new TEveStraightLineSet("cluster directions");
      set.list3D->AddElement(set.direction);
    }
    else
      set.direction->DestroyElements();
  }
  TEvePointSet *clustersCenter = set.center;
  TEveStraightLineSet *clustersDirection = set.direction;

  if (drawDirection) {
    clustersCenter->SetMarkerStyle(4);
    clustersCenter->SetMarkerSize(5);
    clustersDirection->SetLineColor(kWhite);
    clustersDirection->SetLineWidth(5);
  }

  for (unsigned int i = 0; i < nClusters; i++) {
    float E = (*clusters.energy)[i];
    if (E < ClusterEnergyThreshold)
      continue;

    float x0(0.), y0(0.), z0(0.);
    float theta(0.), phi(0.);
    // topo cluster positions are in cm while calo clusters and hits/cells in mm ... ?
    x0 = (*clusters.x)[i] * mm;
    y0 = (*clusters.y)[i] * mm;
    z0 = (*clusters.z)[i] * mm;
    if (drawDirection) {
      theta = (*clusters.theta)[i];
      phi = (*clusters.phi)[i];
    }

    clustersCenter->SetNextPoint(x0, y0, z0);
    if (theta!=0.) {
      float t1(0), t2(0);
      if (geomReader->isInECalBarrel(std::sqrt(x0*x0 + y0*y0), z0)) {
        // for cluster in barrel
        float a = std::sin(theta)*std::sin(theta);
        float b = 2*std::sin(theta)*(x0*std::cos(phi) + y0*std::sin(phi));
        float c = x0*x0 + y0*y0 - rMin*rMin;
        t1 = MathUtils::solveQuadratic(a,b,c,1);
        double rMaxCalo = doHCal ? geomReader->rMaxHCal : geomReader->rMax;
        c = x0*x0 + y0*y0 - rMaxCalo*rMaxCalo;
        t2 = MathUtils::solveQuadratic(a,b,c,1);
      }
      else if (geomReader->isInECalEndCap(std::sqrt(x0*x0 + y0*y0), z0)) {
        // for cluster in endcap
        double zMin = zMinEC;
        double zMax = zMaxEC;
        if (z0<0) {
          zMin = -zMin;
          zMax = -zMax;
        }
        t1 = (zMin-z0)/std::cos(theta);
        t2 = (zMax-z0)/std::cos(theta);
      }
      else {
        //
        std::cout << "Cluster barycentre is not in ECal barrel nor endcap!" << std::endl;
        std::cout << "  x0 = " << x0 << ", y0 = " << y0 << ", z0 = " << z0 << std::endl;
      }
      if (t1!=0.)
        clustersDirection->AddLine(x0 + t1*std::sin(theta)*std::cos(phi),
                                   y0 + t1*std::sin(theta)*std::sin(phi),
                                   z0 + t1*std::cos(theta),
                                   x0 + t2*std::sin(theta)*std::cos(phi),
                                   y0 + t2*std::sin(theta)*std::sin(phi),
                                   z0 + t2*std::cos(theta));
    }
  }

  if (debug)
  {
    std::cout << "  Creating the cluster graphic objects" << std::endl;
  }

  // now create the visual representation of the clustered cells
  std::vector<TEveQuadSet *> qs_rhoz;
  std::vector<TEveQuadSet *> qs_rhophi;
  std::vector<TEveBoxSet *> bs;

  // now loop over the clusters and fill the containers
  // in 3D this is done in this first loop
  // for 2D projection, in this loop integrate energy over projected dimension
  // and fill map of 2D cell ID vs energy
  // then the filling is done in a second loop later
  if (debug)
  {
    std::cout << "  Looping over clusters to fill 3D projections" << std::endl;
  }
  int iclAboveThreshold = 0;
  for (unsigned int i = 0; i < nClusters; i++)
  {
    float energy = (*clusters.energy)[i];
    float xcl = (*clusters.x)[i];
    float ycl = (*clusters.y)[i];
    float zcl = (*clusters.z)[i];
    float rcl = sqrt(xcl * xcl + ycl * ycl);
    float phicl = atan2(ycl, xcl);
    float thetacl = atan2(rcl, zcl);

    if (energy < ClusterEnergyThreshold)
    {
      qs_rhoz.push_back(nullptr);
      qs_rhophi.push_back(nullptr);
      bs.push_back(nullptr);
    }
    else
    {
      TEveQuadSet *aqs = new TEveQuadSet(TEveQuadSet::kQT_FreeQuad, false, 32,
                                         Form("%s cluster %d", clusterType.c_str(), (int)i));
      aqs->SetMainTransparency(80);
      // by calling SetOwnIds(kTRUE) the digit-set becomes
      // the owner of the assigned objects and deletes
      // them on destruction.
      aqs->SetOwnIds(kTRUE);
      aqs->SetPalette(cellPaletteRhoZ);
      aqs->SetTitle(Form("E = %f GeV\nR = %f mm\ntheta = %f\nphi = %f",
                         energy,
                         rcl,
                         thetacl,
                         phicl));
      qs_rhoz.push_back(aqs);
      clusters_rhoz->AddElement(aqs);

      TEveQuadSet *aqs2 = new TEveQuadSet(TEveQuadSet::kQT_FreeQuad, false, 32,
                                          Form("%s cluster %d", clusterType.c_str(), (int)i));
      aqs2->SetMainTransparency(80);
      aqs2->SetOwnIds(kTRUE);
      aqs2->SetPalette(cellPaletteRhoPhi);
      aqs2->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",
                          energy,
                          rcl,
                          thetacl,
                          phicl));
      qs_rhophi.push_back(aqs2);
      clusters_rhophi->AddElement(aqs2);
      //}

      // window of the sliding-window clusters
      if (set.drawWindow)
        drawSlidingWindow(set.window, xcl, ycl, zcl, Form("%s cluster %d window", clusterType.c_str(), (int)i),
                          clusters_rhoz, clusters_rhophi);

      TEveBoxSet *_bs = new TEveBoxSet(Form("%s cluster %d", clusterType.c_str(), (int)i), "");
      _bs->Reset(TEveBoxSet::kBT_FreeBox, false, 32);
      _bs->SetMainTransparency(80);
      _bs->SetOwnIds(kTRUE);
      _bs->SetPalette(cellPalette3D);
      _bs->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",
                         energy,
                         rcl,
                         thetacl,
                         phicl));
      bs.push_back(_bs);
      clusters_3D->AddElement(_bs);

      // try showing the barycenter vs layer
      if (drawClustersBarycenterVsLayer)
      {
        TEvePointSet *layerBarycenters = new TEvePointSet();
        layerBarycenters->SetName(Form("%s cluster %d barycenters", clusterType.c_str(), (int)i));
        layerBarycenters->SetMarkerColor(kGreen);
        layerBarycenters->SetMarkerStyle(4);
        layerBarycenters->SetMarkerSize(3);
        clusters_3D->AddElement(layerBarycenters);
        CaloCluster *clus = (*clusterData)[iclAboveThreshold].get();
        // (not for the muon tagger)
        for (CaloCluster::Part part : {CaloCluster::kECalBarrel, CaloCluster::kECalEndcap,
                                       CaloCluster::kHCalBarrel, CaloCluster::kHCalEndcap})
        {
          for (unsigned int iLayer = 0; iLayer < clus->getNLayers(part); iLayer++)
          {
            if (clus->getEnergyInLayer(part, iLayer) > 0)
            {
              TVector3 barycenter = clus->getBarycenterInLayer(part, iLayer);
              layerBarycenters->SetNextPoint(barycenter.X() / cm, barycenter.Y() / cm, barycenter.Z() / cm);
            }
          }
        }
      }
      iclAboveThreshold++;
    }
  }
  // loop over cells
  std::map<ProjectedCellKey, double> cellEnergies_rhoz;
  std::map<ProjectedCellKey, double> cellEnergies_rhophi;
  // if the cells of the clusters are not available, the clusters are
  // drawn without their cells (barycenter and direction only)
  unsigned int nCells = cells.size();

  // index of the cluster each cell belongs to (-1 if none), built once:
  // the cells of cluster j are [hits_begin[j], hits_end[j])
  std::vector<int> cellCluster(nCells, -1);
  for (unsigned int j = 0; j < nClusters; j++)
  {
    unsigned int first_hit = (*clusters.hits_begin)[j];
    unsigned int last_hit = (*clusters.hits_end)[j];
    // in case of overlapping ranges, keep the first cluster
    for (unsigned int i = first_hit; i < std::min(last_hit, nCells); i++)
      if (cellCluster[i] == -1)
        cellCluster[i] = j;
  }

  // check that the cells of the sliding-window clusters are inside their window
  if (debug && set.drawWindow && set.window.nTheta > 0)
  {
    std::vector<int> nOutside(nClusters, 0);
    for (unsigned int i = 0; i < nCells; i++)
    {
      int icl = cellCluster[i];
      if (icl == -1)
        continue;
      double thetaLow, thetaHigh, phiLow, phiHigh;
      getSlidingWindowRange(set.window, (*clusters.x)[icl],
                            (*clusters.y)[icl],
                            (*clusters.z)[icl],
                            thetaLow, thetaHigh, phiLow, phiHigh);
      const double theta = std::atan2(std::sqrt(cells.x[i] * cells.x[i] + cells.y[i] * cells.y[i]), cells.z[i]);
      // phi of the cell in the range of the window
      double phi = std::atan2(cells.y[i], cells.x[i]);
      while (phi < phiLow)
        phi += 2 * TMath::Pi();
      while (phi > phiHigh)
        phi -= 2 * TMath::Pi();
      if (theta < thetaLow || theta > thetaHigh || phi < phiLow)
        nOutside[icl]++;
    }
    for (unsigned int j = 0; j < nClusters; j++)
      std::cout << "  sw cluster " << j << ": " << nOutside[j] << " cells outside the window" << std::endl;
  }

  for (unsigned int i = 0; i < nCells; i++)
  {
    int icl = cellCluster[i];
    if (icl == -1)
      continue; // should never happen..
    ULong_t cellID;
    float energy, x_center, y_center, z_center;
    cellID = cells.cellID[i];
    if (cellID == 0)
      continue; // cell not available (see EventReader::fillClusterCells)
    energy = cells.energy[i];
    x_center = cells.x[i] * mm;
    y_center = cells.y[i] * mm;
    z_center = cells.z[i] * mm;
    float r_center = sqrt(x_center * x_center + y_center * y_center);
    int systemID = (int)DetectorGeometry::SystemID(cellID);
    // for 2D projections, identify the projected cell by
    // ECAL: layer-theta for rho-z, layer-module for rho-phi
    // HCAL: layer-theta for rho-z, layer-phi for rho-phi
    // ECAL endcap: z bin (the layer) and rho bin for rho-z, not drawn in rho-phi and 3D
    // HCAL endcap: layer-theta for rho-z, not drawn in rho-phi
    // (in rho-z also the sign of y, as the cells above and below the beam
    // axis are drawn separately)
    float r_in(0.), r_out(0.);
    int layer(0), thetaID(0), moduleID(0), phiID(0);
    ProjectedCellKey rhophiID{}, rhozID{};
    bool inRhoPhi = true, in3D = true;

    if (systemID == DetectorGeometry::kECalBarrel)
    {
      layer = (int)DetectorGeometry::ECalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::ECalBarrelThetaBin(cellID);
      moduleID = (int)DetectorGeometry::ECalBarrelModule(cellID);
      if (layer < 0 || layer >= geomReader->nLayers)
      {
        std::cout << "WARNING: layer " << layer << " of system " << systemID << " out of range [0, " << geomReader->nLayers << "), skipping cell" << std::endl;
        continue;
      }
      rhophiID = {systemID, icl, layer, moduleID};
      rhozID = {systemID, icl, layer, thetaID * 2 + (y_center > 0)};
      r_in = geomReader->r[layer];
      r_out = geomReader->r[layer + 1];
    }
    else if (systemID == DetectorGeometry::kHCalBarrel)
    {
      layer = (int)DetectorGeometry::HCalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::HCalBarrelThetaBin(cellID);
      phiID = (int)DetectorGeometry::HCalBarrelPhiBin(cellID);
      if (layer < 0 || layer >= geomReader->nLayersHCal)
      {
        std::cout << "WARNING: layer " << layer << " of system " << systemID << " out of range [0, " << geomReader->nLayersHCal << "), skipping cell" << std::endl;
        continue;
      }
      rhophiID = {systemID, icl, layer, phiID};
      rhozID = {systemID, icl, layer, thetaID * 2 + (y_center > 0)};
      r_in = geomReader->rHCal[layer];
      r_out = geomReader->rHCal[layer + 1];
    }
    else if (systemID == DetectorGeometry::kECalEndcap)
    {
      layer = (int)DetectorGeometry::ECalEndCapZBin(cellID);
      rhozID = {systemID, icl, layer, ecalEndcapRhoZBin(cellID, y_center)};
      inRhoPhi = in3D = false;
    }
    else if (systemID == DetectorGeometry::kHCalEndcap)
    {
      layer = (int)DetectorGeometry::HCalEndCapLayer(cellID);
      thetaID = (int)DetectorGeometry::HCalEndCapThetaBin(cellID);
      rhozID = {systemID, icl, layer, thetaID * 2 + (y_center > 0)};
      inRhoPhi = false;
    }
    else if (systemID == DetectorGeometry::kMuonBarrel || systemID == DetectorGeometry::kMuonEndcap)
    {
      // cells of MUON barrel and endcap are not drawn
      continue;
    }
    else
    {
      std::cout << "WARNING: unknown system " << systemID << " for cell " << cellID << ", skipping cell" << std::endl;
      continue;
    }

    float theta_center = atan2(r_center, z_center);
    // float phi_center = atan2(y_center, x_center);

    // cluster cells in rho-z projection
    // skip clusters with energy below threshold
    if (qs_rhoz[icl] != nullptr)
    {
      // if cell ID (in rho-z) already inserted in map, sum energy of this cell
      if (cellEnergies_rhoz.count(rhozID))
      {
        // debug
        // std::cout << "Cell with layer, theta = " << layer << " " << thetaID << " already in map with energy " << cellEnergies_rhoz[rhozID] << " , adding energy " << energy << std::endl;
        cellEnergies_rhoz[rhozID] += energy;
      }
      // otherwise insert new cell in map
      else
      {
        // debug
        // std::cout << "Cell with layer, theta = " << layer << " " << thetaID << " NOT in map, inserting energy " << energy << std::endl;
        cellEnergies_rhoz[rhozID] = energy;
      }
      // if cell ID (in rho-phi) already inserted in map, sum energy of this cell
      if (!inRhoPhi)
      {
        // not drawn in rho-phi
      }
      else if (cellEnergies_rhophi.count(rhophiID))
      {
        // debug
        // std::cout << "Cell with layer, module = " << layer << " " << moduleID << " already in map with energy " << cellEnergies_rhophi[rhophiID] << " , adding energy " << energy << std::endl;
        cellEnergies_rhophi[rhophiID] += energy;
      }
      // otherwise insert new cell in map
      else
      {
        // std::cout << "Cell with layer, module = " << layer << " " << moduleID << " NOT in map, inserting energy " << energy << std::endl;
        cellEnergies_rhophi[rhophiID] = energy;
      }
    }

    // cluster cells in 3D
    if (bs[icl] != nullptr && in3D)
    {
      float verts3D[24];
      if (systemID == DetectorGeometry::kECalBarrel)
      {
        double Lin = geomReader->getL(alpha, rMin, r_in);
        double Lout = geomReader->getL(alpha, rMin, r_out);
        double deltaL = rMin * sin(gridPhi / 2.0) * sin(alpha);
        for (int j = 0; j < geomReader->mergedModules[layer]; j++)
        {
          int iModule = moduleID + j;
          double phi0 = phiMin + iModule * gridPhi;

          verts3D[0] = rMin * cos(phi0) + (Lin + deltaL) * cos(phi0 + alpha);
          verts3D[1] = rMin * sin(phi0) + (Lin + deltaL) * sin(phi0 + alpha);
          verts3D[2] = r_in / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);

          verts3D[3] = rMin * cos(phi0 + gridPhi) + (Lin - deltaL) * cos(phi0 + gridPhi + alpha);
          verts3D[4] = rMin * sin(phi0 + gridPhi) + (Lin - deltaL) * sin(phi0 + gridPhi + alpha);
          verts3D[5] = r_in / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);

          verts3D[6] = rMin * cos(phi0 + gridPhi) + (Lin - deltaL) * cos(phi0 + gridPhi + alpha);
          verts3D[7] = rMin * sin(phi0 + gridPhi) + (Lin - deltaL) * sin(phi0 + gridPhi + alpha);
          verts3D[8] = r_in / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);

          verts3D[9] = rMin * cos(phi0) + (Lin + deltaL) * cos(phi0 + alpha);
          verts3D[10] = rMin * sin(phi0) + (Lin + deltaL) * sin(phi0 + alpha);
          verts3D[11] = r_in / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);

          verts3D[12] = rMin * cos(phi0) + (Lout + deltaL) * cos(phi0 + alpha);
          verts3D[13] = rMin * sin(phi0) + (Lout + deltaL) * sin(phi0 + alpha);
          verts3D[14] = r_out / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);

          verts3D[15] = rMin * cos(phi0 + gridPhi) + (Lout - deltaL) * cos(phi0 + gridPhi + alpha);
          verts3D[16] = rMin * sin(phi0 + gridPhi) + (Lout - deltaL) * sin(phi0 + gridPhi + alpha);
          verts3D[17] = r_out / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);

          verts3D[18] = rMin * cos(phi0 + gridPhi) + (Lout - deltaL) * cos(phi0 + gridPhi + alpha);
          verts3D[19] = rMin * sin(phi0 + gridPhi) + (Lout - deltaL) * sin(phi0 + gridPhi + alpha);
          verts3D[20] = r_out / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);

          verts3D[21] = rMin * cos(phi0) + (Lout + deltaL) * cos(phi0 + alpha);
          verts3D[22] = rMin * sin(phi0) + (Lout + deltaL) * sin(phi0 + alpha);
          verts3D[23] = r_out / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);

          bs[icl]->AddBox(verts3D);
          bs[icl]->DigitValue((int)(1000 * energy));
          // bs->BoxId(new TNamed(Form("Cell %lu", cellID), "Dong!"));
        }
      }
      else if (systemID == DetectorGeometry::kHCalBarrel)
      {
        // HCAL
        float phiLow = phiMinHCal + gridPhiHCal * phiID;
        float phiHigh = phiMinHCal + gridPhiHCal * (phiID + 1);
        float thetaMin = theta_center - thetaGridHCal / 2.;
        float thetaMax = theta_center + thetaGridHCal / 2.;

        verts3D[0] = r_in * cos(phiLow);
        verts3D[1] = r_in * sin(phiLow);
        verts3D[2] = r_in / tan(thetaMin);

        verts3D[3] = r_in * cos(phiHigh);
        verts3D[4] = r_in * sin(phiHigh);
        verts3D[5] = r_in / tan(thetaMin);

        verts3D[6] = r_in * cos(phiHigh);
        verts3D[7] = r_in * sin(phiHigh);
        verts3D[8] = r_in / tan(thetaMax);

        verts3D[9] = r_in * cos(phiLow);
        verts3D[10] = r_in * sin(phiLow);
        verts3D[11] = r_in / tan(thetaMax);

        verts3D[12] = r_out * cos(phiLow);
        verts3D[13] = r_out * sin(phiLow);
        verts3D[14] = r_out / tan(thetaMin);

        verts3D[15] = r_out * cos(phiHigh);
        verts3D[16] = r_out * sin(phiHigh);
        verts3D[17] = r_out / tan(thetaMin);

        verts3D[18] = r_out * cos(phiHigh);
        verts3D[19] = r_out * sin(phiHigh);
        verts3D[20] = r_out / tan(thetaMax);

        verts3D[21] = r_out * cos(phiLow);
        verts3D[22] = r_out * sin(phiLow);
        verts3D[23] = r_out / tan(thetaMax);

        bs[icl]->AddBox(verts3D);
        bs[icl]->DigitValue((int)(1000 * energy));
        // bs->BoxId(new TNamed(Form("Cell %lu", cellID), "Dong!"));
      }
      else if (systemID == DetectorGeometry::kHCalEndcap)
      {
        // a sector of ring between two radii and two z planes, as the HCAL
        // barrel cells with the z of thetaMin and thetaMax
        double rhoMin, rhoMax, zLow, zHigh, phiLow, phiHigh;
        if (geomReader->HCalEndCapCellBounds(cellID, rhoMin, rhoMax, zLow, zHigh, phiLow, phiHigh))
        {
          const double rs[2] = {rhoMin, rhoMax};
          for (int j = 0; j < 2; j++)
          {
            const double r = rs[j];
            float *v = verts3D + 12 * j;
            v[0] = r * cos(phiLow);  v[1] = r * sin(phiLow);  v[2] = zHigh;
            v[3] = r * cos(phiHigh); v[4] = r * sin(phiHigh); v[5] = zHigh;
            v[6] = r * cos(phiHigh); v[7] = r * sin(phiHigh); v[8] = zLow;
            v[9] = r * cos(phiLow);  v[10] = r * sin(phiLow); v[11] = zLow;
          }
          bs[icl]->AddBox(verts3D);
          bs[icl]->DigitValue((int)(1000 * energy));
        }
      }
    }
  }

  // and now draw the 2D projections
  float verts[12];
  // cluster cells in rho-z projection
  // for this, unfortunately, I have to loop over all cells
  // because I need to retrieve their position
  // I cannot just loop over the map of the energies :(
  // in alternative I should add methods to determine the
  // center of the cells from theta bin, layer, and alpha...

  if (debug)
    std::cout << "  Draw rho-z view of clustered cells" << std::endl;
  for (unsigned int i = 0; i < nCells; i++)
  {
    int icl = cellCluster[i];
    if (icl == -1)
      continue; // should never happen..
    ULong_t cellID;
    float energy;
    cellID = cells.cellID[i];

    int system = (int)DetectorGeometry::SystemID(cellID);
    int layer(0), thetaID(0);
    ProjectedCellKey rhozID{};
    float r_in(0.), r_out(0.), thetaMin(0.), thetaMax(0.);
    if (system == DetectorGeometry::kECalBarrel)
    {
      layer = (int)DetectorGeometry::ECalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::ECalBarrelThetaBin(cellID);
      rhozID = {system, icl, layer, thetaID * 2 + (cells.y[i] > 0)};
      r_in = geomReader->r[layer];
      r_out = geomReader->r[layer + 1];
    }
    else if (system == DetectorGeometry::kHCalBarrel)
    {
      layer = (int)DetectorGeometry::HCalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::HCalBarrelThetaBin(cellID);
      rhozID = {system, icl, layer, thetaID * 2 + (cells.y[i] > 0)};
      r_in = geomReader->rHCal[layer];
      r_out = geomReader->rHCal[layer + 1];
    }
    else if (system == DetectorGeometry::kECalEndcap)
    {
      layer = (int)DetectorGeometry::ECalEndCapZBin(cellID);
      rhozID = {system, icl, layer, ecalEndcapRhoZBin(cellID, cells.y[i])};
    }
    else if (system == DetectorGeometry::kHCalEndcap)
    {
      layer = (int)DetectorGeometry::HCalEndCapLayer(cellID);
      thetaID = (int)DetectorGeometry::HCalEndCapThetaBin(cellID);
      rhozID = {system, icl, layer, thetaID * 2 + (cells.y[i] > 0)};
    }
    else
    {
      // cells of other systems are not drawn
      continue;
    }
    if (cellEnergies_rhoz.count(rhozID) == 0)
      continue;
    energy = cellEnergies_rhoz[rhozID];
    cellEnergies_rhoz.erase(rhozID);
    float x_center, y_center, z_center;
    x_center = cells.x[i] * mm;
    y_center = cells.y[i] * mm;
    z_center = cells.z[i] * mm;
    float r_center = sqrt(x_center * x_center + y_center * y_center);
    if (system == DetectorGeometry::kECalEndcap)
    {
      // the cell is a rectangle in rho-z, from the segmentation (the position
      // of the cell in the file is displaced in phi, along the blades, so its
      // rho is not the centre of the cell in rho)
      double rho, z, dRho, dZ;
      if (!geomReader->ECalEndCapCellRhoZ(cellID, rho, z, dRho, dZ))
      {
        std::cout << "WARNING: ECAL endcap cell " << cellID << " outside of the segmentation, skipping cell" << std::endl;
        continue;
      }
      // (warn once if the position is outside of the cell, e.g. for a file
      // produced with a different version of the geometry)
      static bool warnedPosition = false;
      if (debug && !warnedPosition && (std::fabs(rho - r_center) > dRho / 2. || std::fabs(z - z_center) > dZ / 2.))
      {
        std::cout << "  ECAL endcap cell " << cellID << " at rho, z = " << r_center / mm << ", " << z_center / mm
                  << " mm, outside of the cell from the segmentation (rho, z = " << rho / mm << " +- " << dRho / 2. / mm
                  << ", " << z / mm << " +- " << dZ / 2. / mm << " mm): different version of the geometry?" << std::endl;
        warnedPosition = true;
      }
      float sign = MathUtils::sign(y_center);
      verts[0] = z - dZ / 2.;
      verts[1] = (rho - dRho / 2.) * sign;
      verts[3] = z + dZ / 2.;
      verts[4] = (rho - dRho / 2.) * sign;
      verts[6] = z + dZ / 2.;
      verts[7] = (rho + dRho / 2.) * sign;
      verts[9] = z - dZ / 2.;
      verts[10] = (rho + dRho / 2.) * sign;
      verts[2] = verts[5] = verts[8] = verts[11] = projectedEventDepth;
      qs_rhoz[icl]->AddQuad(verts);
      qs_rhoz[icl]->QuadValue((int)(1000 * energy));
      qs_rhoz[icl]->QuadId(new TNamed(Form("%s %d cell W%lu Rho%lu Z%d", clusterType.c_str(), icl,
                                           DetectorGeometry::ECalEndCapWheel(cellID),
                                           DetectorGeometry::ECalEndCapRhoBin(cellID), layer), "Dong!"));
      continue;
    }
    if (system == DetectorGeometry::kHCalEndcap)
    {
      // the cell is a rectangle in rho-z, from the segmentation
      double rhoMin, rhoMax, zLow, zHigh, phiLow, phiHigh;
      if (!geomReader->HCalEndCapCellBounds(cellID, rhoMin, rhoMax, zLow, zHigh, phiLow, phiHigh))
      {
        std::cout << "WARNING: HCAL endcap cell " << cellID << " outside of the segmentation, skipping cell" << std::endl;
        continue;
      }
      float sign = MathUtils::sign(y_center);
      verts[0] = zLow;
      verts[1] = rhoMin * sign;
      verts[3] = zHigh;
      verts[4] = rhoMin * sign;
      verts[6] = zHigh;
      verts[7] = rhoMax * sign;
      verts[9] = zLow;
      verts[10] = rhoMax * sign;
      verts[2] = verts[5] = verts[8] = verts[11] = projectedEventDepth;
      qs_rhoz[icl]->AddQuad(verts);
      qs_rhoz[icl]->QuadValue((int)(1000 * energy));
      qs_rhoz[icl]->QuadId(new TNamed(Form("%s %d cell L%d Th%d", clusterType.c_str(), icl, layer, thetaID), "Dong!"));
      continue;
    }
    float theta_center = atan2(r_center, z_center);
    if (system == DetectorGeometry::kECalBarrel)
    {
      thetaMin = theta_center - thetaGrid * geomReader->mergedCells_Theta[layer] / 2.;
      thetaMax = theta_center + thetaGrid * geomReader->mergedCells_Theta[layer] / 2.;
    }
    else if (system == DetectorGeometry::kHCalBarrel)
    {
      thetaMin = theta_center - thetaGridHCal / 2.;
      thetaMax = theta_center + thetaGridHCal / 2.;
    }
    else
      continue;

    verts[0] = r_in / tan(thetaMin);
    verts[1] = r_in * MathUtils::sign(y_center);
    verts[3] = r_in / tan(thetaMax);
    verts[4] = r_in * MathUtils::sign(y_center);
    verts[6] = r_out / tan(thetaMax);
    verts[7] = r_out * MathUtils::sign(y_center);
    verts[9] = r_out / tan(thetaMin);
    verts[10] = r_out * MathUtils::sign(y_center);
    verts[2] = verts[5] = verts[8] = verts[11] = projectedEventDepth;
    qs_rhoz[icl]->AddQuad(verts);
    qs_rhoz[icl]->QuadValue((int)(1000 * energy));
    qs_rhoz[icl]->QuadId(new TNamed(Form("%s %d cell L%d Th%d", clusterType.c_str(), icl, layer, thetaID), "Dong!"));
  }

  // cluster cells in rho-phi projection
  if (debug)
    std::cout << "  Draw rho-phi view of clustered cells" << std::endl;
  for (auto &[key, energy] : cellEnergies_rhophi)
  {
    int system = key.system;
    int icl = key.cluster;
    if (system == DetectorGeometry::kECalBarrel)
    {
      int layer = key.layer;
      int moduleID = key.bin;
      // std::cout << "cluster layer module = " << icl << " " << layer << " " << moduleID << std::endl;
      float r_in = geomReader->r[layer];
      float r_out = geomReader->r[layer + 1];
      double Lin = geomReader->getL(alpha, rMin, r_in);
      double Lout = geomReader->getL(alpha, rMin, r_out);
      double deltaL = rMin * sin(gridPhi / 2.0) * sin(alpha);
      for (int j = 0; j < geomReader->mergedModules[layer]; j++)
      {
        int iModule = moduleID + j;
        double phi0 = phiMin + iModule * gridPhi;
        verts[0] = rMin * cos(phi0) + (Lin + deltaL) * cos(phi0 + alpha);
        verts[1] = rMin * sin(phi0) + (Lin + deltaL) * sin(phi0 + alpha);
        verts[3] = rMin * cos(phi0) + (Lout + deltaL) * cos(phi0 + alpha);
        verts[4] = rMin * sin(phi0) + (Lout + deltaL) * sin(phi0 + alpha);
        verts[6] = rMin * cos(phi0 + gridPhi) + (Lout - deltaL) * cos(phi0 + gridPhi + alpha);
        verts[7] = rMin * sin(phi0 + gridPhi) + (Lout - deltaL) * sin(phi0 + gridPhi + alpha);
        verts[9] = rMin * cos(phi0 + gridPhi) + (Lin - deltaL) * cos(phi0 + gridPhi + alpha);
        verts[10] = rMin * sin(phi0 + gridPhi) + (Lin - deltaL) * sin(phi0 + gridPhi + alpha);
        verts[2] = verts[5] = verts[8] = verts[11] = projectedEventDepth;
        qs_rhophi[icl]->AddQuad(verts);
        qs_rhophi[icl]->QuadValue((int)(1000 * energy));
        // qs_rhophi[icl]->QuadId( new TNamed(Form("%s cell %lu", clusterType.c_str(), cellID), "Dong!"));
        qs_rhophi[icl]->QuadId(new TNamed(Form("%s %d cell L%d M%d", clusterType.c_str(), icl, layer, moduleID), "Dong!"));
      }
    }
    else if (system == DetectorGeometry::kHCalBarrel)
    {
      int layer = key.layer;
      int phiID = key.bin;
      // std::cout << "cluster layer module = " << icl << " " << layer << " " << moduleID << std::endl;
      float r_in = geomReader->rHCal[layer];
      float r_out = geomReader->rHCal[layer + 1];
      float phiLow = phiMinHCal + gridPhiHCal * phiID;
      float phiHigh = phiMinHCal + gridPhiHCal * (phiID + 1);
      verts[0] = r_in * cos(phiLow);
      verts[1] = r_in * sin(phiLow);
      verts[3] = r_in * cos(phiHigh);
      verts[4] = r_in * sin(phiHigh);
      verts[6] = r_out * cos(phiHigh);
      verts[7] = r_out * sin(phiHigh);
      verts[9] = r_out * cos(phiLow);
      verts[10] = r_out * sin(phiLow);
      verts[2] = verts[5] = verts[8] = verts[11] = projectedEventDepth;
      qs_rhophi[icl]->AddQuad(verts);
      qs_rhophi[icl]->QuadValue((int)(1000 * energy));
      qs_rhophi[icl]->QuadId(new TNamed(Form("%s %d cell L%d P%d", clusterType.c_str(), icl, layer, phiID), "Dong!"));
    }
  }

  // largest cell energy in each view, for the range of the palettes
  for (TEveBoxSet *b : bs)
    maxCellValue3D = std::max(maxCellValue3D, maxDigitValue(b));
  for (TEveQuadSet *q : qs_rhophi)
    maxCellValueRhoPhi = std::max(maxCellValueRhoPhi, maxDigitValue(q));
  for (TEveQuadSet *q : qs_rhoz)
    maxCellValueRhoZ = std::max(maxCellValueRhoZ, maxDigitValue(q));
}

void EventDisplay::drawEvent(int event)
{

  if (event != -1)
    eventId = event;

  printf("Loading event %d ...\n", eventId);

  // write event number in status bar
  browser->SetStatusText(Form("Loading event %d ...", eventId),0);

  // read event from file
  eventReader->loadEvent(eventId);

  // find out pdgID and momentum of primary particle
  // not OK for processes like Z->ee....
  // maybe can read process info from metadata instead?
  // (gen particle branches are not set up if drawGenParticles is false or
  // the collection is missing from the file)
  double pmax = 0.0;
  int ipmax = -1;
  unsigned int nGenParticles = eventReader->genParticles.generatorStatus ? eventReader->genParticles.generatorStatus->GetSize() : 0;
  for (unsigned int i = 0; i < nGenParticles; i++)
  {
    if ((*eventReader->genParticles.generatorStatus)[i] != 1)
      continue;
    float _px = (*eventReader->genParticles.momentum_x)[i];
    float _py = (*eventReader->genParticles.momentum_y)[i];
    float _pz = (*eventReader->genParticles.momentum_z)[i];
    float _p = sqrt(_px * _px + _py * _py + _pz * _pz);
    // std::cout << p << std::endl;
    if (_p == 0.)
      continue;
    if (_p > pmax)
    {
      pmax = _p;
      ipmax = i;
    }
  }
  int pdgID = (ipmax >= 0) ? (*eventReader->genParticles.PDG)[ipmax] : 0;

  // draw the various event elements
  drawParticles();
  drawHits();
  drawDigis();
  drawTracks();
  maxCellValue3D = maxCellValueRhoPhi = maxCellValueRhoZ = 0;
  for (ClusterSet &set : clusterSets)
  {
    fillClusters(set);
    drawClusters(set);
  }
  // adapt the range of the cell palettes to the event
  setCellPaletteRange(cellPalette3D, maxCellValue3D);
  setCellPaletteRange(cellPaletteRhoPhi, maxCellValueRhoPhi);
  setCellPaletteRange(cellPaletteRhoZ, maxCellValueRhoZ);

  // add process and event label
  std::wstring processStr = L"no primary particle";
  if (ipmax >= 0)
  {
    wchar_t energy[32];
    swprintf(energy, 32, L", %.1f GeV", pmax);
    processStr = ParticleUtils::name(pdgID) + energy;
  }
  if (processLabel == nullptr)
    processLabel = new TGLTextLabel(mainGLView, processStr, 0.05, 0.95, 0.05); // size: % of window height
  else
    processLabel->SetText(processStr);
  std::wstring eventStr = L"Event " + std::to_wstring(eventId);
  if (eventLabel == nullptr)
    eventLabel = new TGLTextLabel(mainGLView, eventStr, 0.05, 0.90, 0.05);
  else
    eventLabel->SetText(eventStr);
  // the labels are not elements of the scene: changing them does not trigger a redraw
  mainGLView->RequestDraw();

  TEveElement *top = (TEveElement *)gEve->GetCurrentEvent();
  // if nothing is drawn, top is null
  if (top)
  {
    // the event is projected in front of the geometry (see projectedEventDepth)
    for (TEveProjectionManager *manager : {rhoPhiProjManager, rhoZProjManager})
      manager->SetCurrentDepth(projectedEventDepth);

    rhoPhiEventScene->DestroyElements();
    rhoPhiProjManager->ImportElements(top, rhoPhiEventScene);

    rhoZEventScene->DestroyElements();
    rhoZProjManager->ImportElements(top, rhoZEventScene);

    for (TEveProjectionManager *manager : {rhoPhiProjManager, rhoZProjManager})
      manager->SetCurrentDepth(0.);

    attachClusterProjections();
  }

  std::cout << "Done" << std::endl
            << std::endl;

  browser->SetStatusText(Form("Event %d loaded", eventId),0);
  updateEventEntries();
}

// do not implement, use compiler-generated one
// EventDisplay::EventDisplay()
// {
// }

void EventDisplay::startDisplay(int initialEvent)
{
  // read the display setup from the configuration
  readConfiguration();

  // create magnetic field
  createMagneticField();

  std::cout << "******************************************************************************" << std::endl;
  std::cout << "Displaying the geometry" << std::endl;
  std::cout << "******************************************************************************" << std::endl
            << std::endl;

  // create the eve manager, the main window and the 3D view
  createMainWindow();

  // create the geometry and the readouts
  if (!buildGeometry())
    return;

  // create the projections; the geometry is projected when they are first
  // shown, but it is needed now to write the geometry cache (see below)
  createRhoPhiView();
  createRhoZView();

  gEve->DoRedraw3D();
  gEve->FullRedraw3D(true);

  // set the 3D view as the active tab, rotate the camera and set the clipping
  // (after the first drawing, which sets the camera from the bounding box of
  // the scene), so that the final 3D view is shown also while the
  // projections below are computed
  setInitialCamera();

  if (geometryExtract)
  {
    buildRhoPhiProjection();
    buildRhoZProjection();
  }

  // write the geometry cache, now that the volumes that can be drawn are known
  writeGeometryCache();

  // create the gui for event navigation
  makeGui();

  //
  // display the events
  //
  if (debug)
    std::cout << "evtFile: " << evtFile << std::endl;
  if (evtFile != "")
  {
    openEventFile();

    // load and display the requested event
    std::cout << std::endl;
    std::cout << "******************************************************************************" << std::endl;
    std::cout << "Reading the events" << std::endl;
    std::cout << "******************************************************************************" << std::endl
         << std::endl;
    drawEvent(initialEvent);
  }

  // legends of the cluster cell energies
  addPaletteLegends();
}

void EventDisplay::readConfiguration()
{
  // read some display setup from config file
  ClusterEnergyThreshold = displayConfig.getFloatConfig("energyThresholdClusters");
  ParticleEnergyThreshold = displayConfig.getFloatConfig("energyThresholdParticles");
  HitEnergyThreshold = displayConfig.getFloatConfig("energyThresholdHits");
  CellEnergyThreshold = displayConfig.getFloatConfig("energyThresholdCells");
  showReadout = displayConfig.getBoolConfig("showReadout");

  // calculate the geometry parameters
  // detectorVersion: 2 = ECAL barrel with 12 layers (ALLEGRO o1_v02)
  //                  3 = ECAL barrel with 11 layers (ALLEGRO o1_v03, o2_v01)
  if (geomFile.find("o1_v03") != std::string::npos || geomFile.find("o2_v01") != std::string::npos)
    detectorVersion = 3;
  else
    detectorVersion = 2;
  geomReader = std::make_unique<DetectorGeometry>(detectorVersion);
}

void EventDisplay::createMagneticField()
{
  // create magnetic field
  // bField: GlobalSolenoid for an ideal solenoid with the default field, or
  // GlobalSolenoid_<B>T for a field of <B> tesla inside the solenoid;
  // otherwise the path of the field map, or empty for the default field map
  std::string bField = displayConfig.getStringConfig("bField");
  const std::string solenoidPrefix = "GlobalSolenoid";
  const bool useSolenoid = (bField.compare(0, solenoidPrefix.size(), solenoidPrefix) == 0);
  double bSolenoid = geomReader->Bin;
  if (useSolenoid && bField != solenoidPrefix) {
    const std::string value = bField.substr(solenoidPrefix.size());
    const char *begin = value.c_str() + 1;
    char *end = nullptr;
    double b = 0.;
    if (value.size() > 2 && value.front() == '_' && value.back() == 'T')
      b = std::strtod(begin, &end);
    if (end != nullptr && end != begin && end == value.c_str() + value.size() - 1)
      bSolenoid = b;
    else
      std::cerr << "ERROR: invalid bField " << bField << " (expected GlobalSolenoid or GlobalSolenoid_<B>T, e.g. GlobalSolenoid_2.5T), "
                << "using a solenoid with B = " << bSolenoid << " T" << std::endl;
  }
  if (!useSolenoid) {
    if (bField == "")
      bField = "data/ALLEGRO_fieldmap_2T_20260424.root";
    std::string fieldMapFile = PathUtils::findFile(bField);
    std::cout << "Using magnetic field map " << fieldMapFile << std::endl;
    try {
      CylField* fieldMap = new CylField(fieldMapFile, "fieldmap");
      if (displayConfig.getBoolConfig("printFieldMap")) {
        fieldMap->PrintFieldMap();
        fieldMap->DrawFieldMap("field.pdf", 20, 20, 20.0);
      }
      // Note: field is reversed in EVE (see https://root.cern/doc/v636/classTEveTrackPropagator.html)
      fieldMap->setReverseState(true);
      magField.reset(fieldMap);
    }
    catch (const std::exception& e) {
      std::cerr << "ERROR: cannot read the magnetic field map " << fieldMapFile << " (" << e.what()
                << "), using the ideal solenoid magnetic field instead" << std::endl;
    }
  }
  if (!magField) {
    // the return field outside the solenoid scales with the field inside
    // (+ 0. avoids printing -0 when there is no field)
    const double bOut = bSolenoid * geomReader->Bout / geomReader->Bin + 0.;
    std::cout << "Using ideal solenoid magnetic field, B = " << bSolenoid << " T inside, "
              << bOut << " T outside" << std::endl;
    MagField* solenoid = new MagField(bSolenoid, bOut, geomReader->zMax, geomReader->rMax);
    // Note: field is reversed in EVE (see https://root.cern/doc/v636/classTEveTrackPropagator.html)
    solenoid->setReverseState(true);
    magField.reset(solenoid);
  }
  std::cout << "\n" << std::endl;
}

void EventDisplay::createMainWindow()
{
  // create the eve manageer
  TEveManager::Create();
  browser = gEve->GetBrowser();

  // Set title of main window
  browser->SetWindowName("FCC-ee ALLEGRO detector event display");

  // see palettes here: https://root.cern.ch/doc/master/classTColor.html
  // gStyle->SetPalette(kAvocado);
  gStyle->SetPalette(kSienna);

  // palettes for the cluster cells in each view (the energies summed in the
  // projected cells are larger than those of the cells in 3D)
  cellPalette3D = makeCellPalette();
  cellPaletteRhoPhi = makeCellPalette();
  cellPaletteRhoZ = makeCellPalette();

  // first tab
  mainGLView = gEve->GetDefaultGLViewer();
  int bkgColor = displayConfig.getIntConfig("bkgColor3D");
  mainGLView->ColorSet().Background().SetColor(bkgColor>=0 ? bkgColor : kBlack);
  mainGLView->SetGuideState(TGLUtil::kAxesOrigin, false, false, 0);
  mainGLView->DrawGuides();
  gEve->GetDefaultViewer()->SetElementName("3D view");
  mainGLView->CurrentCamera().RotateRad(-.8, 0.5);
  browser->GetTabRight()->Connect("Selected(Int_t)", "EventDisplay", this, "onTabSelected(Int_t)");
  // (emitted when the window is closed, before TEve is terminated)
  browser->Connect("CloseWindow()", "EventDisplay", this, "onCloseWindow()");
}

bool EventDisplay::buildGeometry()
{
  // Create the geometry and the readouts
  geom = new TEveElementList("Geometry");
  readout = new TEveElementList("Readout");
  readout->SetRnrSelfChildren(showReadout, showReadout);
  if (useG4geom)
  {
    if (!buildG4Geometry())
      return false;
  }
  else
    buildSimplifiedGeometry();
  gEve->AddToListTree(readout, true);
  gEve->FullRedraw3D(true);
  return true;
}

bool EventDisplay::buildG4Geometry()
{
  gEve->AddGlobalElement(geom);
  gEve->AddToListTree(geom, true);

  // the geometry without the volumes that are never drawn, from the cache
  // if available, otherwise the full geometry from the geometry file
  TEveGeoShapeExtract *gse = nullptr;
  if (useGeometryCache)
    gse = GeometryCache::read(geomFile, geometryCacheTag(), geometryCodeVersion);
  const bool fromCache = (gse != nullptr);
  if (!fromCache)
  {
    std::cout << "Reading Geant4 geometry from file " << geomFile << std::endl;
    std::unique_ptr<TFile> fGeom(TFile::Open(geomFile.c_str(), "READ"));
    if (!fGeom)
      return false;
    gse = dynamic_cast<TEveGeoShapeExtract *>(fGeom->Get("world"));
    if (!gse)
      return false;
    fGeom->Close();
  }

  TEveGeoShape *world = TEveGeoShape::ImportShapeExtract(gse, 0);
  world->SetMainTransparency(100);
  if (useGeometryCache && !fromCache)
  {
    // keep the full geometry until the display is set up, to write the cache
    geometryExtract = gse;
    GeometryCache::mapExtracts(gse, world, geometryExtracts);
  }
  else
    GeometryCache::deleteExtract(gse);

  // check that the number of ECAL barrel layers in the geometry matches detectorVersion
  int nLayersInGeom = countECalBarrelLayers(world);
  if (nLayersInGeom < 0)
    std::cout << "WARNING: could not find the ECAL barrel layers in the geometry, cannot check detectorVersion" << std::endl;
  else if (nLayersInGeom != geomReader->nLayers)
    std::cout << "WARNING: geometry has " << nLayersInGeom << " ECAL barrel layers but detectorVersion "
              << detectorVersion << " (from the geometry file name) assumes " << geomReader->nLayers << std::endl;

  if (not showFullDetector) {
    rMax = doHCal ? geomReader->rMaxHCal : geomReader->rMax;
  }
  else {
    rMax = 5000.*mm;
  }

  if (showFullDetector)
    buildFullDetector(world);
  else if (!buildCalorimeters(world))
    return false;
  return true;
}

void EventDisplay::buildFullDetector(TEveElement *world)
{
  m_subDetectors.clear();
  m_subDetectors.reserve(20);
  auto simple = [this](TEveElement* e, SubDetector& d)
  {
    addSimpleElement(e, d);
  };
  m_subDetectors.push_back({"MDI", "", {"ScreenSol","CompSol"}, nullptr, false,
      simple});
  m_subDetectors.push_back({"Beampipe", "", {"BeamPipe","Beampipe"}, nullptr, false,
      simple});
  m_subDetectors.push_back({"LumiCal", "LumiCal_*", {}, nullptr, false,
        simple});
  m_subDetectors.push_back({"Vertex", "", {"Vertex"}, nullptr, true,
    [this](TEveElement* e, SubDetector& d) { addVertex(e, d); }});
  m_subDetectors.push_back({"Drift chamber", "DCH_*", {}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addDCH(e, d); }});
  m_subDetectors.push_back({"Straw tube tracker", "STT_*", {}, nullptr, true,
                [this](TEveElement* e, SubDetector& d) { addDCH(e, d); }});
  m_subDetectors.push_back({"Si wrapper barrel", "SiWrB_*", {}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addWrapper(e, d); }});
  m_subDetectors.push_back({"Si wrapper endcap", "SiWrD_*", {}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addWrapper(e, d); }});
  m_subDetectors.push_back({"ECal barrel", "ECalBarrel*", {}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addECalB(e, d); }});
  m_subDetectors.push_back({"ECal endcap", "", {"ECalEndcaps"}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addECalEC(e, d); }});
  m_subDetectors.push_back({"HCal barrel", "HCalEnvelopeVolume*", {}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addHCalB(e, d); }});
  m_subDetectors.push_back({"HCal endcap", "HCal(\\w+)Endcap(\\w+)", {}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addHCalEC(e, d); }});
  m_subDetectors.push_back({"Muon tagger barrel", "MuonTaggerBarrel*", {}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addMuonB(e, d); }});
  m_subDetectors.push_back({"Muon tagger endcap", "MuonTaggerEndcap*", {}, nullptr, true,
        [this](TEveElement* e, SubDetector& d) { addMuonEC(e, d); }});

  for (TEveElement::List_i itr = world->BeginChildren(); itr != world->EndChildren(); itr++)
  {
    TEveElement *el = *itr;
    TString elName(el->GetElementName());
    SubDetector* det = findSubDetector(elName);
    if (det) {
      if (det->handler)
        det->handler(el, *det);
    }
  }
}

bool EventDisplay::buildCalorimeters(TEveElement *world)
{
  // lists of the elements in the LAr bath of the ECAL barrel, to turn them on/off together
  TEveElementList *PCBs = new TEveElementList("PCBs");
  TEveElementList *actives = new TEveElementList("Active elements");
  TEveElementList *passives = new TEveElementList("Passive elements");
  //
  // ECAL barrel
  //
  TPRegexp re("ECalBarrel*");
  TEveElement *ecalBarrelEl = world->FindChild(re);
  if (!ecalBarrelEl)
  {
    std::cerr << "ERROR: ECalBarrel volume not found in geometry file " << geomFile << std::endl;
    return false;
  }
  ecalBarrelEl->SetPickableRecursively(kTRUE);
  ((TEveGeoShape *)ecalBarrelEl)->SetNSegments(256);
  geom->AddElement(ecalBarrelEl);
  // hide the envelope, only show physical volumes in it
  ecalBarrelEl->SetRnrSelfChildren(false, true);

  // set transparency of the subvolumes of the bath
  re = TPRegexp("LAr_bath*");
  TEveElement *bath = ecalBarrelEl->FindChild(re);
  TEveElement::List_t cryoMatches;
  re = TPRegexp("ECAL_Cryo*");
  ecalBarrelEl->FindChildren(cryoMatches, re);
  for (TEveElement *a : cryoMatches)
  {
    a->SetMainTransparency(70);
    ((TEveGeoShape *)a)->SetNSegments(256);
  }
  TEveElement::List_t servicesMatches;
  re = TPRegexp("services*");
  ecalBarrelEl->FindChildren(servicesMatches, re);
  for (auto a : servicesMatches)
  {
    a->SetMainTransparency(70);
    ((TEveGeoShape *)a)->SetNSegments(256);
  }

  // turn off Cryo side in 3d, weird rendering
  re = TPRegexp("ECAL_Cryo_side_(\\w+)");
  TEveElement::List_t matches2;
  ecalBarrelEl->FindChildren(matches2, re);
  for (auto a : matches2)
  {
    a->SetRnrSelf(false);
  }
  // make lists of elements inside bath to turn on/off simultaneously
  if (bath)
  {
    TEveElementList *newbath = new TEveElementList("LAr_bath");
    ecalBarrelEl->AddElement(newbath);
    TEveElement::List_t bathMatches;
    re = TPRegexp("PCB(\\w+)");
    bath->FindChildren(bathMatches, re);
    for (auto a : bathMatches)
      PCBs->AddElement(a);
    newbath->AddElement(PCBs);
    TEveElement::List_t bathMatches2;
    re = TPRegexp("active*");
    bath->FindChildren(bathMatches2, re);
    for (TEveElement *a : bathMatches2)
    {
      actives->AddElement(a);
      ((TEveGeoShape *)a)->SetMainTransparency(100);
      // increase line width of active layers from 1.0 to 5.0
      // make them transparent so that passive elements are not covered
      TEveElement::List_t _matches;
      re = TPRegexp("layer*");
      a->FindChildren(_matches, re);
      for (TEveElement *l : _matches)
      {
        ((TEveGeoShape *)l)->SetLineWidth(5.0);
        ((TEveGeoShape *)l)->SetMainTransparency(100);
      }
    }
    newbath->AddElement(actives);
    TEveElement::List_t matches3;
    re = TPRegexp("passive*");
    bath->FindChildren(matches3, re);
    for (TEveElement *a : matches3)
    {
      passives->AddElement(a);
      // decrease transparency from 60 to 20
      //  ((TEveGeoShape*) a)->SetMainTransparency(20);
      //  do not draw individual constituents of passive elements
      //  a->SetRnrSelfChildren(true, false);
    }
    newbath->AddElement(passives);
    ecalBarrelEl->RemoveElement(bath);
    // hide elements inside bath by default because they are slow in 3D
    newbath->SetRnrSelfChildren(true, false);
  }

  //
  // ECAL endcaps
  //
  if (doEndcaps)
  {
    ecalendcap = new TEveElementList("ECal endcap");
    geom->AddElement(ecalendcap);

    TPRegexp re_ecalec("(\\w+)ECalEndcap(\\w+)");
    TEveElement::List_t matches_ec;
    world->FindChildren(matches_ec, re_ecalec);
    for (TEveElement *a : matches_ec)
    {
      ecalendcap->AddElement(a);
    }

    TPRegexp re_ecalec2("EMEC(\\w+)");
    TEveElement::List_t matches_ec2;
    world->FindChildren(matches_ec2, re_ecalec2);
    for (TEveElement *a : matches_ec2)
    {
      ecalendcap->AddElement(a);
      // do not draw envelope, only physical elements inside
      a->SetRnrSelfChildren(false, true);
      // hide elements inside layerEnvelope by default because they are slow in 3D
      // tree is EMEC_vol->EMEC_positive/negative_vol->layerEnvelope
      for (TEveElement::List_i itr = a->BeginChildren(); itr != a->EndChildren(); itr++)
      {
        TEveElement *emec_posneg_vol = (*itr);
        emec_posneg_vol->SetRnrSelfChildren(false, true);
        for (TEveElement::List_i itr2 = emec_posneg_vol->BeginChildren(); itr2 != emec_posneg_vol->EndChildren(); itr2++)
        {
          (*itr2)->SetRnrSelfChildren(true, false);
        }
      }
    }
  }

  // HCAL
  TEveElement *hcalBarrelEl = nullptr;
  if (doHCal)
  {
    re = TPRegexp("HCalEnvelopeVolume*");
    hcalBarrelEl = world->FindChild(re);
    if (!hcalBarrelEl)
      std::cout << "WARNING: HCalEnvelopeVolume not found in geometry file, not drawing the HCAL barrel" << std::endl;
  }
  if (hcalBarrelEl)
  {
    // add HCal envelope and subvolumes
    hcalBarrelEl->SetPickableRecursively(kTRUE);
    ((TEveGeoShape *)hcalBarrelEl)->SetNSegments(256);
    geom->AddElement(hcalBarrelEl);
    // do not draw envelope, only physical elements inside
    hcalBarrelEl->SetRnrSelfChildren(false, true);

    // re = TPRegexp("HCalLayerVol*");
    // hcalBarrelEl->FindChildren(matches, re);
    // for (auto a : matches) a->SetRnrSelfChildren(true, false);

    // set transparency and number of segments in 2D projection of plate faces and steel
    {
      re = TPRegexp("HCal(\\w+)PlateVol(\\w+)");
      TEveElement::List_t plateMatches;
      hcalBarrelEl->FindChildren(plateMatches, re);
      for (TEveElement *a : plateMatches)
      {
        a->SetMainTransparency(70);
        ((TEveGeoShape *)a)->SetNSegments(256);
      }
    }

    re = TPRegexp("HCal*Steel*");
    TEveElement::List_t steelMatches;
    hcalBarrelEl->FindChildren(steelMatches, re);
    for (auto a : steelMatches)
    {
      a->SetMainTransparency(70);
      ((TEveGeoShape *)a)->SetNSegments(256);
    }

    // group together the layers so that they can be turned on/off together
    // set number of segments of 2D projection of layer envelopes
    TEveElementList *hcalLayers = new TEveElementList("HCalLayers");
    hcalBarrelEl->AddElement(hcalLayers);
    re = TPRegexp("HCalLayerVol*");
    TEveElement::List_t matches4;
    hcalBarrelEl->FindChildren(matches4, re);
    for (TEveElement *a : matches4)
    {
      hcalLayers->AddElement(a);
      hcalBarrelEl->RemoveElement(a);
      ((TEveGeoShape *)a)->SetNSegments(256);
      // set number of segments of 2D projection of TileSequenceVol and TileVol
      // slow, don't do it!
      // instead, hide by default the single tiles because they are slow in 3D
      a->SetRnrSelfChildren(true, false);
      /*
      TPRegexp re2("HCalTileSequenceVol*");
      TEveElement::List_t matches_tsv;
      a->FindChildren(matches_tsv, re2);
      for (TEveElement* tsv : matches_tsv)
      {
        ((TEveGeoShape*) tsv)->SetNSegments(256);
        TPRegexp re3("HCalTileVol*");
        TEveElement::List_t matches_tv;
        tsv->FindChildren(matches_tv, re3);
        for (TEveElement* tv : matches_tv)
        {
          ((TEveGeoShape*) tv)->SetNSegments(256);
        }
      }
      */
    }
  }

  //
  // HCAL endcaps (also if the HCAL barrel is not found)
  //
  if (doHCal && doEndcaps)
  {
    hcalendcap = new TEveElementList("HCal endcap");
    geom->AddElement(hcalendcap);

    TPRegexp re_hcalec("HCal(\\w+)Endcap(\\w+)");
    TEveElement::List_t matches_ec;
    world->FindChildren(matches_ec, re_hcalec);
    for (TEveElement *a : matches_ec)
    {
      hcalendcap->AddElement(a);
      // do not draw envelope, only physical elements inside
      a->SetRnrSelfChildren(false, true);
      TEveElementList *hcalECLayers = new TEveElementList("HCalECLayers");
      hcalendcap->AddElement(hcalECLayers);
      re = TPRegexp("HCalECLayerVol*");
      TEveElement::List_t matches_layer;
      a->FindChildren(matches_layer, re);
      for (TEveElement *l : matches_layer)
      {
        hcalECLayers->AddElement(l);
        a->RemoveElement(l);
        l->SetRnrSelfChildren(true, false);
      }
    }
  }
  return true;
}

void EventDisplay::buildSimplifiedGeometry()
{
  std::cout << "Creating simplified geometry based on calculated dimensions " << std::endl;

  // the ECAL barrel envelope
  ecalbarrel = new TEveGeoShape("ECAL barrel");
  ecalbarrel->SetShape(new TGeoTube(geomReader->rMin, geomReader->rMax, geomReader->zMax));
  ecalbarrel->SetMainColor(kCyan);
  ecalbarrel->SetMainTransparency(90);
  ecalbarrel->SetNSegments(128);
  // geom->AddElement(ecalbarrel);
  gEve->AddGlobalElement(ecalbarrel);

  // the ECAL barrel layers
  TEveElementList *layers = new TEveElementList("layers");
  ecalbarrel->AddElement(layers);
  TEveGeoShape *b;
  for (int iLayer = 0; iLayer < geomReader->nLayers; iLayer++)
  {
    b = new TEveGeoShape(Form("ECAL barrel layer %d", iLayer));
    b->SetShape(new TGeoTube(geomReader->r[iLayer], geomReader->r[iLayer + 1], geomReader->zMax));
    b->SetMainColor(kCyan);
    b->SetMainTransparency(90);
    b->SetNSegments(128);
    layers->AddElement(b);
  }

  // the ECAL electrodes
  TEveElementList *modules = new TEveElementList("modules");
  ecalbarrel->AddElement(modules);
  for (int iModule = 0; iModule < geomReader->nModules; iModule++)
  {
    // double phi0 = iModule * geomReader->gridPhi - geomReader->gridPhi / 12.; // small extra shift is due to finite width of element (?)
    double phi0 = iModule * geomReader->gridPhi;
    double phi = phi0 + geomReader->dPhiAvg;
    b = new TEveGeoShape(Form("Module %d", iModule));
    b->SetShape(new TGeoBBox(geomReader->Ltot / 2, 0.01, (geomReader->zMax - geomReader->zMin) / 2.0));
    b->SetMainColor(kGray);
    b->SetMainTransparency(95);
    TGeoRotation *rot = new TGeoRotation();
    rot->SetAngles((geomReader->alpha + iModule * geomReader->gridPhi) * 180. / TMath::Pi(), 0., 0.);
    TGeoCombiTrans *c1 = new TGeoCombiTrans(geomReader->rAvg * cos(phi), geomReader->rAvg * sin(phi), 0.0, rot);
    b->SetTransMatrix(*c1);
    modules->AddElement(b);
  }
  gEve->AddToListTree(ecalbarrel, true);

  // the HCAL barrel envelope
  if (doHCal)
  {
    hcalbarrel = new TEveGeoShape("HCAL barrel");
    hcalbarrel->SetShape(new TGeoTube(geomReader->rMinHCal, geomReader->rMaxHCal, geomReader->zMaxHCal));
    hcalbarrel->SetMainColor(kRed);
    hcalbarrel->SetMainTransparency(80);
    hcalbarrel->SetNSegments(128);
    // geom->AddElement(hcalbarrel);
    gEve->AddGlobalElement(hcalbarrel);

    // the HCAL barrel layers
    TEveElementList *hcallayers = new TEveElementList("layers");
    hcalbarrel->AddElement(hcallayers);
    for (int iLayer = 0; iLayer < geomReader->nLayersHCal; iLayer++)
    {
      b = new TEveGeoShape(Form("HCAL barrel layer %d", iLayer));
      b->SetShape(new TGeoTube(geomReader->rHCal[iLayer], geomReader->rHCal[iLayer + 1], geomReader->zMaxHCal));
      b->SetMainColor(kRed);
      b->SetMainTransparency(80);
      b->SetNSegments(128);
      hcallayers->AddElement(b);
    }
    gEve->AddToListTree(hcalbarrel, true);
  }
}

void EventDisplay::createRhoPhiView()
{
  // create second tab (R-phi view)
  rhoPhiView = gEve->SpawnNewViewer("Projection Rho-Phi");

  // two scenes, for geometry and event
  rhoPhiScene = gEve->SpawnNewScene("Rho-Phi geometry",
                                    "Scene holding projected geometry data for the Rho-Phi view.");
  rhoPhiView->AddScene(rhoPhiScene);
  if (evtFile != "")
  {
    // (added to the view when the geometry is projected, see buildProjection)
    rhoPhiEventScene = gEve->SpawnNewScene("Rho-Phi Event Data",
                                           "Scene holding projected event-data for the Rho-Phi view.");
  }

  rhoPhiGLView = rhoPhiView->GetGLViewer();
  // set camera orientation
  rhoPhiGLView->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  // create 3D->2D projection manager for rho-phi
  rhoPhiProjManager = new TEveProjectionManager();
  rhoPhiProjManager->SetProjection(TEveProjection::kPT_RPhi);
  auto axes = new TEveProjectionAxes(rhoPhiProjManager);
  axes->SetElementName("Rho-Phi projection axes");
  axes->SetDrawOrigin(true);
  rhoPhiScene->AddElement(axes);
  // the geometry and the readout are projected when the view is first shown
  // (see buildRhoPhiProjection)
}

void EventDisplay::buildRhoPhiProjection()
{
  if (rhoPhiProjectionBuilt)
    return;
  rhoPhiProjectionBuilt = true;
  buildProjection("Rho-Phi", rhoPhiView, rhoPhiScene, rhoPhiEventScene, [this]() {
    projectGeometryRhoPhi();
    drawRhoPhiReadout();
  });
}

void EventDisplay::buildRhoZProjection()
{
  if (rhoZProjectionBuilt)
    return;
  rhoZProjectionBuilt = true;
  buildProjection("Rho-Z", rhoZView, rhoZScene, rhoZEventScene, [this]() {
    projectGeometryRhoZ();
    drawRhoZReadout();
  });
}

void EventDisplay::buildProjection(const char *view, TEveViewer *viewer, TEveScene *scene, TEveScene *eventScene,
                                   const std::function<void()> &build)
{
  std::cout << "Projecting the geometry in the " << view << " view (slow, only first event)" << std::endl;
  if (browser)
  {
    browser->SetStatusText(Form("Projecting the geometry in the %s view ...", view), 0);
    // redraw only the status bar now: processing the events here (e.g. with
    // gSystem->ProcessEvents) could run other commands in the middle
    gClient->NeedRedraw(browser->GetStatusBar(), kTRUE);
    gVirtualX->Update(1);
  }
  // the readouts added now follow the state of the readout button
  const int nReadouts = readout->NumChildren();
  build();
  int i = 0;
  for (TEveElement::List_i itr = readout->BeginChildren(); itr != readout->EndChildren(); ++itr, ++i)
    if (i >= nReadouts && !showReadout)
      (*itr)->SetRnrSelf(false);
  // rebuild the GL representation of the scene (redrawing the viewer is not
  // enough), as done for all the scenes by FullRedraw3D at startup
  scene->Repaint(kTRUE);
  // show the projected event only now, not without the geometry
  if (eventScene)
    viewer->AddScene(eventScene);
  if (browser)
    browser->SetStatusText(evtFile != "" ? Form("Event %d loaded", eventId) : "", 0);
}

void EventDisplay::projectGeometryRhoPhi()
{
  if (useG4geom)
    rhoPhiProjManager->ImportElements(geom, rhoPhiScene);
  else
  {
    rhoPhiProjManager->ImportElements(ecalbarrel, rhoPhiScene);
    if (doHCal)
      rhoPhiProjManager->ImportElements(hcalbarrel, rhoPhiScene);
  }

  TEveElement *element = rhoPhiScene->FindChild("Geometry [P]");
  if (!element)
    return;

  // do not draw the cryo side in the rho-phi view as it renders poorly
  // and turn on the active elements so we can see the cell edges
  showECalBathCells(element, true);

  // when drawing the full detector - do not draw endcap volumes in rho-phi view
  // and show active volumes in ecal barrel
  if (!showFullDetector)
    return;
  for (TEveElement::List_i itr = element->BeginChildren(); itr != element->EndChildren(); itr++)
  {
    TEveElement *a = *itr;
    TString s(a->GetElementName());
    // do not draw endcap volumes
    if (s.Contains("endcap"))
      a->SetRnrSelfChildren(false, false);

    // turn off the envelopes of the vertex barrel
    else if (s.BeginsWith("Vertex barrel"))
      showVertexBarrelLayers(a);

    // show active volumes in ecal barrel
    else if (s.BeginsWith("ECal barrel"))
      showECalBarrelContent(a, true);

    // show layers in HCAL barrel
    else if (s.BeginsWith("HCal barrel"))
    {
      if (a->NumChildren() == 0)
        continue;
      showEnvelopeContent(*(a->BeginChildren()), [](TEveElement *el, const TString &elname) {
        if (elname.Contains("EndPlate"))
          el->SetRnrSelfChildren(false, false);
        else {
          el->SetRnrSelfChildren(true, false);
          if (elname.Contains("LayerVol"))
            el->SetMainColor(kAzure-7);
          el->SetMainTransparency(0); // check if this does not hide the hits
        }
      });
    }

    // show layers in MUON barrel
    else if (s.BeginsWith("Muon tagger barrel"))
    {
      if (a->NumChildren() == 0)
        continue;
      showMuonLayers(*(a->BeginChildren()));
    }
  }
}

void EventDisplay::drawRhoPhiReadout()
{
  // draw the merged ECAL readout segmentation in rho-phi
  TEveStraightLineSet *gridmod = new TEveStraightLineSet("ECAL phi readout merged");
  gridmod->SetLineColor(kViolet);
  gridmod->SetLineWidth(3);
  for (int iLayer = 0; iLayer < geomReader->nLayers; iLayer++)
  {
    double Lin = geomReader->getL(geomReader->alpha, geomReader->rMin, geomReader->r[iLayer]);
    double Lout = geomReader->getL(geomReader->alpha, geomReader->rMin, geomReader->r[iLayer + 1]);
    for (int iModule = 0; iModule < geomReader->nModules; iModule++)
    {
      if (iModule % geomReader->mergedModules[iLayer] != 0)
        continue;
      double phi0 = geomReader->phiMin + iModule * geomReader->gridPhi;
      double x1 = geomReader->rMin * cos(phi0) + Lin * cos(phi0 + geomReader->alpha);
      double y1 = geomReader->rMin * sin(phi0) + Lin * sin(phi0 + geomReader->alpha);
      double x2 = geomReader->rMin * cos(phi0) + Lout * cos(phi0 + geomReader->alpha);
      double y2 = geomReader->rMin * sin(phi0) + Lout * sin(phi0 + geomReader->alpha);
      gridmod->AddLine(x1, y1, 0.0, x2, y2, 0.0);
    }
  }
  rhoPhiScene->AddElement(gridmod);
  readout->AddElement(gridmod);

  // draw the HCAL readout segmentation in rho-phi
  if (doHCal)
  {
    TEveStraightLineSet *hcalPhiReadout = new TEveStraightLineSet("HCAL phi readout");
    hcalPhiReadout->SetLineColor(kViolet);
    hcalPhiReadout->SetLineWidth(5);
    double r1 = geomReader->rMinHCal;
    double r2 = geomReader->rMaxHCal;
    for (int iLayer = 0; iLayer < geomReader->nLayersHCal; iLayer++)
    {
      for (int i = 0; i < geomReader->nPhiBinsHCal; i++)
      {
        double phi = geomReader->phiMinHCal + i * geomReader->gridPhiHCal;
        double x1 = r1 * cos(phi);
        double y1 = r1 * sin(phi);
        double x2 = r2 * cos(phi);
        double y2 = r2 * sin(phi);
        hcalPhiReadout->AddLine(x1, y1, 0., x2, y2, 0.);
      }
    }
    rhoPhiScene->AddElement(hcalPhiReadout);
    readout->AddElement(hcalPhiReadout);
  }

  // the muon tagger readout in the barrel (with the full detector, where it is
  // drawn): the edges of the phi bins
  if (showFullDetector)
  {
    TEveStraightLineSet *muonPhiReadout = new TEveStraightLineSet("Muon tagger phi readout");
    muonPhiReadout->SetLineColor(kViolet);
    const double r1 = geomReader->rMinMuonBarrel;
    const double r2 = geomReader->rMaxMuonBarrel;
    const double gridPhi = TMath::TwoPi() / geomReader->nPhiBinsMuon;
    for (int i = 0; i < geomReader->nPhiBinsMuon; i++)
    {
      const double phi = geomReader->phiOffsetMuon + (i - 0.5) * gridPhi;
      muonPhiReadout->AddLine(r1 * cos(phi), r1 * sin(phi), 0., r2 * cos(phi), r2 * sin(phi), 0.);
    }
    rhoPhiScene->AddElement(muonPhiReadout);
    readout->AddElement(muonPhiReadout);
  }
}

void EventDisplay::createRhoZView()
{
  // third tab (R-z view)
  rhoZView = gEve->SpawnNewViewer("Projection Rho-Z");
  rhoZScene = gEve->SpawnNewScene("Rho-Z geometry",
                                  "Scene holding projected geometry data for the Rho-Z view.");
  rhoZView->AddScene(rhoZScene);
  if (evtFile != "")
  {
    // (added to the view when the geometry is projected, see buildProjection)
    rhoZEventScene = gEve->SpawnNewScene("Rho-Z Event Data",
                                         "Scene holding projected event-data for the Rho-Z view.");
  }
  rhoZGLView = rhoZView->GetGLViewer();
  // set camera orientation
  rhoZGLView->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  // create 3D->2D projection manager for rho-z
  rhoZProjManager = new TEveProjectionManager();
  rhoZProjManager->SetProjection(TEveProjection::kPT_RhoZ);
  auto axes2 = new TEveProjectionAxes(rhoZProjManager);
  axes2->SetElementName("Rho-Z projection axes");
  axes2->SetDrawOrigin(true);
  rhoZScene->AddElement(axes2);
  // the geometry and the readout are projected when the view is first shown
  // (see buildRhoZProjection)
}

void EventDisplay::projectGeometryRhoZ()
{
  if (useG4geom)
    rhoZProjManager->ImportElements(geom, rhoZScene);
  else
  {
    rhoZProjManager->ImportElements(ecalbarrel, rhoZScene);
    if (doHCal)
      rhoZProjManager->ImportElements(hcalbarrel, rhoZScene);
  }

  TEveElement *element = rhoZScene->FindChild("Geometry [P]");
  if (!element)
    return;

  // turn on the active elements in rho-Z so we can see the cell edges
  showECalBathCells(element, false);
  // show the content of the LAr bath of the ECAL endcap
  for (TEveElement::List_i itr = element->BeginChildren(); itr != element->EndChildren(); itr++)
    if (TString((*itr)->GetElementName()).BeginsWith("ECal endcap"))
      showECalEndcapContent(*itr);

  if (!showFullDetector)
    return;
  for (TEveElement::List_i itr = element->BeginChildren(); itr != element->EndChildren(); itr++) {
    TEveElement *a = *itr;
    TString s(a->GetElementName());

    // MDI, beampipe, lumical can be enabled
    if (s.BeginsWith("MDI") || s.BeginsWith("Beampipe") || s.BeginsWith("LumiCal")) {
      a->SetRnrSelfChildren(true, true);
    }
    // turn off the envelopes of the vertex barrel
    else if (s.BeginsWith("Vertex barrel")) {
      showVertexBarrelLayers(a);
    }
    // show vessel and gas for DCH
    else if (s.BeginsWith("Drift chamber")) {
      a->SetRnrSelfChildren(false, true);
      for (TEveElement::List_i itrChild = a->BeginChildren(); itrChild != a->EndChildren(); itrChild++) {
        TEveElement *b = *itrChild;
        b->SetRnrSelfChildren(true, false);
        if (overrideDetectorColors) b->SetMainColor(kGreen-5);
        if (TEveShape *shape = dynamic_cast<TEveShape *>(b))
          shape->SetDrawFrame(false);
      }
    }
    // turn off the Si wrapper endcap envelope since it does not render well
    // instead show its daughter volume instead
    else if (s.BeginsWith("Si wrapper endcap")) {
      a->SetRnrSelfChildren(false, true);
      for (TEveElement::List_i itrChild = a->BeginChildren(); itrChild != a->EndChildren(); itrChild++) {
        TEveElement *b = *itrChild;
        TString childName(b->GetElementName());
        if (childName.BeginsWith("SiWrD_envelope")) {
          b->SetRnrSelfChildren(false, true);
          for (TEveElement::List_i itrGrandChild = b->BeginChildren(); itrGrandChild != b->EndChildren(); itrGrandChild++) {
            TEveElement *c = *itrGrandChild;
            c->SetRnrSelfChildren(true, false);
          }
        }
      }
    }
    // ECAL barrel
    else if (s.BeginsWith("ECal barrel")) {
      showECalBarrelContent(a, false);
    }
    // HCAL barrel - do not draw the envelope, rather show the layers
    else if (s.BeginsWith("HCal barrel"))
    {
      if (a->NumChildren() == 0)
        continue;
      showEnvelopeContent(*(a->BeginChildren()), [this](TEveElement *el, const TString &elname) {
        el->SetRnrSelfChildren(true, false);
        if (elname.Contains("HCalLayerVol")) {
          if (overrideDetectorColors) el->SetMainColor(kAzure-7);
          el->SetMainTransparency(0);
          // note: can do el->SetRnrSelfChildren(false, true); to see detailed row structure (tile/abs/..)
        }
      });
    }
    // HCAL endcap - do not draw the envelope
    else if (s.BeginsWith("HCal endcap"))
    {
      TPRegexp re("HCalThreePartsEndcap*");
      TEveElement *envelope = a->FindChild(re);
      if (!envelope)
        continue;
      showEnvelopeContent(envelope, [this](TEveElement *el, const TString &elname) {
        el->SetRnrSelfChildren(true, false);
        if (elname.Contains("HCalECLayerVol")) {
          if (overrideDetectorColors) el->SetMainColor(kAzure-7);
          el->SetMainTransparency(0);
          // note: can do el->SetRnrSelfChildren(false, true); to see detailed row structure (tile/abs/..)
        }
      });
    }
    // MUON barrel - do not draw the envelopes, rather show the layers
    else if (s.BeginsWith("Muon tagger barrel"))
    {
      if (a->NumChildren() == 0)
        continue;
      showMuonLayers(*(a->BeginChildren()));
    }
    // MUON endcap - do not draw the envelope, rather show the layers
    else if (s.BeginsWith("Muon tagger endcap"))
    {
      for (TEveElement::List_i itrenv = a->BeginChildren(); itrenv != a->EndChildren(); itrenv++)
        showMuonLayers(*itrenv);
    }
  }
}

void EventDisplay::showEnvelopeContent(TEveElement *envelope,
                                       const std::function<void(TEveElement *, const TString &)> &setChild)
{
  envelope->SetRnrSelfChildren(false, true);
  if (TEveShape *shape = dynamic_cast<TEveShape *>(envelope))
    shape->SetDrawFrame(false);
  for (TEveElement::List_i itr = envelope->BeginChildren(); itr != envelope->EndChildren(); itr++)
    setChild(*itr, TString((*itr)->GetElementName()));
}

void EventDisplay::showVertexBarrelLayers(TEveElement *vertexBarrel)
{
  // turn off the envelopes of the vertex barrel
  vertexBarrel->SetRnrSelfChildren(false, true);
  for (TEveElement::List_i itr = vertexBarrel->BeginChildren(); itr != vertexBarrel->EndChildren(); itr++)
  {
    TEveElement *b = *itr;
    TString childName(b->GetElementName());
    if (childName.Contains("assembly"))
      b->SetRnrSelfChildren(false, true);
  }
}

void EventDisplay::showECalBathCells(TEveElement *projectedGeometry, bool hideCryoSide)
{
  // the lists of the elements in the LAr bath created by buildCalorimeters
  TPRegexp re("ECalBarrel_vol_(\\w+)");
  TEveElement *projbarrel = projectedGeometry->FindChild(re);
  if (!projbarrel)
    return;
  if (hideCryoSide)
  {
    TPRegexp reCryo("ECAL_Cryo_side_(\\w+)");
    TEveElement *projcryoside = projbarrel->FindChild(reCryo);
    if (projcryoside)
      projcryoside->SetRnrSelf(false);
  }
  TEveElement *projbath = projbarrel->FindChild("LAr_bath [P]");
  if (!projbath)
    return;
  projbath->SetRnrSelfChildren(true, true);
  TEveElement *projpassives = projbath->FindChild("Passive elements [P]");
  if (projpassives)
    projpassives->SetRnrSelfChildren(false, false);
  TEveElement *projpcbs = projbath->FindChild("PCBs [P]");
  if (projpcbs)
    projpcbs->SetRnrSelfChildren(false, false);
}

void EventDisplay::showECalEndcapContent(TEveElement *ecalEndcap)
{
  // tree: one envelope per side -> cryostat walls and LAr bath -> active and
  // passive volumes between the blades of the turbine, support tubes of the
  // wheels and mechanical supports
  TPRegexp rebath("LAr_bath(\\w+)");
  // the active volumes of a wheel, between the blades, all have the same
  // projection on each side of the beam axis: draw only one active volume for
  // each projection, identified by its bounding box (rounded to 1 mm); the
  // others are not drawn, and are pruned from the geometry cache.
  // The layers inside the active volumes (calibration layers) are not drawn:
  // the readout cells are defined by the segmentation in rho and z (see
  // drawRhoZReadout)
  std::set<std::array<long, 4>> projections;
  for (TEveElement::List_i itr = ecalEndcap->BeginChildren(); itr != ecalEndcap->EndChildren(); itr++)
  {
    TEveElement *envelope = *itr;
    TEveElement *bath = envelope->FindChild(rebath);
    if (!bath)
      continue;
    // do not draw the envelope and the bath, which would cover their content
    envelope->SetRnrSelfChildren(false, true);
    bath->SetRnrSelfChildren(false, true);
    for (TEveElement::List_i itr2 = bath->BeginChildren(); itr2 != bath->EndChildren(); itr2++)
    {
      TEveElement *el = *itr2;
      TString name(el->GetElementName());
      if (name.BeginsWith("active_"))
      {
        bool newProjection = true;
        if (TEvePolygonSetProjected *projected = dynamic_cast<TEvePolygonSetProjected *>(el))
        {
          projected->ComputeBBox();
          if (const Float_t *bbox = projected->GetBBox())
            newProjection = projections.insert({std::lround(bbox[0] * 10), std::lround(bbox[1] * 10),
                                                std::lround(bbox[2] * 10), std::lround(bbox[3] * 10)}).second;
        }
        el->SetRnrSelfChildren(newProjection, false);
      }
      else
        el->SetRnrSelfChildren(name.BeginsWith("supportTube_"), false);
    }
  }
}

void EventDisplay::showECalBarrelContent(TEveElement *ecalBarrel, bool transparentActiveLayers)
{
  TPRegexp re("ECalBarrel_vol_(\\w+)");
  TEveElement *projbarrel = ecalBarrel->FindChild(re);
  if (!projbarrel)
    return;
  projbarrel->SetMainTransparency(useTransparencies ? 60 : 0);
  // do not draw the envelope, which would cover the layers, but its content
  projbarrel->SetRnrSelfChildren(false, true);
  TPRegexp rebath("LAr_bath(\\w+)");
  TEveElement *projbath = projbarrel->FindChild(rebath);
  if (!projbath)
    return;
  projbath->SetRnrSelfChildren(false, true);
  for (TEveElement::List_i itr = projbath->BeginChildren(); itr != projbath->EndChildren(); itr++)
  {
    TEveElement *el = *itr;
    TString name(el->GetElementName());
    if (name.BeginsWith("passive_") || name.BeginsWith("PCB_"))
      el->SetRnrSelfChildren(false, false);
    else if (transparentActiveLayers && name.BeginsWith("active_"))
    {
      for (TEveElement::List_i itr2 = el->BeginChildren(); itr2 != el->EndChildren(); itr2++)
        (*itr2)->SetMainTransparency(100);
    }
  }
}

void EventDisplay::showMuonLayers(TEveElement *envelope)
{
  showEnvelopeContent(envelope, [this](TEveElement *el, const TString &) {
    el->SetRnrSelfChildren(true, false);
    if (overrideDetectorColors) el->SetMainColor(kOrange);
    el->SetMainTransparency(0);
    if (TEveShape *shape = dynamic_cast<TEveShape *>(el))
      shape->SetDrawFrame(true);
  });
}

void EventDisplay::drawRhoZReadout()
{
  // the theta readout grid
  TEveStraightLineSet *grid = new TEveStraightLineSet("ECAL theta readout");
  grid->SetLineColor(kViolet);
  for (int iTheta = 0; iTheta <= geomReader->nThetaBins; iTheta++)
  {
    double theta = geomReader->thetaMin + iTheta * geomReader->thetaGrid;
    double r1 = geomReader->rMin;
    double r2 = geomReader->rMax;
    double z1 = r1 * cos(theta) / sin(theta);
    double z2 = r2 * cos(theta) / sin(theta);
    if (z1 < geomReader->zMax && z1 > geomReader->zMin)
    {
      if (z2 > geomReader->zMax)
      {
        z2 = geomReader->zMax;
        r2 = z2 * sin(theta) / cos(theta);
      }
      if (z2 < geomReader->zMin)
      {
        z2 = geomReader->zMin;
        r2 = z2 * sin(theta) / cos(theta);
      }
      grid->AddLine(z1, r1, 0., z2, r2, 0.);
      grid->AddLine(z1, -r1, 0., z2, -r2, 0.);
    }
  }
  rhoZScene->AddElement(grid);
  readout->AddElement(grid);
  // by default we do not show the unmerged readout
  grid->SetRnrSelf(false);
  readoutsHiddenByDefault.insert(grid);

  // the merged grid
  TEveStraightLineSet *grid2 = new TEveStraightLineSet("ECAL theta readout merged");
  grid2->SetLineColor(kViolet);
  grid2->SetLineWidth(3);
  for (int iLayer = 0; iLayer < geomReader->nLayers; iLayer++)
  {
    double r1 = geomReader->r[iLayer];
    double r2 = geomReader->r[iLayer + 1];
    for (int iTheta = 0; iTheta <= geomReader->nThetaBins; iTheta++)
    {
      if (iTheta % geomReader->mergedCells_Theta[iLayer] != 0)
        continue;
      double theta = geomReader->thetaMin + iTheta * geomReader->thetaGrid;
      double z1 = r1 * cos(theta) / sin(theta);
      double z2 = r2 * cos(theta) / sin(theta);
      double r2tmp = r2;
      if (z1 < geomReader->zMax && z1 > geomReader->zMin)
      {
        if (z2 > geomReader->zMax)
        {
          z2 = geomReader->zMax;
          r2tmp = z2 * sin(theta) / cos(theta);
        }
        if (z2 < geomReader->zMin)
        {
          z2 = geomReader->zMin;
          r2tmp = z2 * sin(theta) / cos(theta);
        }
        grid2->AddLine(z1, r1, 0.,
                       z2, r2tmp, 0.);
        grid2->AddLine(z1, -r1, 0.,
                       z2, -r2tmp, 0.);
      }
    }
  }
  rhoZScene->AddElement(grid2);
  readout->AddElement(grid2);

  // the ECAL endcap readout: bins in rho and z in each wheel, on both sides
  // and above and below the beam axis
  if (doEndcaps)
  {
    TEveStraightLineSet *ecalEndcapGrid = new TEveStraightLineSet("ECAL endcap rho-z readout");
    ecalEndcapGrid->SetLineColor(kViolet);
    for (int wheel = 0; wheel < geomReader->nWheelsECalEndCap; wheel++)
    {
      const int nRho = geomReader->nRhoBinsECalEndCap[wheel];
      const int nZ = geomReader->nZBinsECalEndCap[wheel];
      const double gridRho = geomReader->gridSizeRhoECalEndCap[wheel];
      const double gridZ = geomReader->gridSizeZECalEndCap[wheel];
      const double rhoMin = geomReader->offsetRhoECalEndCap[wheel];
      const double rhoMax = rhoMin + nRho * gridRho;
      const double zMin = geomReader->offsetZECalEndCap[wheel];
      const double zMax = zMin + nZ * gridZ;
      for (double side : {-1., 1.})
      {
        for (double sign : {-1., 1.})
        {
          for (int i = 0; i <= nRho; i++)
          {
            const double rho = sign * (rhoMin + i * gridRho);
            ecalEndcapGrid->AddLine(side * zMin, rho, 0., side * zMax, rho, 0.);
          }
          for (int i = 0; i <= nZ; i++)
          {
            const double z = side * (zMin + i * gridZ);
            ecalEndcapGrid->AddLine(z, sign * rhoMin, 0., z, sign * rhoMax, 0.);
          }
        }
      }
    }
    rhoZScene->AddElement(ecalEndcapGrid);
    readout->AddElement(ecalEndcapGrid);
  }

  // the muon tagger readout (with the full detector, where it is drawn): the
  // edges of the theta bins, in the barrel and in the endcaps
  if (showFullDetector)
  {
    TEveStraightLineSet *muonGrid = new TEveStraightLineSet("Muon tagger theta readout");
    muonGrid->SetLineColor(kViolet);
    // the segment of the line at theta, from the origin, with rho in
    // [rhoMin, rhoMax] and |z| in [zMin, zMax], above and below the beam axis
    auto addSegment = [muonGrid](double theta, double rhoMin, double rhoMax, double zMin, double zMax) {
      const double s = std::sin(theta), c = std::fabs(std::cos(theta));
      // the distance from the origin along the line
      double tMin = rhoMin / s, tMax = rhoMax / s;
      if (c > 0.)
      {
        tMin = std::max(tMin, zMin / c);
        tMax = std::min(tMax, zMax / c);
      }
      else if (zMin > 0.)
        return;
      if (tMin >= tMax)
        return;
      for (double sign : {-1., 1.})
        muonGrid->AddLine(tMin * std::cos(theta), sign * tMin * s, 0., tMax * std::cos(theta), sign * tMax * s, 0.);
    };
    for (int i = 0;; i++)
    {
      const double theta = geomReader->thetaOffsetMuon + (i - 0.5) * geomReader->thetaGridMuon;
      if (theta >= TMath::Pi())
        break;
      if (theta <= 0.)
        continue;
      addSegment(theta, geomReader->rMinMuonBarrel, geomReader->rMaxMuonBarrel, 0., geomReader->zMaxMuonBarrel);
      addSegment(theta, geomReader->rMinMuonEndcap, geomReader->rMaxMuonEndcap,
                 geomReader->zMinMuonEndcap, geomReader->zMaxMuonEndcap);
    }
    rhoZScene->AddElement(muonGrid);
    readout->AddElement(muonGrid);
  }

  // draw the HCAL readout segmentation in eta or theta (rho-z view)
  if (doHCal)
  {
    TEveStraightLineSet *hcalRhoZReadout = new TEveStraightLineSet("HCAL theta readout");
    hcalRhoZReadout->SetLineColor(kViolet);
    hcalRhoZReadout->SetLineWidth(5);
    for (int iTheta = 0; iTheta <= geomReader->nThetaBinsHCal; iTheta++)
    {
      double theta = geomReader->thetaMinHCal + iTheta * geomReader->thetaGridHCal;

      double r1 = geomReader->rMinHCal;
      double r2 = geomReader->rMaxHCal;
      double z1 = r1 * cos(theta) / sin(theta);
      double z2 = r2 * cos(theta) / sin(theta);
      if (z1 < geomReader->zMaxHCal && z1 > geomReader->zMinHCal)
      {
        if (z2 > geomReader->zMaxHCal)
        {
          z2 = geomReader->zMaxHCal;
          r2 = z2 * sin(theta) / cos(theta);
        }
        if (z2 < geomReader->zMinHCal)
        {
          z2 = geomReader->zMinHCal;
          r2 = z2 * sin(theta) / cos(theta);
        }
        hcalRhoZReadout->AddLine(z1, r1, 0., z2, r2, 0.);
        hcalRhoZReadout->AddLine(z1, -r1, 0., z2, -r2, 0.);
      }
    }
    rhoZScene->AddElement(hcalRhoZReadout);
    readout->AddElement(hcalRhoZReadout);
  }
}

void EventDisplay::writeGeometryCache()
{
  // write the geometry cache, now that the volumes that can be drawn are known
  if (geometryExtract)
  {
    GeometryCache::write(geomFile, geometryCacheTag(), geometryCodeVersion, geometryExtract,
                         geometryExtracts, {geom, rhoPhiScene, rhoZScene});
    GeometryCache::deleteExtract(geometryExtract, true);
    geometryExtract = nullptr;
    geometryExtracts.clear();
  }
}

void EventDisplay::openEventFile()
{
  // setup the event reader
  std::cout << std::endl;
  std::cout << "******************************************************************************" << std::endl;
  std::cout << "Setting up the event reader" << std::endl;
  std::cout << "******************************************************************************" << std::endl
            << std::endl;

  std::cout << "Reading event data from file " << evtFile << std::endl
            << std::endl;

  // open the file with the events and create the reader
  eventFile.reset(TFile::Open(evtFile.c_str(), "READ"));
  eventReader = std::make_unique<EventReader>(eventFile.get());

  // setup the branches
  eventReader->setBranches();

  // decode the cellIDs with the CellIDEncoding stored in the file
  setCellIDEncodings();

  // a set of clusters for each collection of clusters, with the size of the
  // window of the sliding-window clusters
  // (with several collections of a type, their names are added to the type)
  auto addClusterSets = [this](const std::vector<EventReader::Clusters> &clustersList, const std::string &type,
                               const std::string &label, Color_t centerColor, bool slidingWindow)
  {
    for (const EventReader::Clusters &clusters : clustersList)
    {
      const std::string suffix = clustersList.size() > 1 ? " " + clusters.name : "";
      ClusterSet set{type + suffix, label + suffix, centerColor, slidingWindow};
      set.clusters = &clusters;
      if (slidingWindow)
        setSlidingWindow(set);
      clusterSets.push_back(std::move(set));
    }
  };
  addClusterSets(eventReader->topoClusters, "topo", "topo", kWhite, false);
  addClusterSets(eventReader->swClusters, "sw", "SW", kGray, true);

  // print updated draw settings based on info found in ROOT file
  displayConfig.print();

  // read the number of events in the file
  nEvents = eventReader->nEvents;
  updateEventEntries();
}

void EventDisplay::addPaletteLegends()
{
  // vertical legend of the cluster cell energies in the projections, on the
  // right, inside the axis labels (palette values in MeV, shown in GeV),
  // shown at startup if set in the configuration (not without events)
  showPaletteLegends = displayConfig.getBoolConfig("showCellEnergyLegend") && evtFile != "";
  const std::vector<std::pair<TGLViewer *, TEveRGBAPalette *>> legends = {
      {rhoPhiGLView, cellPaletteRhoPhi}, {rhoZGLView, cellPaletteRhoZ}};
  for (const auto &[viewer, palette] : legends)
  {
    if (!viewer)
      continue;
    PaletteLegend *legend = new PaletteLegend(palette, 0.82, 0.30, 0.02, 0.45, 0.001);
    viewer->AddOverlayElement(legend);
    TGLConstAnnotation *title = new TGLConstAnnotation(viewer, "E [GeV]", 0.81, 0.80);
    title->SetTextSize(0.03);
    title->SetAllowClose(false);
    paletteLegends.push_back(legend);
    paletteLegends.push_back(title);
  }
  for (TGLOverlayElement *el : paletteLegends)
    el->SetBinaryState(showPaletteLegends);
  updateToggleButtons();
}

void EventDisplay::setInitialCamera()
{
  // Set the 3D view as the active tab and rotate the camera
  browser->GetTabRight()->SetTab(0);

  // rotate camera
  mainGLView->CurrentCamera().RotateRad(-0.5, -2.4);
  mainGLView->DoDraw();

  // clipping
  mainGLView->GetClipSet()->SetClipType(TGLClip::EType(2));
  double data[6];
  mainGLView->GetClipSet()->GetClipState(TGLClip::EType(2), data);
  // set center of bounding box at z=0 and deltaZ to full world size
  data[2]=0.;
  data[5] = data[5]*2;
  mainGLView->GetClipSet()->SetClipState(TGLClip::EType(2), data);
  mainGLView->RefreshPadEditor(mainGLView);
  mainGLView->DoDraw();

  activeGLViewer = mainGLView;
}

/******************************************************************************/
// GUI
/******************************************************************************/

void EventDisplay::fwd()
{
  if (evtFile == "")
  {
    printf("\nNo events loaded\n");
    int ret;
    new TGMsgBox(gClient->GetRoot(), browser,
                          "Error", "No events loaded",
                          kMBIconStop, kMBOk, &ret);
  }
  else if (eventId < nEvents - 1)
  {
    ++eventId;
    drawEvent(eventId);
  }
  else
  {
    printf("\nAlready at last event\n");

    int ret;
    new TGMsgBox(gClient->GetRoot(), browser,
                          "Error", "Already at last event",
                          kMBIconStop, kMBOk, &ret);
  }
}

void EventDisplay::bck()
{
  if (evtFile == "")
  {
    printf("\nNo events loaded\n");
    int ret;
    new TGMsgBox(gClient->GetRoot(), browser,
                          "Error", "No events loaded",
                          kMBIconStop, kMBOk, &ret);
  }
  else if (eventId > 0)
  {
    --eventId;
    drawEvent(eventId);
  }
  else
  {
    printf("\nAlready at first event\n");

    int ret;
    new TGMsgBox(gClient->GetRoot(), browser,
                          "Error", "Already at first event",
                          kMBIconStop, kMBOk, &ret);
  }
}

void EventDisplay::goToEvent()
{
  if (evtFile == "" || eventEntries.empty())
    return;

  // the entry whose text field emitted the signal (one per tab)
  auto it = std::find_if(eventEntries.begin(), eventEntries.end(),
                         [](TGNumberEntry *e) { return e->GetNumberEntry() == gTQSender; });
  TGNumberEntry *entry = (it != eventEntries.end()) ? *it : eventEntries.front();
  Long_t event = entry->GetIntNumber();

  if (event < 0 || event >= nEvents)
  {
    printf("\nEvent %ld out of range [0, %d]\n", event, nEvents - 1);
    int ret;
    new TGMsgBox(gClient->GetRoot(), browser,
                 "Error", Form("Event %ld out of range [0, %d]", event, nEvents - 1),
                 kMBIconStop, kMBOk, &ret);
    updateEventEntries();
    return;
  }
  if (event == eventId)
    return;
  drawEvent(event);
}

void EventDisplay::stepEvent(Long_t val)
{
  // TGNumberEntry::ValueChanged: val / 10000 != 0 for the down arrow
  bool up = (val / 10000 == 0);

  // the arrows are repeat-fire buttons, which keep firing until the mouse
  // button is released. If fwd/bck is going to open a (modal) message box,
  // release them first: otherwise the release goes to the message box, and
  // the arrow keeps firing, opening a new message box each time
  if (evtFile == "" || (up && eventId >= nEvents - 1) || (!up && eventId <= 0))
  {
    Event_t release{};
    release.fType = kButtonRelease;
    for (TGNumberEntry *entry : eventEntries)
    {
      entry->GetButtonUp()->HandleButton(&release);
      entry->GetButtonDown()->HandleButton(&release);
    }
  }

  if (up)
    fwd();
  else
    bck();
}

void EventDisplay::updateEventEntries()
{
  for (TGNumberEntry *entry : eventEntries)
    entry->SetIntNumber(eventId, kFALSE);
  for (TGLabel *label : eventRangeLabels)
  {
    label->SetText(nEvents > 0 ? Form("(0-%d)", nEvents - 1) : "(no events)");
    label->Resize(label->GetDefaultSize());
    if (label->GetParent())
      ((TGCompositeFrame *) label->GetParent())->Layout();
  }
}

void EventDisplay::toggleReadout()
{
  if (showReadout) {
    showReadout = false;
  }
  else {
    showReadout = true;
  }
  readout->SetRnrSelfChildren(showReadout, showReadout);
  for (TEveElement::List_i itr = readout->BeginChildren(); itr != readout->EndChildren(); itr++) {
    TEveElement *el = *itr;
    el->SetRnrSelf(showReadout && readoutsHiddenByDefault.count(el) == 0);
  }
  updateToggleButtons();
  // update now the scenes with the changed elements (faster than a full
  // redraw, or a repaint of the projected geometry scenes, with the full
  // detector), then redraw the views of the projections, where the readouts
  // are, which are not redrawn by it
  gEve->DoRedraw3D();
  for (TGLViewer *viewer : {rhoPhiGLView, rhoZGLView})
    if (viewer)
      viewer->RequestDraw(TGLRnrCtx::kLODHigh);
}

void EventDisplay::togglePaletteLegends()
{
  showPaletteLegends = !showPaletteLegends;
  for (TGLOverlayElement *el : paletteLegends)
    el->SetBinaryState(showPaletteLegends);
  updateToggleButtons();
  // (at high quality: the default one may not draw all the volumes)
  for (TGLViewer *viewer : {rhoPhiGLView, rhoZGLView})
    if (viewer)
      viewer->RequestDraw(TGLRnrCtx::kLODHigh);
}

void EventDisplay::takeScreenshot()
{
  TString view;
  if (activeGLViewer == rhoPhiGLView)
    view = "rhophi";
  else if (activeGLViewer == rhoZGLView)
    view = "rhoz";
  else if (activeGLViewer == mainGLView)
    view = "3d";
  else
    std::cout << "Unknown view!!" << std::endl;

  TTimeStamp ts;
  TString tss(ts.AsString("s"));
  tss.ReplaceAll(" ", "_");

  if (!std::filesystem::exists("screenshots/")) {
    std::filesystem::create_directory("screenshots/");
  }

  TString filename = "screenshots/calodisplay_screenshot_" + view + "_" + tss;
  std::cout << "Saving a screenshot of the selected view to " << filename << " ... ";
  activeGLViewer->SavePictureScale((filename+".png").Data(), 10.);
  activeGLViewer->SavePicture((filename+".pdf").Data());
  std::cout << "DONE" << std::endl << std::endl;
}

// Create minimal GUI for event navigation.
void EventDisplay::makeGui()
{
  // add buttons to go through events
  TGTab* tab = browser->GetTabRight();
  Int_t nTabs = tab->GetNumberOfTabs();
  for (Int_t i = 0; i < nTabs; ++i) {

    TGCompositeFrame* frm = tab->GetTabContainer(i);

    // (the buttons and the event number entry are expanded to the height of
    // the bar, i.e. of the tallest of them, the buttons with the icons)
    TGHorizontalFrame* hf = new TGHorizontalFrame(frm);
    // (its buttons, labels and layouts are deleted with it)
    hf->SetCleanup(kDeepCleanup);
    buttonBars.push_back(hf);
    // picture button with the icon from the installation directory, or a text
    // button if the icon is not found
    auto makeButton = [hf](const char* icon, const char* label) -> TGButton* {
      const TGPicture* picture = gClient->GetPicture(PathUtils::findFile(std::string("icons/") + icon).c_str());
      if (picture)
        return new TGPictureButton(hf, picture);
      return new TGTextButton(hf, label);
    };
    TGButton* b = nullptr;

    b = makeButton("TakeScreenshot.png", "Screenshot");
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 10, 2, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "takeScreenshot()");
    b->SetToolTipText("Take a screenshot of the selected view");

    // event number entry: type the number and press Enter to go to that
    // event, or use its up/down arrows to go to the next/previous event
    // (no limits, so that out-of-range numbers can be reported)
    TGNumberEntry *entry = new TGNumberEntry(hf, 0, 6, -1, TGNumberFormat::kNESInteger,
                                             TGNumberFormat::kNEAAnyNumber);
    // the arrows do not change the number themselves, they emit ValueChanged
    entry->SetButtonToNum(kFALSE);
    entry->Connect("ValueChanged(Long_t)", "EventDisplay", this, "stepEvent(Long_t)");
    entry->GetNumberEntry()->Connect("ReturnPressed()", "EventDisplay", this, "goToEvent()");
    entry->GetNumberEntry()->SetToolTipText("Event number: type it and press Enter to go to that event");
    entry->GetButtonUp()->SetToolTipText("Go to next event");
    entry->GetButtonDown()->SetToolTipText("Go to previous event");
    entry->SetState(evtFile != "");
    hf->AddFrame(entry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 10, 2, 10, 10));
    eventEntries.push_back(entry);
    TGLabel *rangeLabel = new TGLabel(hf, "(no events)");
    hf->AddFrame(rangeLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 10, 10, 10));
    eventRangeLabels.push_back(rangeLabel);

    b = makeButton("grid_icon.png", "Readout");
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 10, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "toggleReadout()");
    b->SetToolTipText("Hide/Show readouts");
    b->AllowStayDown(kTRUE);
    readoutButtons.push_back(b);

    b = makeButton("palette_icon.png", "Legend");
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 10, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "togglePaletteLegends()");
    b->SetToolTipText("Hide/Show the legend of the cluster cell energies in the projections");
    b->AllowStayDown(kTRUE);
    // (no legend without events)
    if (evtFile == "")
      b->SetState(kButtonDisabled);
    paletteLegendButtons.push_back(b);

    b = new TGTextButton(hf, "Camera");
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 10, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "printCamera()");
    b->SetToolTipText("Print the camera of the view, as a --camera option for the screenshots");

    frm->AddFrame(hf);
    frm->Layout();
    frm->MapSubwindows();
  }

  updateToggleButtons();

  // hide command tab
  gEve->GetBrowser()->HideBottomTab();
}

void EventDisplay::onCloseWindow()
{
  for (TGCompositeFrame *bar : buttonBars)
  {
    static_cast<TGCompositeFrame *>(const_cast<TGWindow *>(bar->GetParent()))->RemoveFrame(bar);
    bar->DestroyWindow();
    delete bar;
  }
  buttonBars.clear();
  eventEntries.clear();
  eventRangeLabels.clear();
  readoutButtons.clear();
  paletteLegendButtons.clear();
  // (after TEveManager::Terminate, called next, back in the event loop;
  // ending the line of the ROOT prompt, for the one of the shell)
  std::cout << std::endl;
  TTimer::SingleShot(0, "TApplication", gApplication, "Terminate(=0)");
}

void EventDisplay::updateToggleButtons()
{
  // in all the tabs, also the ones of the buttons not clicked: a pressed
  // stay-down button is engaged (not just down, or its next click would not
  // be signalled); the buttons already in the right state are not changed, as
  // the clicked one, which is still handling its click, and the disabled ones
  auto setPressed = [](const std::vector<TGButton *> &buttons, bool pressed) {
    for (TGButton *b : buttons)
      if (b->GetState() != kButtonDisabled && b->IsDown() != pressed)
        b->SetState(pressed ? kButtonEngaged : kButtonUp);
  };
  setPressed(readoutButtons, showReadout);
  setPressed(paletteLegendButtons, showPaletteLegends);
}

void EventDisplay::onTabSelected(Int_t tab)
{
  if (tab==0) {
    activeGLViewer = mainGLView;
  }
  else if (tab==1) {
    activeGLViewer = rhoPhiGLView;
    if (!initRhoPhiView) {
      buildRhoPhiProjection();
      rhoPhiView->Redraw(true);
      initRhoPhiView = true;
    }
  }
  else if (tab==2) {
    activeGLViewer = rhoZGLView;
    if (!initRhoZView) {
      buildRhoZProjection();
      rhoZView->Redraw(true);
      initRhoZView = true;
    }
  }
}


/******************************************************************************/
// other helpers
/******************************************************************************/

std::vector<TEvePathMarkD> EventDisplay::getElectronBremsstrahlungMarks(
    unsigned int electronIndex,
    const TTreeReaderArray<Int_t>& daughters,
    const TTreeReaderArray<Int_t>& pdg,
    const TTreeReaderArray<Double_t>& vx,
    const TTreeReaderArray<Double_t>& vy,
    const TTreeReaderArray<Double_t>& vz,
    const TTreeReaderArray<Double_t>& px,
    const TTreeReaderArray<Double_t>& py,
    const TTreeReaderArray<Double_t>& pz,
    const TTreeReaderArray<UInt_t>& daughters_begin,
    const TTreeReaderArray<UInt_t>& daughters_end)
{
  std::vector<TEvePathMarkD> marks;
  if(electronIndex >= pdg.GetSize()) return marks;
  if(std::abs(pdg[electronIndex]) != 11)
    return marks;

  //
  // Initial electron momentum and vertex
  //
  TVector3 electronMomentum(
                            px[electronIndex],
                            py[electronIndex],
                            pz[electronIndex]);

  TVector3 electronVertex(
                          vx[electronIndex],
                          vy[electronIndex],
                          vz[electronIndex]);

  //
  // Find daughters
  //
  std::vector<unsigned int> photons;
  bool hasElectronDaughter = false;
  for(unsigned int i=daughters_begin[electronIndex];
      i<daughters_end[electronIndex];
      ++i)
  {
    if(i>=daughters.GetSize()) continue;

    int child = daughters[i];
    if(child < 0) continue;

    unsigned int childIndex = static_cast<unsigned int>(child);
    if(childIndex >= pdg.GetSize()) continue;

    if(std::abs(pdg[childIndex]) == 11)
      hasElectronDaughter = true;

    if(pdg[childIndex] == 22)
      photons.push_back(childIndex);
  }


  //
  // If Geant produced a new electron particle,
  // the daughter electron already has the correct
  // momentum history.
  //
  if(hasElectronDaughter)
    return marks;

  if(photons.empty())
    return marks;

  //
  // Sort photons along electron flight direction
  //
  TVector3 direction = electronMomentum.Unit();

  std::sort(photons.begin(), photons.end(),
            [&](unsigned int a, unsigned int b)
            {
              TVector3 va(vx[a], vy[a], vz[a]);
              TVector3 vb(vx[b], vy[b], vz[b]);
              double da = (va-electronVertex).Dot(direction);
              double db = (vb-electronVertex).Dot(direction);
              return da < db;
            });

  //
  // Follow the electron momentum after each photon
  //
  TVector3 currentMomentum = electronMomentum;
  for(unsigned int gamma : photons)
  {
    if(gamma >= static_cast<unsigned int>(px.GetSize()))
      continue;
    TVector3 photonMomentum(px[gamma], py[gamma], pz[gamma]);

    //
    // Bremsstrahlung:
    // p_e(new) = p_e(old) - p_gamma
    //
    currentMomentum -= photonMomentum;
    double p = currentMomentum.Mag();
    if(p==0)
      continue;

    TEveVectorD position(vx[gamma]*mm, vy[gamma]*mm, vz[gamma]*mm);
    TEveVectorD momentum(currentMomentum.X(), currentMomentum.Y(), currentMomentum.Z());
    marks.emplace_back(TEvePathMarkD::kReference,position,momentum);
  }
  return marks;
}

std::vector<ElectronSegment> EventDisplay::getElectronBremsstrahlungSegments(
    unsigned int electronIndex,
    const TTreeReaderArray<Int_t>& daughters,
    const TTreeReaderArray<Int_t>& pdg,
    const TTreeReaderArray<Double_t>& vx,
    const TTreeReaderArray<Double_t>& vy,
    const TTreeReaderArray<Double_t>& vz,
    const TTreeReaderArray<Double_t>& px,
    const TTreeReaderArray<Double_t>& py,
    const TTreeReaderArray<Double_t>& pz,
    const TTreeReaderArray<UInt_t>& daughters_begin,
    const TTreeReaderArray<UInt_t>& daughters_end)
{
  // 1. Pre-calculate search bounds and safely reserve space
  const size_t d_begin = daughters_begin[electronIndex];
  const size_t d_end = daughters_end[electronIndex];
  const size_t max_daughters = (d_end > d_begin) ? (d_end - d_begin) : 0;

  // Use a pair to store: <Distance Along Direction, Photon Index>
  std::vector<std::pair<double, unsigned int>> photonsWithDist;
  photonsWithDist.reserve(max_daughters);

  // Initial momentum and vertex
  const double el_px = px[electronIndex];
  const double el_py = py[electronIndex];
  const double el_pz = pz[electronIndex];
  const double el_vx = vx[electronIndex]*mm;
  const double el_vy = vy[electronIndex]*mm;
  const double el_vz = vz[electronIndex]*mm;

  TVector3 electronMomentum(el_px, el_py, el_pz);
  TVector3 direction = electronMomentum.Unit();

  // Find photon daughters and pre-calculate projection distances
  const size_t daughters_size = daughters.GetSize();
  for(unsigned int i=d_begin; i<d_end; ++i)
  {
    if(i >= daughters_size) continue;

    int child = daughters[i];
    if(child < 0) continue;

    if(pdg[child] == 22)
    {
      // Pre-calculate projection along direction to avoid heavy lifting inside std::sort
      double dx = vx[child]*mm - el_vx;
      double dy = vy[child]*mm - el_vy;
      double dz = vz[child]*mm - el_vz;
      double dist = dx * direction.X() + dy * direction.Y() + dz * direction.Z();
      photonsWithDist.emplace_back(dist, static_cast<unsigned int>(child));
    }
  }

  //
  // Sort photons by distance along electron direction
  // (rather sort by time...)
  std::sort(photonsWithDist.begin(), photonsWithDist.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

  // Create electron segments
  std::vector<ElectronSegment> segments;
  segments.reserve(photonsWithDist.size() + 1);

  TVector3 currentMomentum = electronMomentum;
  TVector3 currentVertex(el_vx, el_vy, el_vz);
  int charge = - MathUtils::sign(pdg[electronIndex]);

  // Create propagator (local, only used to extrapolate the electron between
  // the photon vertices; it must not own the field, which is shared)
  TEveTrackPropagator propagator;
  propagator.SetMagFieldObj(magField.get(), false);
  // propagator.SetStepper(TEveTrackPropagator::kHelix);

  for(const auto& [dist, gamma] : photonsWithDist)
  {
    const double g_vx = vx[gamma] * mm;
    const double g_vy = vy[gamma] * mm;
    const double g_vz = vz[gamma] * mm;

    TVector3 start(currentVertex.X(), currentVertex.Y(), currentVertex.Z());
    TVector3 end(g_vx, g_vy, g_vz);

    segments.push_back({start, end, currentMomentum});

    // Extrapolate electron to gamma vertex
    TEveVectorD start_pos(start.X(), start.Y(), start.Z());
    propagator.InitTrack(start_pos, charge);
    TEveVectorD target_pos(g_vx, g_vy, g_vz);

    TEveVectorD propagated_mom(currentMomentum.X(), currentMomentum.Y(), currentMomentum.Z());
    Bool_t success = propagator.GoToVertex(target_pos, propagated_mom);
    if (not success) {
      std::cout << "ERROR: failed to propagate electron!" << std::endl;
    }

    // Update electron momentum after photon emission
    currentMomentum.SetXYZ(propagated_mom.fX - px[gamma],
                           propagated_mom.fY - py[gamma],
                           propagated_mom.fZ - pz[gamma]);

    // Set new vertex
    currentVertex.SetXYZ(g_vx, g_vy, g_vz);
  }

  // Add final electron segment
  TVector3 endpoint((*eventReader->genParticles.endpoint_x)[electronIndex]*mm,
                    (*eventReader->genParticles.endpoint_y)[electronIndex]*mm,
                    (*eventReader->genParticles.endpoint_z)[electronIndex]*mm);
  segments.push_back({
      TVector3(currentVertex.X(), currentVertex.Y(), currentVertex.Z()),
      TVector3(endpoint.X(), endpoint.Y(), endpoint.Z()),
      currentMomentum
    });
  return segments;
}

std::string EventDisplay::geometryCacheTag() const
{
  return Form("fulldet%d_hcal%d_endcap%d", showFullDetector, doHCal, doEndcaps);
}

void EventDisplay::saveScreenshots(const std::string &prefix, int firstEvent, int nEventsToSave,
                                   const std::vector<std::string> &views, int width)
{
  // the tab and the GL viewer of each view
  const std::map<std::string, std::pair<int, TGLViewer *>> viewers = {
      {"3d", {0, mainGLView}}, {"rhophi", {1, rhoPhiGLView}}, {"rhoz", {2, rhoZGLView}}};
  // (without event file, only the geometry)
  const int lastEvent = eventReader ? std::min(firstEvent + nEventsToSave, nEvents) : firstEvent + 1;
  for (int event = firstEvent; event < lastEvent; event++)
  {
    if (eventReader && event != eventId)
      drawEvent(event);
    for (const std::string &view : views)
    {
      auto it = viewers.find(view);
      if (it == viewers.end())
      {
        std::cout << "WARNING: unknown view " << view << " (3d, rhophi or rhoz), no screenshot" << std::endl;
        continue;
      }
      // show the view (which projects the geometry the first time), let it be
      // drawn and frame it
      browser->GetTabRight()->SetTab(it->second.first);
      gSystem->ProcessEvents();
      setCamera(view);
      const std::string file = Form("%s_ev%d_%s.png", prefix.c_str(), event, view.c_str());
      if (it->second.second->SavePictureWidth(file.c_str(), width))
      {
        std::cout << "Screenshot saved to " << file << ", with camera ";
        printCamera();
      }
      else
        std::cout << "WARNING: cannot save the screenshot " << file << std::endl;
    }
  }
}

bool EventDisplay::getFocusBox(const std::string &view, double boxMin[3], double boxMax[3])
{
  // the 3D points of the objects to frame
  std::vector<std::array<double, 3>> points;
  if (cameraFocus == "cluster" || cameraFocus == "track+cluster")
  {
    // the cells of the most energetic cluster (or its position, if they are
    // not available)
    const ClusterSet *bestSet = nullptr;
    int best = -1;
    for (const ClusterSet &set : clusterSets)
      for (unsigned int i = 0; i < set.clusters->size(); i++)
        if (best < 0 || (*set.clusters->energy)[i] > (*bestSet->clusters->energy)[best])
        {
          bestSet = &set;
          best = i;
        }
    if (best >= 0)
    {
      const EventReader::Clusters &clusters = *bestSet->clusters;
      const EventReader::ClusterCells &cells = clusters.cells;
      bool cellFound = false;
      for (unsigned int i = (*clusters.hits_begin)[best]; i < (*clusters.hits_end)[best] && i < cells.size(); i++)
        if (cells.cellID[i] != 0)
        {
          points.push_back({cells.x[i] * mm, cells.y[i] * mm, cells.z[i] * mm});
          cellFound = true;
        }
      if (!cellFound)
      {
        const double size = 25. * cm;
        const double x = (*clusters.x)[best] * mm, y = (*clusters.y)[best] * mm, z = (*clusters.z)[best] * mm;
        points.push_back({x - size, y - size, z - size});
        points.push_back({x + size, y + size, z + size});
      }
    }
  }
  if (cameraFocus == "track" || cameraFocus == "track+cluster")
  {
    // the trajectory of the most energetic charged particle, from its vertex
    // (the electrons are drawn as segments between the bremsstrahlung photons:
    // all the tracks of the same type with at least half of its momentum)
    std::vector<TEveTrack *> charged;
    TEveTrack *best = nullptr;
    if (particles)
      for (TEveElement::List_i it = particles->BeginChildren(); it != particles->EndChildren(); ++it)
      {
        TEveTrack *track = dynamic_cast<TEveTrack *>(*it);
        if (!track || track->GetCharge() == 0)
          continue;
        charged.push_back(track);
        if (!best || track->GetMomentum().Mag() > best->GetMomentum().Mag())
          best = track;
      }
    for (TEveTrack *track : charged)
    {
      if (track->GetPdg() != best->GetPdg() || track->GetMomentum().Mag() < 0.5 * best->GetMomentum().Mag())
        continue;
      const TEveVectorD &vertex = track->GetVertex();
      points.push_back({vertex.fX, vertex.fY, vertex.fZ});
      for (int i = 0; i < track->Size(); i++)
      {
        float x, y, z;
        track->GetPoint(i, x, y, z);
        points.push_back({x, y, z});
      }
    }
  }
  if (points.empty())
    return false;

  // the box in the coordinates of the view
  for (int j = 0; j < 3; j++)
  {
    boxMin[j] = std::numeric_limits<double>::max();
    boxMax[j] = -std::numeric_limits<double>::max();
  }
  for (const auto &p : points)
  {
    double q[3] = {p[0], p[1], p[2]};
    if (view == "rhophi")
      q[2] = 0.;
    else if (view == "rhoz")
    {
      q[0] = p[2];
      q[1] = MathUtils::sign(p[1]) * std::sqrt(p[0] * p[0] + p[1] * p[1]);
      q[2] = 0.;
    }
    for (int j = 0; j < 3; j++)
    {
      boxMin[j] = std::min(boxMin[j], q[j]);
      boxMax[j] = std::max(boxMax[j], q[j]);
    }
  }
  return true;
}

void EventDisplay::setOrthoCamera(TGLViewer *viewer, double x, double y, double zoom)
{
  TGLOrthoCamera &camera = dynamic_cast<TGLOrthoCamera &>(viewer->CurrentCamera());
  camera.Reset();
  camera.SetZoom(zoom);
  // move the centre of the view (in the plane of the camera) from the default one
  const TGLVector3 center = (camera.RefCamBase() * camera.RefCamTrans()).GetTranslation();
  camera.Truck(x - center.X(), y - center.Y());
  camera.IncTimeStamp();
}

void EventDisplay::setPerspectiveCamera(TGLViewer *viewer, const double center[3], const double eye[3],
                                        const double up[3], double fov)
{
  TGLPerspectiveCamera &camera = dynamic_cast<TGLPerspectiveCamera &>(viewer->CurrentCamera());
  camera.SetFOV(fov);
  camera.SetExternalCenter(kTRUE);
  camera.SetCenterVec(center[0], center[1], center[2]);
  // the camera frame: axis 1 from the centre to the eye, axis 3 vertical,
  // at the position of the eye
  TGLVector3 axis1(eye[0] - center[0], eye[1] - center[1], eye[2] - center[2]);
  axis1.Normalise();
  TGLVector3 axis3(up[0], up[1], up[2]);
  axis3 = axis3 - Dot(axis3, axis1) * axis1;
  axis3.Normalise();
  TGLMatrix frame;
  frame.SetBaseVec(1, axis1);
  frame.SetBaseVec(2, Cross(axis3, axis1));
  frame.SetBaseVec(3, axis3);
  frame.SetBaseVec(4, eye[0], eye[1], eye[2]);
  TGLMatrix baseInverse = camera.RefCamBase();
  baseInverse.Invert();
  camera.RefCamTrans() = baseInverse * frame;
  camera.IncTimeStamp();
}

void EventDisplay::setCamera(const std::string &view)
{
  TGLViewer *viewer = view == "3d" ? mainGLView : view == "rhophi" ? rhoPhiGLView : rhoZGLView;

  // set explicitly
  auto it = cameraSettings.find(view);
  if (it != cameraSettings.end())
  {
    const std::vector<double> &v = it->second;
    if (view == "3d" && v.size() == 10)
      setPerspectiveCamera(viewer, &v[1], &v[4], &v[7], v[0]);
    else if (view != "3d" && v.size() == 3)
      setOrthoCamera(viewer, v[0], v[1], v[2]);
    else
      std::cout << "WARNING: wrong number of values in --camera " << view << " (10 for 3d, 3 for the projections)" << std::endl;
    return;
  }

  // automatic, only for the projections (the 3D view shows the full
  // detector): the box of the object to frame, scaled, fills the image
  double boxMin[3], boxMax[3];
  if (view == "3d" || !getFocusBox(view, boxMin, boxMax))
    return;
  const double scale = cameraFocusScale > 0. ? cameraFocusScale : (cameraFocus == "cluster" ? 2. : 1.2);
  double center[2], size[2];
  for (int j = 0; j < 2; j++)
  {
    center[j] = 0.5 * (boxMin[j] + boxMax[j]);
    // (at least 10 cm, e.g. for a cluster with a single cell)
    size[j] = std::max(scale * (boxMax[j] - boxMin[j]), 10. * cm);
  }
  // the zoom at which the box fits, from the region shown by the viewport at
  // the current zoom
  TGLOrthoCamera &camera = dynamic_cast<TGLOrthoCamera &>(viewer->CurrentCamera());
  const TGLRect &viewport = camera.RefViewport();
  const TGLVector3 shown = camera.ViewportDeltaToWorld(TGLVertex3(0., 0., 0.), viewport.Width(), viewport.Height());
  const double zoom = camera.GetZoom() * std::min(std::fabs(shown.X()) / size[0], std::fabs(shown.Y()) / size[1]);
  setOrthoCamera(viewer, center[0], center[1], zoom);
}

void EventDisplay::printCamera()
{
  if (activeGLViewer == mainGLView)
  {
    TGLPerspectiveCamera &camera = dynamic_cast<TGLPerspectiveCamera &>(mainGLView->CurrentCamera());
    const TGLMatrix frame = camera.RefCamBase() * camera.RefCamTrans();
    const TGLVector3 eye = frame.GetTranslation();
    const TGLVector3 axis1 = frame.GetBaseVec(1);
    const TGLVector3 axis3 = frame.GetBaseVec(3);
    // the centre on the axis of the camera, at the distance of the rotation centre
    const Double_t *c = camera.GetCenterVec();
    const double distance = std::sqrt((eye.X() - c[0]) * (eye.X() - c[0]) + (eye.Y() - c[1]) * (eye.Y() - c[1]) +
                                      (eye.Z() - c[2]) * (eye.Z() - c[2]));
    std::cout << Form("--camera 3d=%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g", camera.GetFOV(),
                      eye.X() - distance * axis1.X(), eye.Y() - distance * axis1.Y(), eye.Z() - distance * axis1.Z(),
                      eye.X(), eye.Y(), eye.Z(), axis3.X(), axis3.Y(), axis3.Z())
              << std::endl;
  }
  else
  {
    TGLOrthoCamera &camera = dynamic_cast<TGLOrthoCamera &>(activeGLViewer->CurrentCamera());
    const TGLVector3 center = (camera.RefCamBase() * camera.RefCamTrans()).GetTranslation();
    std::cout << Form("--camera %s=%.6g,%.6g,%.6g", activeGLViewer == rhoPhiGLView ? "rhophi" : "rhoz",
                      center.X(), center.Y(), camera.GetZoom())
              << std::endl;
  }
}

void EventDisplay::attachClusterProjections()
{
  // the projections of the clusters are drawn by hand in the projected
  // coordinates, so they are not in the event, which is projected
  // automatically: they are added to the projected replicas of the 3D lists of
  // the clusters, to be shown or hidden with them (and with the event)
  for (ClusterSet &set : clusterSets)
  {
    if (set.list3D == nullptr)
      continue;
    for (auto it = set.list3D->BeginProjecteds(); it != set.list3D->EndProjecteds(); ++it)
    {
      TEveElement *replica = (*it)->GetProjectedAsElement();
      if ((*it)->GetManager() == rhoZProjManager)
        replica->AddElement(set.listRhoZ);
      else if ((*it)->GetManager() == rhoPhiProjManager)
        replica->AddElement(set.listRhoPhi);
    }
  }
}

void EventDisplay::createElementList(TEveElementList*& list, const std::string& listName)
{
  if (list == nullptr)
  {
    list = new TEveElementList(listName.c_str());
    gEve->AddElement(list);
  }
}

void EventDisplay::setSlidingWindow(ClusterSet &set)
{
  SlidingWindow &swWindow = set.window;
  swWindow = SlidingWindow();
  const std::string &collection = set.clusters->name;
  std::map<std::string, std::string> options = eventReader->readGaudiOptions();
  if (options.empty())
  {
    std::cout << "WARNING: no job options in the event file, cannot draw the windows of the sliding-window clusters" << std::endl;
    return;
  }

  // name of the algorithm (or tool) with a property <property> = <value>
  auto findOwner = [&options](const std::string &property, const std::string &value) -> std::string {
    const std::string suffix = "." + property;
    for (const auto &[key, val] : options)
      if (val == value && key.size() > suffix.size() &&
          key.compare(key.size() - suffix.size(), suffix.size(), suffix) == 0)
        return key.substr(0, key.size() - suffix.size());
    return "";
  };
  auto getOption = [&options](const std::string &key) -> std::string {
    auto it = options.find(key);
    return it == options.end() ? "" : it->second;
  };

  // go back from the clusters in the event file (e.g. after calibration) to
  // the ones created by the sliding-window algorithm
  std::string clusters = collection;
  std::string algorithm = findOwner("clusters", clusters);
  for (int i = 0; i < 10 && (algorithm == "" || getOption(algorithm + ".nThetaFinal") == ""); i++)
  {
    std::string previous = findOwner("outClusters", clusters);
    if (previous == "")
      break;
    clusters = getOption(previous + ".inClusters");
    algorithm = findOwner("clusters", clusters);
  }
  std::string towerTool = getOption(algorithm + ".towerTool");
  // the tool is given as type/name
  towerTool = towerTool.substr(towerTool.find('/') + 1);
  const std::string tool = "ToolSvc." + towerTool + ".";
  try
  {
    swWindow.nTheta = std::stoi(getOption(algorithm + ".nThetaFinal"));
    swWindow.nPhi = std::stoi(getOption(algorithm + ".nPhiFinal"));
    swWindow.deltaTheta = std::stod(getOption(tool + "deltaThetaTower"));
    swWindow.deltaPhi = std::stod(getOption(tool + "deltaPhiTower"));
    swWindow.thetaMax = std::stod(getOption(tool + "thetaMax"));
    swWindow.phiMin = std::stod(getOption(tool + "phiMin"));
  }
  catch (const std::exception &)
  {
    std::cout << "WARNING: cannot find the window size of the sliding-window clusters " << collection
              << " in the job options of the event file, not drawing the windows" << std::endl;
    swWindow = SlidingWindow();
    return;
  }

  // radial range of the barrel calorimeters and |z| range of the endcap
  // calorimeters used for the towers (the calorimeterIDs are listed as [ 4 , 8 ])
  std::string ids = getOption(tool + "calorimeterIDs");
  std::replace_if(ids.begin(), ids.end(), [](char c) { return !std::isdigit(static_cast<unsigned char>(c)); }, ' ');
  std::istringstream idStream(ids);
  int id;
  while (idStream >> id)
  {
    double rMin = 0., rMax = 0.;
    if (id == DetectorGeometry::kECalBarrel)
    {
      rMin = geomReader->rMin;
      rMax = geomReader->rMax;
    }
    else if (id == DetectorGeometry::kHCalBarrel)
    {
      rMin = geomReader->rMinHCal;
      rMax = geomReader->rMaxHCal;
    }
    else if (id == DetectorGeometry::kECalEndcap)
    {
      // the active part of the ECAL endcap, i.e. the z bins of its readout
      for (int wheel = 0; wheel < geomReader->nWheelsECalEndCap; wheel++)
      {
        const double zMin = geomReader->offsetZECalEndCap[wheel];
        const double zMax = zMin + geomReader->nZBinsECalEndCap[wheel] * geomReader->gridSizeZECalEndCap[wheel];
        swWindow.zMin = (swWindow.zMin == 0.) ? zMin : std::min(swWindow.zMin, zMin);
        swWindow.zMax = std::max(swWindow.zMax, zMax);
      }
      continue;
    }
    else
      continue;
    swWindow.rMin = (swWindow.rMin == 0.) ? rMin : std::min(swWindow.rMin, rMin);
    swWindow.rMax = std::max(swWindow.rMax, rMax);
  }
  if (swWindow.rMax == 0. && swWindow.zMax == 0.)
  {
    std::cout << "WARNING: the sliding-window clusters " << collection
              << " do not use the barrel calorimeters or the ECAL endcap, not drawing their windows" << std::endl;
    swWindow = SlidingWindow();
    return;
  }
  std::cout << "Window of the sliding-window clusters " << collection << " (from " << algorithm << "): "
            << swWindow.nTheta << " x " << swWindow.nPhi << " towers of "
            << swWindow.deltaTheta << " x " << swWindow.deltaPhi << " rad in theta x phi" << std::endl;
}

void EventDisplay::getSlidingWindowRange(const SlidingWindow &swWindow, float xcl, float ycl, float zcl,
                                         double &thetaLow, double &thetaHigh,
                                         double &phiLow, double &phiHigh) const
{
  const double theta = std::atan2(std::sqrt(xcl * xcl + ycl * ycl), zcl);
  const double phi = std::atan2(ycl, xcl);
  // towers are numbered from thetaMax down, and from phiMin up
  const int idTheta = (int)std::floor((swWindow.thetaMax - theta) / swWindow.deltaTheta);
  const int idPhi = (int)std::floor((phi - swWindow.phiMin) / swWindow.deltaPhi);
  thetaHigh = swWindow.thetaMax - (idTheta - swWindow.nTheta / 2) * swWindow.deltaTheta;
  thetaLow = swWindow.thetaMax - (idTheta + swWindow.nTheta / 2 + 1) * swWindow.deltaTheta;
  phiLow = swWindow.phiMin + (idPhi - swWindow.nPhi / 2) * swWindow.deltaPhi;
  phiHigh = swWindow.phiMin + (idPhi + swWindow.nPhi / 2 + 1) * swWindow.deltaPhi;
}

void EventDisplay::drawSlidingWindow(const SlidingWindow &swWindow, float xcl, float ycl, float zcl, const char *name,
                                     TEveElement *rhoz, TEveElement *rhophi)
{
  if (swWindow.nTheta == 0)
    return;
  // the window is drawn in the barrel or in the endcap, depending on the
  // position of the cluster and on the calorimeters used for the towers
  const bool inBarrel = std::fabs(zcl * mm) <= geomReader->zMax;
  if ((inBarrel && swWindow.rMax == 0.) || (!inBarrel && swWindow.zMax == 0.))
    return;
  double thetaLow, thetaHigh, phiLow, phiHigh;
  getSlidingWindowRange(swWindow, xcl, ycl, zcl, thetaLow, thetaHigh, phiLow, phiHigh);

  // rho-z: between the lines at constant theta, on the side of the cluster
  const double sign = MathUtils::sign(ycl);
  TEveStraightLineSet *lines = new TEveStraightLineSet(name);
  lines->SetLineColor(kYellow);
  lines->SetLineWidth(2);
  // on top of the projected detector, which is in the same plane
  lines->SetDepthTest(kFALSE);
  if (!inBarrel)
  {
    // endcap: between the inner and outer |z|, on the side of the cluster
    // (rho = z tan(theta) on both sides); not drawn in rho-phi, as the cells
    // of the endcap
    const double zIn = MathUtils::sign(zcl) * swWindow.zMin;
    const double zOut = MathUtils::sign(zcl) * swWindow.zMax;
    lines->AddLine(zIn, sign * zIn * std::tan(thetaLow), 0., zOut, sign * zOut * std::tan(thetaLow), 0.);
    lines->AddLine(zOut, sign * zOut * std::tan(thetaLow), 0., zOut, sign * zOut * std::tan(thetaHigh), 0.);
    lines->AddLine(zOut, sign * zOut * std::tan(thetaHigh), 0., zIn, sign * zIn * std::tan(thetaHigh), 0.);
    lines->AddLine(zIn, sign * zIn * std::tan(thetaHigh), 0., zIn, sign * zIn * std::tan(thetaLow), 0.);
    rhoz->AddElement(lines);
    return;
  }

  // barrel: between the inner and outer radius
  const double rIn = swWindow.rMin;
  const double rOut = swWindow.rMax;
  lines->AddLine(rIn / std::tan(thetaLow), sign * rIn, 0., rOut / std::tan(thetaLow), sign * rOut, 0.);
  lines->AddLine(rOut / std::tan(thetaLow), sign * rOut, 0., rOut / std::tan(thetaHigh), sign * rOut, 0.);
  lines->AddLine(rOut / std::tan(thetaHigh), sign * rOut, 0., rIn / std::tan(thetaHigh), sign * rIn, 0.);
  lines->AddLine(rIn / std::tan(thetaHigh), sign * rIn, 0., rIn / std::tan(thetaLow), sign * rIn, 0.);
  rhoz->AddElement(lines);

  // rho-phi: sector of the ring between the inner and outer radius
  lines = new TEveStraightLineSet(name);
  lines->SetLineColor(kYellow);
  lines->SetLineWidth(2);
  // on top of the projected detector, which is in the same plane
  lines->SetDepthTest(kFALSE);
  const int nSegments = 16;
  for (int i = 0; i < nSegments; i++)
  {
    const double phi1 = phiLow + (phiHigh - phiLow) * i / nSegments;
    const double phi2 = phiLow + (phiHigh - phiLow) * (i + 1) / nSegments;
    lines->AddLine(rIn * std::cos(phi1), rIn * std::sin(phi1), 0., rIn * std::cos(phi2), rIn * std::sin(phi2), 0.);
    lines->AddLine(rOut * std::cos(phi1), rOut * std::sin(phi1), 0., rOut * std::cos(phi2), rOut * std::sin(phi2), 0.);
  }
  for (double phi : {phiLow, phiHigh})
    lines->AddLine(rIn * std::cos(phi), rIn * std::sin(phi), 0., rOut * std::cos(phi), rOut * std::sin(phi), 0.);
  rhophi->AddElement(lines);
}

void EventDisplay::setCellIDEncodings()
{
  std::map<std::string, std::string> encodings = eventReader->readCellIDEncodings();
  if (encodings.empty())
  {
    std::cout << "WARNING: no CellIDEncoding in the metadata of the event file, using the default cellID encodings" << std::endl;
    return;
  }

  // the encoding of each calorimeter readout is taken from its cell collection
  // (also when its cells are not drawn, as they can be in the clusters)
  const std::vector<std::tuple<std::string, std::string, DetectorGeometry::Readout>> readouts = {
    {"ecalBarrelCells", "drawECalBarrelCells", DetectorGeometry::Readout::ECalBarrel},
    {"ecalEndcapCells", "drawECalEndcapCells", DetectorGeometry::Readout::ECalEndcap},
    {"hcalBarrelCells", "drawHCalBarrelCells", DetectorGeometry::Readout::HCalBarrel},
    {"hcalEndcapCells", "drawHCalEndcapCells", DetectorGeometry::Readout::HCalEndcap},
    {"muonBarrelCells", "drawMuonCells", DetectorGeometry::Readout::MuonBarrel},
    {"muonEndcapCells", "drawMuonCells", DetectorGeometry::Readout::MuonEndcap}};
  for (const auto& [key, drawFlag, readoutType] : readouts)
  {
    std::string collection = displayConfig.getStringConfig(key);
    if (collection == "")
      continue;
    auto it = encodings.find(collection);
    // the encoding is also stored for the readoutType, i.e. without "Positioned"
    const std::string positioned = "Positioned";
    if (it == encodings.end() && collection.size() > positioned.size() &&
        collection.compare(collection.size() - positioned.size(), positioned.size(), positioned) == 0)
      it = encodings.find(collection.substr(0, collection.size() - positioned.size()));
    if (it == encodings.end())
    {
      // (the detector might not be in the file: warn only if its cells are drawn)
      if (displayConfig.getBoolConfig(drawFlag))
        std::cout << "WARNING: no CellIDEncoding for " << key << " (" << collection
                  << ") in the metadata of the event file, using the default encoding" << std::endl;
      continue;
    }
    DetectorGeometry::SetCellIDEncoding(readoutType, it->second);
  }
}

EventDisplay::SubDetector* EventDisplay::findSubDetector(TString const& s)
{
  for (auto& g : m_subDetectors)
  {
    if (!g.regex.empty()) {
      if (TPRegexp(g.regex.c_str()).MatchB(s))
        return &g;
    }
    for (auto const& c : g.contains)
      if (s.Contains(c))
        return &g;
  }
  std::cout << "Unexpected volume: " << s << std::endl;
  return nullptr;
}

TEveElementList *EventDisplay::subDetectorList(SubDetector &det)
{
  if (! det.list) {
    det.list = new TEveElementList(det.name.c_str());
    geom->AddElement(det.list);
    if (! det.visibleByDefault)
      det.list->SetRnrSelfChildren(false, false);
  }
  return det.list;
}

void EventDisplay::setVolumeStyle(TEveElement *el, Color_t color, int nSegments)
{
  el->SetRnrSelfChildren(true, false);
  el->SetMainTransparency(useTransparencies ? 60 : 0);
  if (overrideDetectorColors) el->SetMainColor(color);
  if (nSegments > 0)
    ((TEveGeoShape *)el)->SetNSegments(nSegments);
  ((TEveGeoShape *)el)->SetDrawFrame(false);
}

void EventDisplay::addSimpleElement(TEveElement* el, SubDetector& det)
{
  if (debug) {
    std::cout << "Adding " << el->GetElementName() << " to " << det.name << std::endl;
  }
  subDetectorList(det)->AddElement(el);
  ((TEveGeoShape *)el)->SetDrawFrame(false);
}

void EventDisplay::addVertex(TEveElement* el, SubDetector& /*det*/)
{
  TEveElementList* barrel = new TEveElementList("Vertex barrel");
  geom->AddElement(barrel);
  TEveElementList* endcap = new TEveElementList("Vertex endcap");
  geom->AddElement(endcap);

  // for vertex, loop over children to add them either to barrel or endcap
  TString s(el->GetElementName());
  if (debug) {
    std::cout << "Found " << s << ", looping over its children" << std::endl;
  }
  for (TEveElement::List_i itr = el->BeginChildren(); itr != el->EndChildren(); ++itr)
  {
    TEveElement* child = *itr;
    TString name(child->GetElementName());
    if (name.Contains("VertexBarrel") || name.Contains("VertexInnerBarrel"))
    {
      if (debug) std::cout << "Adding " << name << " to vertex barrel" << std::endl;
      // barrel->AddElement(child);
      // child->SetMainColor(kRed);
      // rather than the assembly we add directly the layers
      for (TEveElement::List_i itr2 = child->BeginChildren(); itr2 != child->EndChildren(); itr2++) {
        TEveElement *grandchild = *itr2;
        TString gcName(grandchild->GetElementName());
        barrel->AddElement(grandchild);
        if (overrideDetectorColors) grandchild->SetMainColor(kRed);
        ((TEveGeoShape *)grandchild)->SetDrawFrame(false);
        if (gcName.Contains("VertexInnerBarrel"))
          ((TEveGeoShape *)grandchild)->SetNSegments(128);
        // by default draw only the sensors
        if (!gcName.Contains("sensor")) {
          grandchild->SetRnrSelfChildren(false, false);
        }
        else {
          // sensors are so thin that some of them are not drawn.
          // tried this but it does not fix the problem..
          //((TEveGeoShape *) grandchild)->SetMiniFrame(false);
          //((TEveGeoShape *) grandchild)->SetDrawFrame(true);
          //grandchild->SetRnrSelfChildren(true, false);
        }
      }
    }
    else if (name.Contains("VertexDisks"))
    {
      if (debug) std::cout << "Adding " << name << " to vertex endcap" << std::endl;
      endcap->AddElement(child);
      if (overrideDetectorColors) child->SetMainColor(kRed);
      ((TEveGeoShape *)child)->SetDrawFrame(false);
      // by default draw only the sensors
      if (!name.Contains("sensor"))
        child->SetRnrSelfChildren(false, false);
      else {
        // sensors are so thin that some of them are not drawn.
        // tried this but it does not fix the problem..
        // ((TEveGeoShape *) child)->SetDrawFrame(true);
        // ((TEveGeoShape *) child)->SetMiniFrame(false);
      }
    }
  }
}

void EventDisplay::addDCH(TEveElement* el, SubDetector& det)
{
  addSimpleElement(el, det);
  setVolumeStyle(el, kGreen-5, 128);
}

void EventDisplay::addWrapper(TEveElement* el, SubDetector& det)
{
  // rather than the assembly we add directly the layers
  TEveElementList *list = subDetectorList(det);
  if (debug) {
    std::cout << "Adding " << el->GetElementName() << " to " << det.name << std::endl;
  }
  for (TEveElement::List_i itr2 = el->BeginChildren(); itr2 != el->EndChildren(); itr2++) {
    TEveElement *child = *itr2;
    list->AddElement(child);
    setVolumeStyle(child, kRed, 128);
  }
}

void EventDisplay::addECalB(TEveElement* el, SubDetector& det)
{
  addSimpleElement(el, det);
  setVolumeStyle(el, kAzure+7, 0);

  // daughter properties - essentially for projections...
  for (TEveElement::List_i itr2 = el->BeginChildren(); itr2 != el->EndChildren(); itr2++)
  {
    TEveElement *child = *itr2;
    TString elname(child->GetElementName());
    if (debug) std::cout << elname << std::endl;
    setVolumeStyle(child, kAzure+7, 128);
    if (elname.BeginsWith("ECAL_Cryo_side"))
      child->SetRnrSelfChildren(false, false);
  }
}

void EventDisplay::addECalEC(TEveElement* el, SubDetector& det)
{
  addSimpleElement(el, det);
  setVolumeStyle(el, kAzure+7, 0);
}

void EventDisplay::addHCalB(TEveElement* el, SubDetector& det)
{
  // for the HCal barrel, we draw the envelope in the 3d model,
  // and the subvolumes in the 2D views (see projectGeometryRhoPhi/RhoZ)
  addSimpleElement(el, det);
  setVolumeStyle(el, kAzure-7, 128);

  // setup the properties of the children of the envelope (for the projection)
  // for NSegments I do it here since it does not seem to work in the projection
  for (TEveElement::List_i itr = el->BeginChildren(); itr != el->EndChildren(); itr++)
  {
    TEveElement *child = *itr;
    TString elname(child->GetElementName());
    if (overrideDetectorColors and elname.Contains("LayerVol")) {
      child->SetMainColor(kAzure-7);
    }
    // show the radial layers
    ((TEveGeoShape *)child)->SetDrawFrame(true);
    ((TEveGeoShape *)child)->SetNSegments(128);
  }
}

void EventDisplay::addHCalEC(TEveElement* el, SubDetector& det)
{
  // for the HCal endcap, we draw the envelope in the 3d model,
  // and the subvolumes in the 2D views (see projectGeometryRhoZ)
  addSimpleElement(el, det);
  setVolumeStyle(el, kAzure-7, 128);

  // setup the properties of the children of the envelope
  for (TEveElement::List_i itr = el->BeginChildren(); itr != el->EndChildren(); itr++)
  {
    TEveElement *child = *itr;
    TString elname(child->GetElementName());
    if (elname.Contains("ECLayer")) {
      // show the layers
      child->SetRnrSelfChildren(false, true);
      child->SetMainColor(kAzure-7);
      child->SetMainTransparency(0);
      // set to false if we don't want to see the layer structure
      ((TEveGeoShape *)child)->SetDrawFrame(true);
    }
  }
}

void EventDisplay::addMuonB(TEveElement* el, SubDetector& det)
{
  // draw envelope but not children
  addSimpleElement(el, det);
  setVolumeStyle(el, kOrange, 128);

  // setup the properties of the children of the envelope
  // for the projections (number of segments for rho-phi view)
  for (TEveElement::List_i itr = el->BeginChildren(); itr != el->EndChildren(); itr++)
  {
    TEveElement *child = *itr;
    ((TEveGeoShape *)child)->SetNSegments(128);
    ((TEveGeoShape *)child)->SetDrawFrame(true);
  }
}

void EventDisplay::addMuonEC(TEveElement* el, SubDetector& det)
{
  // draw envelope but not children (the layers are shown in the projections)
  addSimpleElement(el, det);
  setVolumeStyle(el, kOrange, 0);
}
