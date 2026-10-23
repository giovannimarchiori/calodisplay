/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
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

// ROOT core
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTimeStamp.h>

// TG
#include <TGButton.h>
#include <TGTab.h>
#include <TGMsgBox.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>

// TGeo
#include <TGeoTube.h>
#include <TGeoMatrix.h>

// TEve
#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEveProjectionAxes.h>
#include <TEveGeoNode.h>


#include <algorithm>
#include <filesystem>
#include <map>
#include <tuple>
#include <unordered_map>
#include <iostream>

using Units::units;
using Units::cm;
using Units::mm;

namespace
{
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


void EventDisplay::fillClusters(std::string clusterType)
{
  // the cells of the clusters of the current event
  const EventReader::ClusterCells &cells = (clusterType == "topo") ? eventReader->topoClusterCells : eventReader->swClusterCells;
  // read once, not in the loop over the clusters
  const bool calculateBaryCenters = displayConfig.getBoolConfig("calculateClusterBaryCenters");

  std::vector<std::unique_ptr<CaloCluster>> *clusterData = nullptr;
  unsigned int nClusters = 0;

  if (clusterType == "topo")
  {
    std::cout << "Creating topo clusters" << std::endl;
    clusterData = &topoclusterData;
    nClusters = eventReader->CaloTopoClusters_position_x->GetSize();
  }
  else if (clusterType == "sw")
  {
    std::cout << "Creating SW clusters" << std::endl;
    clusterData = &swclusterData;
    nClusters = eventReader->CaloClusters_position_x->GetSize();
  }
  else
  {
    std::cout << "Unknown cluster type " << clusterType << std::endl;
    return;
  }

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

  // loop over clusters and fill the relevant info
  const int nECalBarrelLayers = geomReader->nLayers;
  const int nECalEndCapLayers = geomReader->nLayersECalEndCap;
  const int nHCalBarrelLayers = geomReader->nLayersHCal;
  const int nHCalEndCapLayers = geomReader->nLayersHCalEndCap;
  const int nMuonLayers = geomReader->nLayersMuon;
  for (unsigned int i = 0; i < nClusters; i++)
  {
    // keep only those above threshold
    float E = (clusterType == "topo") ? (*eventReader->CaloTopoClusters_energy)[i] : (*eventReader->CaloClusters_energy)[i];
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
    if (clusterType == "topo")
    {
      cluster->setBarycenterXYZ((*eventReader->CaloTopoClusters_position_x)[i] * mm,
                                (*eventReader->CaloTopoClusters_position_y)[i] * mm,
                                (*eventReader->CaloTopoClusters_position_z)[i] * mm);
    }
    else
    {
      cluster->setBarycenterXYZ((*eventReader->CaloClusters_position_x)[i] * mm,
                                (*eventReader->CaloClusters_position_y)[i] * mm,
                                (*eventReader->CaloClusters_position_z)[i] * mm);
    }

    // set cluster barycenters layer-by-layer

    // - recalculation from cells
    if (calculateBaryCenters) {
      // retrieve index of first and last+1 cell in cell vectors
      if (debug)
      {
        std::cout << "  Retrieving cluster cells and calculating energy per layer .." << std::endl;
      }
      unsigned int first_hit, last_hit;
      // (no cells if the cells of the clusters are not available)
      const unsigned int nCells = cells.size();
      if (clusterType == "topo")
      {
        first_hit = (*eventReader->CaloTopoClusters_hits_begin)[i];
        last_hit = (*eventReader->CaloTopoClusters_hits_end)[i];
      }
      else
      {
        first_hit = (*eventReader->CaloClusters_hits_begin)[i];
        last_hit = (*eventReader->CaloClusters_hits_end)[i];
      }

      // loop over cells to calculate energy per layer
      std::vector<float> energyVsECalBarrelLayer(nECalBarrelLayers, 0.);
      std::vector<float> energyVsECalEndCapLayer(nECalEndCapLayers, 0.);
      std::vector<float> energyVsHCalBarrelLayer(nHCalBarrelLayers, 0.);
      std::vector<float> energyVsHCalEndCapLayer(nHCalEndCapLayers, 0.);
      std::vector<float> energyVsMuonLayer(nMuonLayers, 0.);
      for (unsigned int iCell = first_hit; iCell < std::min(last_hit, nCells); iCell++)
      {
        ULong_t cellID;
        float energy;
        cellID = cells.cellID[iCell];
        energy = cells.energy[iCell];
        if (cellID == 0)
          continue; // cell not available (see EventReader::fillClusterCells)
        ULong_t systemID = geomReader->SystemID(cellID);
        ULong_t layer = 0;
        std::vector<float> *energyVsLayer = nullptr;
        if (systemID == DetectorGeometry::kECalBarrel)
        {
          layer = geomReader->ECalBarrelLayer(cellID);
          energyVsLayer = &energyVsECalBarrelLayer;
        }
        else if (systemID == DetectorGeometry::kHCalBarrel)
        {
          layer = geomReader->HCalBarrelLayer(cellID);
          energyVsLayer = &energyVsHCalBarrelLayer;
        }
        else if (systemID == DetectorGeometry::kECalEndcap)
        {
          layer = geomReader->ECalEndCapZBin(cellID);
          energyVsLayer = &energyVsECalEndCapLayer;
        }
        else if (systemID == DetectorGeometry::kHCalEndcap)
        {
          layer = geomReader->HCalEndCapLayer(cellID);
          energyVsLayer = &energyVsHCalEndCapLayer;
        }
        else if (systemID == DetectorGeometry::kMuonBarrel || systemID == DetectorGeometry::kMuonEndcap)
        {
          // TODO: calculate properly muon barrel and endcap layer
          layer = 0;
          energyVsLayer = &energyVsMuonLayer;
        }
        else
        {
          std::cout << "WARNING: unknown system " << systemID << " for cell " << cellID << ", skipping cell" << std::endl;
          continue;
        }
        if (layer >= energyVsLayer->size())
        {
          std::cout << "WARNING: layer " << layer << " of system " << systemID << " out of range [0, " << energyVsLayer->size() << "), skipping cell" << std::endl;
          continue;
        }
        (*energyVsLayer)[layer] += energy;
      }

      cluster->setEnergyVsECalBarrelLayers(energyVsECalBarrelLayer);
      cluster->setEnergyVsECalEndCapLayers(energyVsECalEndCapLayer);
      cluster->setEnergyVsHCalBarrelLayers(energyVsHCalBarrelLayer);
      cluster->setEnergyVsHCalEndCapLayers(energyVsHCalEndCapLayer);
      cluster->setEnergyVsMuonLayers(energyVsMuonLayer);

      // loop again over cells to calculate barycenter per layer
      if (debug)
      {
        std::cout << "  Calculating cell barycenter per layer .." << std::endl;
      }
      std::vector<TVector3> barycenterVsECalBarrelLayer(nECalBarrelLayers);
      std::vector<float> sumWeightsVsECalBarrelLayer(nECalBarrelLayers, 0.0);

      std::vector<TVector3> barycenterVsECalEndCapLayer(nECalEndCapLayers);
      std::vector<float> sumWeightsVsECalEndCapLayer(nECalEndCapLayers, 0.0);

      std::vector<TVector3> barycenterVsHCalBarrelLayer(nHCalBarrelLayers);
      std::vector<float> sumWeightsVsHCalBarrelLayer(nHCalBarrelLayers, 0.0);

      std::vector<TVector3> barycenterVsHCalEndCapLayer(nHCalEndCapLayers);
      std::vector<float> sumWeightsVsHCalEndCapLayer(nHCalEndCapLayers, 0.0);

      std::vector<TVector3> barycenterVsMuonLayer(nMuonLayers);
      std::vector<float> sumWeightsVsMuonLayer(nMuonLayers, 0.0);

      for (unsigned int iCell = first_hit; iCell < std::min(last_hit, nCells); iCell++)
      {
        ULong_t cellID;
        float energy, x, y, z;
        cellID = cells.cellID[iCell];
        energy = cells.energy[iCell];
        x = cells.x[iCell] * mm;
        y = cells.y[iCell] * mm;
        z = cells.z[iCell] * mm;
        ULong_t systemID = geomReader->SystemID(cellID);
        ULong_t layer = 0;
        std::vector<float> *energyVsLayer = nullptr;
        std::vector<TVector3> *barycenterVsLayer = nullptr;
        std::vector<float> *sumWeightsVsLayer = nullptr;
        if (systemID == DetectorGeometry::kECalBarrel)
        {
          layer = geomReader->ECalBarrelLayer(cellID);
          energyVsLayer = &energyVsECalBarrelLayer;
          barycenterVsLayer = &barycenterVsECalBarrelLayer;
          sumWeightsVsLayer = &sumWeightsVsECalBarrelLayer;
        }
        else if (systemID == DetectorGeometry::kECalEndcap)
        {
          layer = geomReader->ECalEndCapZBin(cellID);
          energyVsLayer = &energyVsECalEndCapLayer;
          barycenterVsLayer = &barycenterVsECalEndCapLayer;
          sumWeightsVsLayer = &sumWeightsVsECalEndCapLayer;
        }
        else if (systemID == DetectorGeometry::kHCalBarrel)
        {
          layer = geomReader->HCalBarrelLayer(cellID);
          energyVsLayer = &energyVsHCalBarrelLayer;
          barycenterVsLayer = &barycenterVsHCalBarrelLayer;
          sumWeightsVsLayer = &sumWeightsVsHCalBarrelLayer;
        }
        else if (systemID == DetectorGeometry::kHCalEndcap)
        {
          layer = geomReader->HCalEndCapLayer(cellID);
          energyVsLayer = &energyVsHCalEndCapLayer;
          barycenterVsLayer = &barycenterVsHCalEndCapLayer;
          sumWeightsVsLayer = &sumWeightsVsHCalEndCapLayer;
        }
        // other systems (muon, unknown) are not used for the barycenter;
        // out-of-range layers were already reported in the first loop
        if (energyVsLayer == nullptr || layer >= energyVsLayer->size())
          continue;
        if ((*energyVsLayer)[layer] > 0)
        {
          float weight = energy / (*energyVsLayer)[layer];
          TVector3 &barycenter = (*barycenterVsLayer)[layer];
          barycenter.SetXYZ(barycenter.X() + x * weight,
                            barycenter.Y() + y * weight,
                            barycenter.Z() + z * weight);
          (*sumWeightsVsLayer)[layer] += weight;
        }
      }

      for (int iLayer = 0; iLayer < nECalBarrelLayers; iLayer++)
      {
        if (energyVsECalBarrelLayer[iLayer] > 0)
        {
          barycenterVsECalBarrelLayer[iLayer].SetXYZ(
            barycenterVsECalBarrelLayer[iLayer].X() / sumWeightsVsECalBarrelLayer[iLayer],
            barycenterVsECalBarrelLayer[iLayer].Y() / sumWeightsVsECalBarrelLayer[iLayer],
            barycenterVsECalBarrelLayer[iLayer].Z() / sumWeightsVsECalBarrelLayer[iLayer]);
        }
      }
      cluster->setBarycenterVsECalBarrelLayers(barycenterVsECalBarrelLayer);

      for (int iLayer = 0; iLayer < nECalEndCapLayers; iLayer++)
      {
        if (energyVsECalEndCapLayer[iLayer] > 0)
        {
          barycenterVsECalEndCapLayer[iLayer].SetXYZ(
            barycenterVsECalEndCapLayer[iLayer].X() / sumWeightsVsECalEndCapLayer[iLayer],
            barycenterVsECalEndCapLayer[iLayer].Y() / sumWeightsVsECalEndCapLayer[iLayer],
            barycenterVsECalEndCapLayer[iLayer].Z() / sumWeightsVsECalEndCapLayer[iLayer]);
        }
      }
      cluster->setBarycenterVsECalEndCapLayers(barycenterVsECalEndCapLayer);

      for (int iLayer = 0; iLayer < nHCalBarrelLayers; iLayer++)
      {
        if (energyVsHCalBarrelLayer[iLayer] > 0.)
        {
          barycenterVsHCalBarrelLayer[iLayer].SetXYZ(
            barycenterVsHCalBarrelLayer[iLayer].X() / sumWeightsVsHCalBarrelLayer[iLayer],
            barycenterVsHCalBarrelLayer[iLayer].Y() / sumWeightsVsHCalBarrelLayer[iLayer],
            barycenterVsHCalBarrelLayer[iLayer].Z() / sumWeightsVsHCalBarrelLayer[iLayer]);
        }
      }
      cluster->setBarycenterVsHCalBarrelLayers(barycenterVsHCalBarrelLayer);

      for (int iLayer = 0; iLayer < nHCalEndCapLayers; iLayer++)
      {
        if (energyVsHCalEndCapLayer[iLayer] > 0.)
        {
          barycenterVsHCalEndCapLayer[iLayer].SetXYZ(
            barycenterVsHCalEndCapLayer[iLayer].X() / sumWeightsVsHCalEndCapLayer[iLayer],
            barycenterVsHCalEndCapLayer[iLayer].Y() / sumWeightsVsHCalEndCapLayer[iLayer],
            barycenterVsHCalEndCapLayer[iLayer].Z() / sumWeightsVsHCalEndCapLayer[iLayer]);
        }
      }
      cluster->setBarycenterVsHCalEndCapLayers(barycenterVsHCalEndCapLayer);

      cluster->setBarycenterVsMuonLayers(barycenterVsMuonLayer);
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
      trkProp->SetMagFieldObj(magField);
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
    for (unsigned int ip = 0; ip < eventReader->genParticles_generatorStatus->GetSize(); ip++)
    {
      float m = (*eventReader->genParticles_mass)[ip];
      float px = (*eventReader->genParticles_momentum_x)[ip];
      float py = (*eventReader->genParticles_momentum_y)[ip];
      float pz = (*eventReader->genParticles_momentum_z)[ip];
      float p = sqrt(px * px + py * py + pz * pz);
      if (p < ParticleEnergyThreshold)
        continue;

      double t = (*eventReader->genParticles_time)[ip];
      double x1 = (*eventReader->genParticles_vertex_x)[ip] * mm;
      double y1 = (*eventReader->genParticles_vertex_y)[ip] * mm;
      double z1 = (*eventReader->genParticles_vertex_z)[ip] * mm;

      TEveMCTrack mct;
      mct.SetPdgCode((*eventReader->genParticles_PDG)[ip]);
      mct.SetMomentum(px, py, pz, sqrt(p * p + m * m));
      mct.SetProductionVertex(x1, y1, z1, t);

      TEveVectorF v;
      v[0] = (*eventReader->genParticles_endpoint_x)[ip] * mm;
      v[1] = (*eventReader->genParticles_endpoint_y)[ip] * mm;
      v[2] = (*eventReader->genParticles_endpoint_z)[ip] * mm;

      // for electrons, for which there could be some hard bremsstrahlung, we
      // split the trajectory into pieces between photon emissions
      if (std::abs((*eventReader->genParticles_PDG)[ip]) == 11)
      {
        // retrieve the various splittings
        auto segments = getElectronBremsstrahlungSegments(
                                                          ip,
                                                          *eventReader->genParticles_daughters_index,
                                                          *eventReader->genParticles_PDG,
                                                          *eventReader->genParticles_vertex_x,
                                                          *eventReader->genParticles_vertex_y,
                                                          *eventReader->genParticles_vertex_z,
                                                          *eventReader->genParticles_momentum_x,
                                                          *eventReader->genParticles_momentum_y,
                                                          *eventReader->genParticles_momentum_z,
                                                          *eventReader->genParticles_daughters_begin,
                                                          *eventReader->genParticles_daughters_end);
        // create a track for each segment
        for(const auto& seg : segments)
        {
          TEveMCTrack mctSegment;
          mctSegment.SetPdgCode((*eventReader->genParticles_PDG)[ip]);
          double pSeg = seg.momentum.Mag();
          double pxSeg = seg.momentum.X();
          double pySeg = seg.momentum.Y();
          double pzSeg = seg.momentum.Z();
          double xSeg = seg.start.X();
          double ySeg = seg.start.Y();
          double zSeg = seg.start.Z();
          double mSeg = (*eventReader->genParticles_mass)[ip];
          mctSegment.SetMomentum(pxSeg, pySeg, pzSeg, sqrt(pSeg*pSeg+mSeg*mSeg));
          mctSegment.SetProductionVertex(xSeg, ySeg, zSeg, 0);
          TEveTrackPropagator* prop = new TEveTrackPropagator("electron_brem_segment",
                                                              "electron brem segment",
                                                              magField,
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
        TEveTrackPropagator* prop = new TEveTrackPropagator(Form("propagator_%d", ip), "per-track propagator", magField, false);
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
      for (unsigned int i = 0; i < eventReader->SimParticleSecondaries_PDG->GetSize(); i++)
      {
        pdgID = (*eventReader->SimParticleSecondaries_PDG)[i];
        // keep only photons
        if (pdgID != 22)
          continue;
        // if ( (* eventReader->SimParticleSecondaries_generatorStatus)[i] != 1 ) continue;
        px = (*eventReader->SimParticleSecondaries_momentum_x)[i];
        py = (*eventReader->SimParticleSecondaries_momentum_y)[i];
        pz = (*eventReader->SimParticleSecondaries_momentum_z)[i];
        m = (*eventReader->SimParticleSecondaries_mass)[i];
        p = sqrt(px * px + py * py + pz * pz);
        float e = sqrt(p * p + m * m);
        // std::cout << "p = "<< p << std::endl;
        if (p < ParticleEnergyThreshold)
          continue;

        // std::cout << "PDG = "<< pdgID << std::endl;
        t = (*eventReader->SimParticleSecondaries_time)[i];
        x1 = (*eventReader->SimParticleSecondaries_vertex_x)[i] * mm;
        y1 = (*eventReader->SimParticleSecondaries_vertex_y)[i] * mm;
        z1 = (*eventReader->SimParticleSecondaries_vertex_z)[i] * mm;
        double r1 = sqrt(x1 * x1 + y1 * y1);
        // the two photons from a pi0 in the origin must come from small R
        // TODO IMPROVE FOR DDSIM FILES - WHERE PARTICLES PARENTS ARE STORED!
        if (r1 > 1.)
          continue;
        double sintheta = sqrt(px * px + py * py) / p;
        x2 = x1 + px / p * rMax / sintheta;
        y2 = y1 + py / p * rMax / sintheta;
        z2 = z1 + pz / p * rMax / sintheta;
        // double x2 = (* eventReader->SimParticleSecondaries_endpoint_x)[i] * mm;
        // double y2 = (* eventReader->SimParticleSecondaries_endpoint_y)[i] * mm;
        // double z2 = (* eventReader->SimParticleSecondaries_endpoint_z)[i] * mm;
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
  for (MyPointSet *hitSet : {vtxHits, dchHits, siwrHits, ecalHits, hcalHits, muonHits})
    if (hitSet)
      hitSet->Reset();

  // hits in vertex
  if (displayConfig.getBoolConfig("drawVertexHits"))
  {
    std::cout << "Creating vertex barrel hits" << std::endl;
    drawHitCollection(vtxHits, "VTX",
                      eventReader->VertexBarrelHits_position_x,
                      eventReader->VertexBarrelHits_position_y,
                      eventReader->VertexBarrelHits_position_z,
                      nullptr, -9999.);
    std::cout << "Creating vertex endcap hits" << std::endl;
    drawHitCollection(vtxHits, "VTX",
                      eventReader->VertexEndcapHits_position_x,
                      eventReader->VertexEndcapHits_position_y,
                      eventReader->VertexEndcapHits_position_z,
                      nullptr, -9999., false, false, true);
  }

  // hits in DCH/STT
  if (displayConfig.getBoolConfig("drawMainTrackerHits"))
  {
    std::cout << "Creating main tracker hits" << std::endl;
    drawHitCollection(dchHits, "DCH",  // TODO modify for STT
                      eventReader->MainTrackerHits_position_x,
                      eventReader->MainTrackerHits_position_y,
                      eventReader->MainTrackerHits_position_z,
                      nullptr, -9999.);
  }

  // hits in Si wrapper
  if (displayConfig.getBoolConfig("drawSiWrapperHits"))
  {
    std::cout << "Creating Si wrapper hits" << std::endl;
    drawHitCollection(siwrHits, "Si wrapper",
                      eventReader->SiWrapperBarrelHits_position_x,
                      eventReader->SiWrapperBarrelHits_position_y,
                      eventReader->SiWrapperBarrelHits_position_z,
                      nullptr, -9999.);
    drawHitCollection(siwrHits, "Si wrapper",
                      eventReader->SiWrapperEndcapHits_position_x,
                      eventReader->SiWrapperEndcapHits_position_y,
                      eventReader->SiWrapperEndcapHits_position_z,
                      nullptr, -9999., false, false, true);
  }

  // hits in ECAL
  if (displayConfig.getBoolConfig("drawECalBarrelHits"))
  {
    std::cout << "Creating ecal barrel hits" << std::endl;
    drawHitCollection(ecalHits, "ECAL",
                      eventReader->ECalBarrelHits_position_x,
                      eventReader->ECalBarrelHits_position_y,
                      eventReader->ECalBarrelHits_position_z,
                      eventReader->ECalBarrelHits_energy,
                      HitEnergyThreshold);
  }

  if (displayConfig.getBoolConfig("drawECalEndcapHits"))
  {
    std::cout << "Creating ecal endcap hits" << std::endl;
    drawHitCollection(ecalHits, "ECAL",
                      (eventReader->ECalEndcapHits_position_x),
                      (eventReader->ECalEndcapHits_position_y),
                      (eventReader->ECalEndcapHits_position_z),
                      (eventReader->ECalEndcapHits_energy),
                      HitEnergyThreshold, false, false, true);
  }

  // hits in HCAL
  if (displayConfig.getBoolConfig("drawHCalBarrelHits"))
  {
    std::cout << "Creating hcal barrel hits" << std::endl;
    drawHitCollection(hcalHits, "HCAL",
                      eventReader->HCalBarrelHits_position_x,
                      eventReader->HCalBarrelHits_position_y,
                      eventReader->HCalBarrelHits_position_z,
                      eventReader->HCalBarrelHits_energy,
                      HitEnergyThreshold);
  }
  if (displayConfig.getBoolConfig("drawHCalEndcapHits"))
  {
    std::cout << "Creating hcal endcap hits" << std::endl;
    drawHitCollection(hcalHits, "HCAL",
                      eventReader->HCalEndcapHits_position_x,
                      eventReader->HCalEndcapHits_position_y,
                      eventReader->HCalEndcapHits_position_z,
                      eventReader->HCalEndcapHits_energy,
                      HitEnergyThreshold, false, false, true);
  }

  // hits in muon detector
  if (displayConfig.getBoolConfig("drawMuonHits"))
  {
    std::cout << "Creating muon barrel hits" << std::endl;
    drawHitCollection(muonHits, "MUON",
                      eventReader->MuonBarrelHits_position_x,
                      eventReader->MuonBarrelHits_position_y,
                      eventReader->MuonBarrelHits_position_z,
                      nullptr, -9999.);
    std::cout << "Creating muon endcap hits" << std::endl;
    drawHitCollection(muonHits, "MUON",
                      eventReader->MuonEndcapHits_position_x,
                      eventReader->MuonEndcapHits_position_y,
                      eventReader->MuonEndcapHits_position_z,
                      nullptr, -9999., false, false, true);
  }
}

void EventDisplay::drawDigis()
{
  createElementList(digis, "digis");

  // clear the digis of the previous event: the endcap collections are appended
  // to the barrel ones, so they would accumulate if the barrel is not drawn
  for (MyPointSet *digiSet : {vtxDigis, dchDigis, siwrDigis, ecalCells, hcalCells, muonCells, cells_merged})
    if (digiSet)
      digiSet->Reset();

  // digis in VTX
  if (displayConfig.getBoolConfig("drawVertexDigis"))
  {
    std::cout << "Creating vertex barrel digis" << std::endl;
    drawHitCollection(vtxDigis, "VTX",
                      eventReader->VertexBarrelDigis_position_x,
                      eventReader->VertexBarrelDigis_position_y,
                      eventReader->VertexBarrelDigis_position_z,
                      nullptr, -9999., true);
    std::cout << "Creating vertex endcap digis" << std::endl;
    drawHitCollection(vtxDigis, "VTX",
                      eventReader->VertexEndcapDigis_position_x,
                      eventReader->VertexEndcapDigis_position_y,
                      eventReader->VertexEndcapDigis_position_z,
                      nullptr, -9999., true, false, true);
  }

  // digis in DCH/STT
  if (displayConfig.getBoolConfig("drawMainTrackerDigis"))
  {
    std::cout << "Creating main tracker digis" << std::endl;
    drawHitCollection(dchDigis, "DCH",  // TODO modify for STT
                      eventReader->MainTrackerDigis_position_x,
                      eventReader->MainTrackerDigis_position_y,
                      eventReader->MainTrackerDigis_position_z,
                      nullptr, -9999., true);
  }

  // digis in Si wrapper
  if (displayConfig.getBoolConfig("drawSiWrapperDigis"))
  {
    std::cout << "Creating Si wrapper digis" << std::endl;
    drawHitCollection(siwrDigis, "Si wrapper",
                      eventReader->SiWrapperBarrelDigis_position_x,
                      eventReader->SiWrapperBarrelDigis_position_y,
                      eventReader->SiWrapperBarrelDigis_position_z,
                      nullptr, -9999., true);
    drawHitCollection(siwrDigis, "Si wrapper",
                      eventReader->SiWrapperEndcapDigis_position_x,
                      eventReader->SiWrapperEndcapDigis_position_y,
                      eventReader->SiWrapperEndcapDigis_position_z,
                      nullptr, -9999., true, false, true);
  }

  // digis (cells) in ECAL
  if (displayConfig.getBoolConfig("drawECalBarrelCells"))
  {
    std::cout << "Creating ecal barrel cells" << std::endl;
    drawHitCollection(ecalCells, "ECAL",
                      eventReader->ECalBarrelCells_position_x,
                      eventReader->ECalBarrelCells_position_y,
                      eventReader->ECalBarrelCells_position_z,
                      eventReader->ECalBarrelCells_energy,
                      CellEnergyThreshold,
                      true, debug, false);
    // add tool tips
    if (debug) {
      for (unsigned int i = 0; i < eventReader->ECalBarrelCells_position_x->GetSize(); i++) {
        float E = (*eventReader->ECalBarrelCells_energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        ULong_t cellID = (*eventReader->ECalBarrelCells_cellID)[i];
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
                      eventReader->ECalEndcapCells_position_x,
                      eventReader->ECalEndcapCells_position_y,
                      eventReader->ECalEndcapCells_position_z,
                      eventReader->ECalEndcapCells_energy,
                      CellEnergyThreshold,
                      true, debug, true);
    // add tool tips
    if (debug) {
      for (unsigned int i = 0; i < eventReader->ECalEndcapCells_position_x->GetSize(); i++)
      {
        float E = (*eventReader->ECalEndcapCells_energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        ULong_t cellID = (*eventReader->ECalEndcapCells_cellID)[i];
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
                      eventReader->HCalBarrelCells_position_x,
                      eventReader->HCalBarrelCells_position_y,
                      eventReader->HCalBarrelCells_position_z,
                      eventReader->HCalBarrelCells_energy,
                      CellEnergyThreshold,
                      true, debug, false);
    // add tool tips
    if (debug) {
      for (unsigned int i = 0; i < eventReader->HCalBarrelCells_position_x->GetSize(); i++)
      {
        float E = (*eventReader->HCalBarrelCells_energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        ULong_t cellID = (*eventReader->HCalBarrelCells_cellID)[i];
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
                      eventReader->HCalEndcapCells_position_x,
                      eventReader->HCalEndcapCells_position_y,
                      eventReader->HCalEndcapCells_position_z,
                      eventReader->HCalEndcapCells_energy,
                      CellEnergyThreshold,
                      true, debug, true);
    // add tool tips
    if (debug) {
      for (unsigned int i = 0; i < eventReader->HCalEndcapCells_position_x->GetSize(); i++)
      {
        float E = (*eventReader->HCalEndcapCells_energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        ULong_t cellID = (*eventReader->HCalEndcapCells_cellID)[i];
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
                      eventReader->MuonBarrelCells_position_x,
                      eventReader->MuonBarrelCells_position_y,
                      eventReader->MuonBarrelCells_position_z,
                      eventReader->MuonBarrelCells_energy,
                      CellEnergyThreshold,
                      true, false);

    std::cout << "Creating muon endcap cells" << std::endl;
    drawHitCollection(muonCells, "MUON",
                      eventReader->MuonEndcapCells_position_x,
                      eventReader->MuonEndcapCells_position_y,
                      eventReader->MuonEndcapCells_position_z,
                      eventReader->MuonEndcapCells_energy,
                      CellEnergyThreshold,
                      true, false, true);
  }

  // cells merged in ECAL barrel
  if (displayConfig.getBoolConfig("drawECalBarrelMergedCells"))
  {
    std::cout << "Creating merged ecal barrel cells" << std::endl;
    drawHitCollection(cells_merged, "ECAL merged",
                      eventReader->ECalBarrelCells2_position_x,
                      eventReader->ECalBarrelCells2_position_y,
                      eventReader->ECalBarrelCells2_position_z,
                      nullptr, -9999.);
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
      trkProp->SetMagFieldObj(magField);
      // trkProp->SetMaxR(rMax);
      // trkProp->SetMaxZ(geomReader->zMaxEndCap);
      // trkProp->SetMaxZ(10000.*mm);
      trkProp->SetRnrReferences(true);
      // trkProp->SetFitReferences(true);
      tracks->SetMainColor(kCyan);
      tracks->SetLineWidth(2.);
      tracks->SetLineStyle(5);
      gEve->AddElement(tracks);
    }
    else
      tracks->DestroyElements();

    for (unsigned int ip = 0; ip < eventReader->Tracks_subdetectorHitNumbers_begin->GetSize(); ip++)
    {
      if (debug) std::cout << "Track " << ip << std::endl;
      unsigned int trackStates_begin = (*eventReader->Tracks_trackStates_begin)[ip];
      unsigned int trackStates_end = (*eventReader->Tracks_trackStates_end)[ip];
      if (debug) std::cout << "  trackStates begin , end = " << trackStates_begin << " , " << trackStates_end << std::endl;
      // store in vectors the track state positions and momenta
      unsigned int nTrackStates = 0;
      if (trackStates_end > trackStates_begin)
        nTrackStates = trackStates_end - trackStates_begin;
      float x[5], y[5], z[5], omega[5], tanLambda[5], phi[5];
      for (unsigned int i=0; i<5; i++) x[i]=-1e9; // set to some large value to check later if track state has been filled
      for (unsigned int its = trackStates_begin; its < trackStates_end; its++) {
        int location = (*eventReader->_Tracks_trackStates_location)[its];
        if (debug) {
          std::cout << "Trackstate " << its << std::endl;
          std::cout << "  location = " << location << std::endl;
        }
        if (location>4) continue;
        x[location] = (*eventReader->_Tracks_trackStates_referencePoint_x)[its];
        y[location] = (*eventReader->_Tracks_trackStates_referencePoint_y)[its];
        z[location] = (*eventReader->_Tracks_trackStates_referencePoint_z)[its];
        omega[location] = (*eventReader->_Tracks_trackStates_omega)[its];
        tanLambda[location] = (*eventReader->_Tracks_trackStates_tanLambda)[its];
        phi[location] = (*eventReader->_Tracks_trackStates_phi)[its];
        if (debug) {
          std::cout << "  x = " << x[location] << std::endl;
          std::cout << "  y = " << y[location] << std::endl;
          std::cout << "  z = " << z[location] << std::endl;
          std::cout << "  r = " << sqrt(x[location]*x[location] + y[location]*y[location]) << std::endl;
          std::cout << "  omega = " << omega[location] << std::endl;
          std::cout << "  phi = " << phi[location] << std::endl;
          std::cout << "  tanLambda = " << tanLambda[location] << std::endl;
          const float FCT = 2.99792458E-4f;
          const float bField = 2.0;
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
      unsigned int hits_begin = (*eventReader->Tracks_subdetectorHitNumbers_begin)[ip];
      unsigned int hits_end = (*eventReader->Tracks_subdetectorHitNumbers_end)[ip];
      if (debug) std::cout << "  subdetectorHitNumbers_begin , end = " << hits_begin << " , " << hits_end << std::endl;
      for (unsigned int ih = hits_begin; ih < hits_end; ih++) {
        nhits[ih-hits_begin] = (*eventReader->_Tracks_subdetectorHitNumbers)[ih];
      }

      // const int tsOrig=2; // which track state to use for track origin: 0=at other, 1=at IP, 2=at first hit, 3=at last hit, 4=at ECAL
      const int tsOrig=1; // which track state to use for track origin: 0=at other, 1=at IP, 2=at first hit, 3=at last hit, 4=at ECAL
      //if (x[tsOrig]<5e5) {  // GM dont remember what is this about..
      if (nTrackStates > 0) {
        TEveRecTrack t;
        float x1 = x[tsOrig]*mm;
        float y1 = y[tsOrig]*mm;
        float z1 = z[tsOrig]*mm;
        const float FCT = 2.99792458E-4f;
        const float bField = 2.0;
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

void EventDisplay::drawClusters(std::string clusterType)
{
  // the cells of the clusters of the current event
  const EventReader::ClusterCells &cells = (clusterType == "topo") ? eventReader->topoClusterCells : eventReader->swClusterCells;
  // read once, not in the loop over the clusters
  const bool drawDirection = displayConfig.getBoolConfig("drawClusterDirection");

  std::vector<std::unique_ptr<CaloCluster>> *clusterData = nullptr;
  if (clusterType == "topo")
  {
    std::cout << "Drawing topo clusters" << std::endl;
    clusterData = &topoclusterData;
  }
  else if (clusterType == "sw")
  {
    std::cout << "Drawing SW clusters" << std::endl;
    clusterData = &swclusterData;
  }
  else
  {
    std::cout << "Unknown cluster type " << clusterType << std::endl;
    return;
  }

  TEvePointSet *clustersCenter = nullptr;
  TEveStraightLineSet *clustersDirection = nullptr;
  TEveElementList *clusters_3D = nullptr;
  TEveElementList *clusters_rhoz = nullptr;
  TEveElementList *clusters_rhophi = nullptr;

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

  // create the graphic containers
  // the palette is ref-counted by the digit sets that use it: hold a reference
  // while filling, released at the end (deletes it if no cluster uses it)
  TEveRGBAPalette *pal = new TEveRGBAPalette(0, 1000);
  pal->IncRefCount();

  // first, create the containers
  // clusters in 3D
  if (debug)
  {
    std::cout << "  Creating the containers to hold the 3D clusters" << std::endl;
  }
  if (clusterType == "topo")
  {
    if (topoclusters_3D == nullptr)
    {
      topoclusters_3D = new TEveElementList(Form("%sclusters (E>%.2f GeV)",
                                                 clusterType.c_str(),
                                                 ClusterEnergyThreshold));
      gEve->AddElement(topoclusters_3D);
    }
    else
    {
      topoclusters_3D->DestroyElements();
      topoclustersCenter = nullptr;
      topoclustersDirection = nullptr;
    }
    clusters_3D = topoclusters_3D;
  }
  else
  {
    if (swclusters_3D == nullptr)
    {
      swclusters_3D = new TEveElementList(Form("%sclusters (E>%.1f GeV)",
                                               clusterType.c_str(),
                                               ClusterEnergyThreshold));
      gEve->AddElement(swclusters_3D);
    }
    else
    {
      swclusters_3D->DestroyElements();
      swclustersCenter = nullptr;
      swclustersDirection = nullptr;
    }
    clusters_3D = swclusters_3D;
  }

  // clusters in 2D
  if (debug)
  {
    std::cout << "  Creating the containers to hold the 2D projections of the clusters" << std::endl;
  }
  if (clusterType == "topo")
  {
    if (topoclusters_rhoz == nullptr)
    {
      topoclusters_rhoz = new TEveElementList(Form("%s clusters in rho-z (E>%.2f GeV)",
                                                   clusterType.c_str(),
                                                   ClusterEnergyThreshold));
      // add to scene (that is not auto-projected!)
      rhoZEventSceneManual->AddElement(topoclusters_rhoz);
      gEve->AddToListTree(topoclusters_rhoz, false);
    }
    else
      topoclusters_rhoz->DestroyElements();

    if (topoclusters_rhophi == nullptr)
    {
      topoclusters_rhophi = new TEveElementList(Form("%s clusters in rho-phi (E>%.2f GeV)",
                                                     clusterType.c_str(),
                                                     ClusterEnergyThreshold));
      // add to scene (that is not auto-projected!)
      rhoPhiEventSceneManual->AddElement(topoclusters_rhophi);
      gEve->AddToListTree(topoclusters_rhophi, false);
    }
    else
      topoclusters_rhophi->DestroyElements();

    clusters_rhoz = topoclusters_rhoz;
    clusters_rhophi = topoclusters_rhophi;
  }
  else
  {
    if (swclusters_rhoz == nullptr)
    {
      swclusters_rhoz = new TEveElementList(Form("%s clusters in rho-z (E>%.2f GeV)",
                                                 clusterType.c_str(),
                                                 ClusterEnergyThreshold));
      // add to scene (that is not auto-projected!)
      rhoZEventSceneManual->AddElement(swclusters_rhoz);
      gEve->AddToListTree(swclusters_rhoz, false);
    }
    else
      swclusters_rhoz->DestroyElements();

    if (swclusters_rhophi == nullptr)
    {
      swclusters_rhophi = new TEveElementList(Form("%s clusters in rho-phi (E>%.2f GeV)",
                                                   clusterType.c_str(),
                                                   ClusterEnergyThreshold));
      // add to scene (that is not auto-projected!)
      rhoPhiEventSceneManual->AddElement(swclusters_rhophi);
      gEve->AddToListTree(swclusters_rhophi, false);
    }
    else
      swclusters_rhophi->DestroyElements();

    clusters_rhoz = swclusters_rhoz;
    clusters_rhophi = swclusters_rhophi;
  }

  if (debug)
  {
    std::cout << "  Filling the containers" << std::endl;
  }
  unsigned int nClusters = (clusterType == "topo") ? eventReader->CaloTopoClusters_position_x->GetSize() : eventReader->CaloClusters_position_x->GetSize();
  if (debug)
  {
    std::cout << "  n(clusters) = " << nClusters << std::endl;
    std::cout << "  Looping over clusters to retrieve barycenter" << std::endl;
  }

  if (clusterType == "topo")
  {
    if (topoclustersCenter == nullptr)
    {
      topoclustersCenter = new TEvePointSet();
      topoclustersCenter->SetName("cluster barycenters");
      topoclustersCenter->SetMarkerColor(kWhite);
      topoclusters_3D->AddElement(topoclustersCenter);
    }
    else
      topoclustersCenter->Reset();
    clustersCenter = topoclustersCenter;

    if (drawDirection) {
      if (topoclustersDirection == nullptr) {
        topoclustersDirection = new TEveStraightLineSet("cluster directions");
        topoclusters_3D->AddElement(topoclustersDirection);
      }
      else
        topoclustersDirection->DestroyElements();
    }
    clustersDirection = topoclustersDirection;
  }
  else
  {
    if (swclustersCenter == nullptr)
    {
      swclustersCenter = new TEvePointSet();
      swclustersCenter->SetName("cluster barycenters");
      swclustersCenter->SetMarkerColor(kGray);
      swclusters_3D->AddElement(swclustersCenter);
    }
    else
      swclustersCenter->Reset();
    clustersCenter = swclustersCenter;

    if (drawDirection) {
      if (swclustersDirection == nullptr) {
        swclustersDirection = new TEveStraightLineSet("cluster directions");
        swclusters_3D->AddElement(swclustersDirection);
      }
      else
        swclustersDirection->DestroyElements();
    }
    clustersDirection = swclustersDirection;

  }
  if (drawDirection) {
    clustersCenter->SetMarkerStyle(4);
    clustersCenter->SetMarkerSize(5);
    clustersDirection->SetLineColor(kWhite);
    clustersDirection->SetLineWidth(5);
  }

  for (unsigned int i = 0; i < nClusters; i++) {
    float E = (clusterType == "topo") ? (*eventReader->CaloTopoClusters_energy)[i] : (*eventReader->CaloClusters_energy)[i];
    if (E < ClusterEnergyThreshold)
      continue;

    float x0(0.), y0(0.), z0(0.);
    float theta(0.), phi(0.);
    // topo cluster positions are in cm while calo clusters and hits/cells in mm ... ?
    if (clusterType == "topo") {
      x0 = (*eventReader->CaloTopoClusters_position_x)[i] * mm;
      y0 = (*eventReader->CaloTopoClusters_position_y)[i] * mm;
      z0 = (*eventReader->CaloTopoClusters_position_z)[i] * mm;
      if (drawDirection) {
        theta = (*eventReader->CaloTopoClusters_theta)[i];
        phi = (*eventReader->CaloTopoClusters_phi)[i];
      }
    }
    else {
      x0 = (*eventReader->CaloClusters_position_x)[i] * mm;
      y0 = (*eventReader->CaloClusters_position_y)[i] * mm;
      z0 = (*eventReader->CaloClusters_position_z)[i] * mm;
      if (drawDirection) {
        theta = (*eventReader->CaloClusters_theta)[i];
        phi = (*eventReader->CaloClusters_phi)[i];
      }
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
    float energy, xcl, ycl, zcl;
    if (clusterType == "topo")
    {
      energy = (*eventReader->CaloTopoClusters_energy)[i];
      xcl = (*eventReader->CaloTopoClusters_position_x)[i];
      ycl = (*eventReader->CaloTopoClusters_position_y)[i];
      zcl = (*eventReader->CaloTopoClusters_position_z)[i];
    }
    else
    {
      energy = (*eventReader->CaloClusters_energy)[i];
      xcl = (*eventReader->CaloClusters_position_x)[i];
      ycl = (*eventReader->CaloClusters_position_y)[i];
      zcl = (*eventReader->CaloClusters_position_z)[i];
    }
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
      aqs->SetPalette(pal);
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
      aqs2->SetPalette(pal);
      aqs2->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",
                          energy,
                          rcl,
                          thetacl,
                          phicl));
      qs_rhophi.push_back(aqs2);
      clusters_rhophi->AddElement(aqs2);
      //}

      TEveBoxSet *_bs = new TEveBoxSet(Form("%s cluster %d", clusterType.c_str(), (int)i), "");
      _bs->Reset(TEveBoxSet::kBT_FreeBox, false, 32);
      _bs->SetMainTransparency(80);
      _bs->SetOwnIds(kTRUE);
      _bs->SetPalette(pal);
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
        for (unsigned int iLayer = 0; iLayer < clus->getNLayersECalBarrel(); iLayer++)
        {
          if (clus->getEnergyInECalBarrelLayer(iLayer) > 0)
          {
            layerBarycenters->SetNextPoint(
                                           clus->getBarycenterInECalBarrelLayer(iLayer).X() / cm,
                                           clus->getBarycenterInECalBarrelLayer(iLayer).Y() / cm,
                                           clus->getBarycenterInECalBarrelLayer(iLayer).Z() / cm);
          }
        }
        for (unsigned int iLayer = 0; iLayer < clus->getNLayersECalEndCap(); iLayer++)
        {
          if (clus->getEnergyInECalEndCapLayer(iLayer) > 0)
          {
            layerBarycenters->SetNextPoint(
                                           clus->getBarycenterInECalEndCapLayer(iLayer).X() / cm,
                                           clus->getBarycenterInECalEndCapLayer(iLayer).Y() / cm,
                                           clus->getBarycenterInECalEndCapLayer(iLayer).Z() / cm);
          }
        }
        for (unsigned int iLayer = 0; iLayer < (*clusterData)[iclAboveThreshold]->getNLayersHCalBarrel(); iLayer++)
        {
          if (clus->getEnergyInHCalBarrelLayer(iLayer) > 0)
          {
            layerBarycenters->SetNextPoint(
                                           clus->getBarycenterInHCalBarrelLayer(iLayer).X() / cm,
                                           clus->getBarycenterInHCalBarrelLayer(iLayer).Y() / cm,
                                           clus->getBarycenterInHCalBarrelLayer(iLayer).Z() / cm);
          }
        }
        for (unsigned int iLayer = 0; iLayer < (*clusterData)[iclAboveThreshold]->getNLayersHCalEndCap(); iLayer++)
        {
          if (clus->getEnergyInHCalEndCapLayer(iLayer) > 0)
          {
            layerBarycenters->SetNextPoint(
                                           clus->getBarycenterInHCalEndCapLayer(iLayer).X() / cm,
                                           clus->getBarycenterInHCalEndCapLayer(iLayer).Y() / cm,
                                           clus->getBarycenterInHCalEndCapLayer(iLayer).Z() / cm);
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
    unsigned int first_hit, last_hit;
    if (clusterType == "topo")
    {
      first_hit = (*eventReader->CaloTopoClusters_hits_begin)[j];
      last_hit = (*eventReader->CaloTopoClusters_hits_end)[j];
    }
    else
    {
      first_hit = (*eventReader->CaloClusters_hits_begin)[j];
      last_hit = (*eventReader->CaloClusters_hits_end)[j];
    }
    // in case of overlapping ranges, keep the first cluster
    for (unsigned int i = first_hit; i < std::min(last_hit, nCells); i++)
      if (cellCluster[i] == -1)
        cellCluster[i] = j;
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
    float r_in(0.), r_out(0.);
    int layer(0), thetaID(0), moduleID(0), phiID(0);
    ProjectedCellKey rhophiID{}, rhozID{};

    if (systemID == DetectorGeometry::kECalBarrel)
    {
      layer = (int)DetectorGeometry::ECalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::ECalBarrelThetaBin(cellID);
      moduleID = (int)DetectorGeometry::ECalBarrelModule(cellID);
      rhophiID = {systemID, icl, layer, moduleID};
      rhozID = {systemID, icl, layer, thetaID};
      r_in = geomReader->r[layer];
      r_out = geomReader->r[layer + 1];
    }
    else if (systemID == DetectorGeometry::kHCalBarrel)
    {
      layer = (int)DetectorGeometry::HCalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::HCalBarrelThetaBin(cellID);
      phiID = (int)DetectorGeometry::HCalBarrelPhiBin(cellID);
      rhophiID = {systemID, icl, layer, phiID};
      rhozID = {systemID, icl, layer, thetaID};
      r_in = geomReader->rHCal[layer];
      r_out = geomReader->rHCal[layer + 1];
    }
    else if (systemID == DetectorGeometry::kECalEndcap || systemID == DetectorGeometry::kHCalEndcap || systemID == DetectorGeometry::kMuonBarrel || systemID == DetectorGeometry::kMuonEndcap)
    {
      // cells of ECAL endcap, HCAL endcap, MUON barrel and endcap are not drawn
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
      if (cellEnergies_rhophi.count(rhophiID))
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
    if (bs[icl] != nullptr)
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
  /*
  for (auto& [key, energy]: cellEnergies_rhoz) {
    int ikey = (int) key;
    int icl = ikey % (512*16);
    int layer = (ikey / 512) % 16;
    int thetaID = ikey / (512*16) ;

    float r_in = geomReader->r[layer];
    float r_out = geomReader->r[layer+1];
    float x_center, y_center, z_center;
    if (clusterType=="topo") {
      x_center = cells.x[i] * mm;
      y_center = cells.y[i] * mm;
      z_center = cells.z[i] * mm;
    }
    else {
      x_center = cells.x[i] * mm;
      y_center = cells.y[i] * mm;
      z_center = cells.z[i] * mm;
    }
    float r_center = sqrt(x_center*x_center + y_center*y_center);

    float theta_center = atan2(r_center, z_center);
    //float phi_center = atan2(y_center, x_center);

    float verts[12];
    verts[0] = r_in / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
    verts[1] = r_in * MathUtils::sign(y_center);
    verts[3] = r_in / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
    verts[4] = r_in * MathUtils::sign(y_center);
    verts[6] = r_out / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
    verts[7] = r_out * MathUtils::sign(y_center);
    verts[9] = r_out / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
    verts[10] = r_out * MathUtils::sign(y_center);
    verts[2] = verts[5] = verts[8] = verts[11] = 0.;
    //if (qs_rhoz[icl]!=nullptr) {
    qs_rhoz[icl]->AddQuad(verts);
    qs_rhoz[icl]->QuadValue( (int) (1000 * energy) );
    qs_rhoz[icl]->QuadId( new TNamed(Form("%s %d cell L%d Th%d", clusterType.c_str(), icl, layer, thetaID), "Dong!") );
    //}
  }
  */

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
      rhozID = {system, icl, layer, thetaID};
      r_in = geomReader->r[layer];
      r_out = geomReader->r[layer + 1];
    }
    else if (system == DetectorGeometry::kHCalBarrel)
    {
      layer = (int)DetectorGeometry::HCalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::HCalBarrelThetaBin(cellID);
      rhozID = {system, icl, layer, thetaID};
      r_in = geomReader->rHCal[layer];
      r_out = geomReader->rHCal[layer + 1];
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
    verts[2] = verts[5] = verts[8] = verts[11] = 0.;
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
        verts[2] = verts[5] = verts[8] = verts[11] = 0.;
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
      verts[2] = verts[5] = verts[8] = verts[11] = 0.;
      qs_rhophi[icl]->AddQuad(verts);
      qs_rhophi[icl]->QuadValue((int)(1000 * energy));
      qs_rhophi[icl]->QuadId(new TNamed(Form("%s %d cell L%d P%d", clusterType.c_str(), icl, layer, phiID), "Dong!"));
    }
  }

  pal->DecRefCount();
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
  unsigned int nGenParticles = eventReader->genParticles_generatorStatus ? eventReader->genParticles_generatorStatus->GetSize() : 0;
  for (unsigned int i = 0; i < nGenParticles; i++)
  {
    if ((*eventReader->genParticles_generatorStatus)[i] != 1)
      continue;
    float _px = (*eventReader->genParticles_momentum_x)[i];
    float _py = (*eventReader->genParticles_momentum_y)[i];
    float _pz = (*eventReader->genParticles_momentum_z)[i];
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
  int pdgID = (ipmax >= 0) ? (*eventReader->genParticles_PDG)[ipmax] : 0;

  // draw the various event elements
  drawParticles();
  drawHits();
  drawDigis();
  drawTracks();
  if (displayConfig.getBoolConfig("drawTopoClusters"))
  {
    fillClusters("topo");
    drawClusters("topo");
  }
  if (displayConfig.getBoolConfig("drawCaloClusters"))
  {
    fillClusters("sw");
    drawClusters("sw");
  }

  // add process and event label
  auto processStr = (ipmax >= 0) ? Form("%s, %.1f GeV", ParticleUtils::name(pdgID), pmax) : "no primary particle";
  if (processLabel == nullptr) {
    float x(0.05), y(0.95);
    processLabel = new TGLConstAnnotation(mainGLView, processStr, x, y);
    processLabel->SetTextSize(0.05); // % of window diagonal
    processLabel->SetAllowClose(false);
  }
  else {
    processLabel->SetText(processStr);
  }
  if (eventLabel == nullptr) {
    auto eventStr = Form("Event %d", eventId);
    float x(0.05), y(0.90);
    eventLabel = new TGLConstAnnotation(mainGLView, eventStr, x, y);
    eventLabel->SetTextSize(0.05); // % of window diagonal
    eventLabel->SetAllowClose(false);
  }
  else {
    eventLabel->SetText(Form("Event %d", eventId));
  }

  TEveElement *top = (TEveElement *)gEve->GetCurrentEvent();
  // if nothing is drawn, top is null
  if (top)
  {
    rhoPhiEventScene->DestroyElements();
    rhoPhiProjManager->ImportElements(top, rhoPhiEventScene);

    // slow??
    // TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(pal, 0.55, 0.1, 0.4, 0.05);
    // rhoPhiGLView->AddOverlayElement(po);

    rhoZEventScene->DestroyElements();
    rhoZProjManager->ImportElements(top, rhoZEventScene);

    // TEveRGBAPaletteOverlay *po2 = new TEveRGBAPaletteOverlay(pal, 0.55, 0.1, 0.4, 0.05);
    // rhoZGLView->AddOverlayElement(po2);
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
  // read some display setup from config file
  ClusterEnergyThreshold = displayConfig.getFloatConfig("energyThresholdClusters");
  ParticleEnergyThreshold = displayConfig.getFloatConfig("energyThresholdParticles");
  HitEnergyThreshold = displayConfig.getFloatConfig("energyThresholdHits");
  CellEnergyThreshold = displayConfig.getFloatConfig("energyThresholdCells");

  // calculate the geometry parameters
  // detectorVersion: 2 = ECAL barrel with 12 layers (ALLEGRO o1_v02)
  //                  3 = ECAL barrel with 11 layers (ALLEGRO o1_v03, o2_v01)
  if (geomFile.find("o1_v03") != std::string::npos || geomFile.find("o2_v01") != std::string::npos)
    detectorVersion = 3;
  else
    detectorVersion = 2;
  geomReader = std::make_unique<DetectorGeometry>(detectorVersion);

  // create magnetic field
  std::string bFieldType = displayConfig.getStringConfig("bField");
  if (bFieldType != "GlobalSolenoid") {
    std::string fieldMapFile = PathUtils::findFile(displayConfig.getStringConfig("bFieldMap"));
    std::cout << "Using magnetic field map " << fieldMapFile << std::endl;
    try {
      CylField* fieldMap = new CylField(fieldMapFile, "fieldmap");
      // Note: field is reversed in EVE (see https://root.cern/doc/v636/classTEveTrackPropagator.html)
      fieldMap->setReverseState(true);
      magField = fieldMap;
    }
    catch (const std::exception& e) {
      std::cerr << "ERROR: cannot read the magnetic field map " << fieldMapFile << " (" << e.what()
                << "), using the ideal solenoid magnetic field instead" << std::endl;
    }
  }
  if (!magField) {
    std::cout << "Using ideal solenoid magnetic field" << std::endl;
    MagField* solenoid = new MagField(geomReader->Bin, geomReader->Bout, geomReader->zMax, geomReader->rMax);
    // Note: field is reversed in EVE (see https://root.cern/doc/v636/classTEveTrackPropagator.html)
    solenoid->setReverseState(true);
    magField = solenoid;
  }
  std::cout << "\n" << std::endl;

  std::cout << "******************************************************************************" << std::endl;
  std::cout << "Displaying the geometry" << std::endl;
  std::cout << "******************************************************************************" << std::endl
            << std::endl;

  // create the eve manageer
  TEveManager::Create();
  browser = gEve->GetBrowser();

  // Set title of main window
  browser->SetWindowName("FCC-ee ALLEGRO detector event display");

  // see palettes here: https://root.cern.ch/doc/master/classTColor.html
  // gStyle->SetPalette(kAvocado);
  gStyle->SetPalette(kSienna);

  // first tab
  mainGLView = gEve->GetDefaultGLViewer();
  int bkgColor = displayConfig.getIntConfig("bkgColor3D");
  mainGLView->ColorSet().Background().SetColor(bkgColor>=0 ? bkgColor : kBlack);
  mainGLView->SetGuideState(TGLUtil::kAxesOrigin, false, false, 0);
  mainGLView->DrawGuides();
  gEve->GetDefaultViewer()->SetElementName("3D view");
  mainGLView->CurrentCamera().RotateRad(-.8, 0.5);
  browser->GetTabRight()->Connect("Selected(Int_t)", "EventDisplay", this, "onTabSelected(Int_t)");

  // Create the geometry and the readouts
  // TEveElementList *
    geom = new TEveElementList("Geometry");
  TEveElementList *PCBs = new TEveElementList("PCBs");
  TEveElementList *actives = new TEveElementList("Active elements");
  TEveElementList *passives = new TEveElementList("Passive elements");
  readout = new TEveElementList("Readout");
  readout->SetRnrSelfChildren(showReadout, showReadout);
  if (useG4geom)
  {
    // auto fGeom = TFile::Open(geomFile.c_str(), "CACHEREAD");
    std::cout << "Reading Geant4 geometry from file " << geomFile << std::endl;
    auto fGeom = TFile::Open(geomFile.c_str(), "READ");
    if (!fGeom)
      return;

    gEve->AddGlobalElement(geom);
    gEve->AddToListTree(geom, true);

    TEveGeoShapeExtract *gse = (TEveGeoShapeExtract *)fGeom->Get("world");
    if (!gse)
      return;

    TEveGeoShape *world = TEveGeoShape::ImportShapeExtract(gse, 0);
    world->SetMainTransparency(100);
    fGeom->Close();
    delete fGeom;

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
    else
    {

      //
      // ECAL barrel
      //
      TPRegexp re("ECalBarrel*");
      TEveElement *ecalBarrelEl = world->FindChild(re);
      if (!ecalBarrelEl)
      {
        std::cerr << "ERROR: ECalBarrel volume not found in geometry file " << geomFile << std::endl;
        return;
      }
      ecalBarrelEl->SetPickableRecursively(kTRUE);
      ((TEveGeoShape *)ecalBarrelEl)->SetNSegments(256);
      geom->AddElement(ecalBarrelEl);
      // hide the envelope, only show physical volumes in it
      ecalBarrelEl->SetRnrSelfChildren(false, true);

      // set transparency of the subvolumes of the bath
      re = TPRegexp("LAr_bath*");
      TEveElement *bath = ecalBarrelEl->FindChild(re);
      TEveElement::List_t matches;
      re = TPRegexp("ECAL_Cryo*");
      ecalBarrelEl->FindChildren(matches, re);
      for (TEveElement *a : matches)
      {
        a->SetMainTransparency(70);
        ((TEveGeoShape *)a)->SetNSegments(256);
      }
      re = TPRegexp("services*");
      ecalBarrelEl->FindChildren(matches, re);
      for (auto a : matches)
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
      }

      // HCAL
      TEveElement *hcalBarrelEl = nullptr;
      if (doHCal)
      {
        re = TPRegexp("HCalEnvelopeVolume*");
        hcalBarrelEl = world->FindChild(re);
        if (!hcalBarrelEl)
          std::cout << "WARNING: HCalEnvelopeVolume not found in geometry file, not drawing HCAL" << std::endl;
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
        hcalBarrelEl->FindChildren(matches, re);
        for (auto a : matches)
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
        //
        // HCAL endcaps
        //
        if (doEndcaps)
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
      }
    }
  }
  else
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
  gEve->AddToListTree(readout, true);
  gEve->FullRedraw3D(true);

  // create second tab (R-phi view)
  rhoPhiView = gEve->SpawnNewViewer("Projection Rho-Phi");

  // two scenes, for geometry and event
  rhoPhiScene = gEve->SpawnNewScene("Rho-Phi geometry",
                                    "Scene holding projected geometry data for the Rho-Phi view.");
  rhoPhiView->AddScene(rhoPhiScene);
  if (evtFile != "")
  {
    rhoPhiEventScene = gEve->SpawnNewScene("Rho-Phi Event Data",
                                           "Scene holding projected event-data for the Rho-Phi view.");
    rhoPhiView->AddScene(rhoPhiEventScene);
  }

  rhoPhiEventSceneManual = gEve->SpawnNewScene("Rho-Phi Event Data 2",
                                               "Scene holding hand-crafted event-data for the Rho-Phi view.");
  rhoPhiView->AddScene(rhoPhiEventSceneManual);
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
  if (useG4geom)
    rhoPhiProjManager->ImportElements(geom, rhoPhiScene);
  else
  {
    rhoPhiProjManager->ImportElements(ecalbarrel, rhoPhiScene);
    if (doHCal)
      rhoPhiProjManager->ImportElements(hcalbarrel, rhoPhiScene);
  }

  {
    // do not draw the cryo side in the rho-phi view as it renders poorly
    // and turn on the active elements so we can see the cell edges
    TEveElement *element = rhoPhiScene->FindChild("Geometry [P]");
    if (element)
    {
      TPRegexp re("ECalBarrel_vol_(\\w+)");
      TEveElement *projbarrel = element->FindChild(re);
      if (projbarrel)
      {
        TPRegexp reCryo("ECAL_Cryo_side_(\\w+)");
        TEveElement *projcryoside = projbarrel->FindChild(reCryo);
        if (projcryoside)
        {
          projcryoside->SetRnrSelf(false);
        }
        TEveElementList *projbath = (TEveElementList *)projbarrel->FindChild("LAr_bath [P]");
        if (projbath)
        {
          projbath->SetRnrSelfChildren(true, true);
          TEveElementList *projpassives = (TEveElementList *)projbath->FindChild("Passive elements [P]");
          if (projpassives)
          {
            projpassives->SetRnrSelfChildren(false, false);
          }
          TEveElementList *projpcbs = (TEveElementList *)projbath->FindChild("PCBs [P]");
          if (projpcbs)
          {
            projpcbs->SetRnrSelfChildren(false, false);
          }
        }
      }
    }

    // when drawing the full detector - do not draw endcap volumes in rho-phi view
    // and show active volumes in ecal barrel
    if (showFullDetector && element)
    {
      for (TEveElement::List_i itr = element->BeginChildren(); itr != element->EndChildren(); itr++)
      {
        TEveElement *a = *itr;
        TString s(a->GetElementName());
        // do not draw endcap volumes
        if (s.Contains("endcap"))
          a->SetRnrSelfChildren(false, false);

        // turn off the envelopes of the vertex barrel
        else if (s.BeginsWith("Vertex barrel")) {
          a->SetRnrSelfChildren(false, true);
          for (TEveElement::List_i itrChild = a->BeginChildren(); itrChild != a->EndChildren(); itrChild++) {
            TEveElement *b = *itrChild;
            TString childName(b->GetElementName());
            if (childName.Contains("assembly")) {
              b->SetRnrSelfChildren(false, true);
            }
          }
        }

        // show active volumes in ecal barrel
        else if (s.BeginsWith("ECal barrel"))
        {
          TPRegexp re("ECalBarrel_vol_(\\w+)");
          TEveElement *projbarrel = a->FindChild(re);
          if (!projbarrel)
            continue;
          projbarrel->SetRnrSelfChildren(false, true);
          TPRegexp rebath("LAr_bath(\\w+)");
          TEveElement *projbath = projbarrel->FindChild(rebath);
          if (projbath)
          {
            projbath->SetRnrSelfChildren(false, true);
            for (TEveElement::List_i itr2 = projbath->BeginChildren(); itr2 != projbath->EndChildren(); itr2++)
            {
              TEveElement *el = *itr2;
              TString s2(el->GetElementName());
              if (s2.BeginsWith("passive_") || s2.BeginsWith("PCB_"))
                el->SetRnrSelfChildren(false, false);
              else if (s2.BeginsWith("active_"))
              {
                for (TEveElement::List_i itr3 = el->BeginChildren(); itr3 != el->EndChildren(); itr3++)
                  (*itr3)->SetMainTransparency(100);
              }
            }
          }
        }

        // show layers in HCAL barrel
        else if (s.BeginsWith("HCal barrel"))
        {
          if (a->NumChildren() == 0)
            continue;
          TEveElement *envelope = *(a->BeginChildren());
          envelope->SetRnrSelfChildren(false, true);
          ((TEveGeoShape *)envelope)->SetDrawFrame(false);
          for (TEveElement::List_i itr2 = envelope->BeginChildren(); itr2 != envelope->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            // std::cout << elname << std::endl;
            if (elname.Contains("EndPlate"))
              el->SetRnrSelfChildren(false, false);
            else {
              el->SetRnrSelfChildren(true, false);
              if (elname.Contains("LayerVol"))
                el->SetMainColor(kAzure-7);
              el->SetMainTransparency(0); // check if this does not hide the hits
            }
          }
        }

        // show layers in MUON barrel
        else if (s.BeginsWith("Muon tagger barrel"))
        {
          if (a->NumChildren() == 0)
            continue;
          TEveElement *envelope = *(a->BeginChildren());
          envelope->SetRnrSelfChildren(false, true);
          ((TEveGeoShape *)envelope)->SetDrawFrame(false);
          for (TEveElement::List_i itr2 = envelope->BeginChildren(); itr2 != envelope->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            el->SetRnrSelfChildren(true, false);
            if (overrideDetectorColors) el->SetMainColor(kOrange);
            el->SetMainTransparency(0);
            ((TEveGeoShape *)el)->SetDrawFrame(true);
          }
        }
      }
    }
  }

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

  // third tab (R-z view)
  rhoZView = gEve->SpawnNewViewer("Projection Rho-Z");
  rhoZScene = gEve->SpawnNewScene("Rho-Z geometry",
                                  "Scene holding projected geometry data for the Rho-Z view.");
  rhoZView->AddScene(rhoZScene);
  if (evtFile != "")
  {
    rhoZEventScene = gEve->SpawnNewScene("Rho-Z Event Data",
                                         "Scene holding projected event-data for the Rho-Z view.");
    rhoZView->AddScene(rhoZEventScene);
  }
  rhoZEventSceneManual = gEve->SpawnNewScene("Rho-Z Event Data 2",
                                             "Scene holding hand-crafted event-data for the Rho-Z view.");
  rhoZView->AddScene(rhoZEventSceneManual);
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
  if (useG4geom)
    rhoZProjManager->ImportElements(geom, rhoZScene);
  else
  {
    rhoZProjManager->ImportElements(ecalbarrel, rhoZScene);
    if (doHCal)
      rhoZProjManager->ImportElements(hcalbarrel, rhoZScene);
  }

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

  // turn on the active elements in rho-Z so we can see the cell edges
  {
    TEveElement *element = rhoZScene->FindChild("Geometry [P]");
    if (element)
    {
      TPRegexp re("ECalBarrel_vol_(\\w+)");
      TEveElement *projbarrel = element->FindChild(re);
      if (projbarrel)
      {
        TEveElementList *projbath = (TEveElementList *)projbarrel->FindChild("LAr_bath [P]");
        if (projbath)
        {
          projbath->SetRnrSelfChildren(true, true);
          TEveElementList *projpassives = (TEveElementList *)projbath->FindChild("Passive elements [P]");
          if (projpassives)
          {
            projpassives->SetRnrSelfChildren(false, false);
          }
          TEveElementList *projpcbs = (TEveElementList *)projbath->FindChild("PCBs [P]");
          if (projpcbs)
          {
            projpcbs->SetRnrSelfChildren(false, false);
          }
        }
      }
    }

    if (showFullDetector && element)
    {
      for (TEveElement::List_i itr = element->BeginChildren(); itr != element->EndChildren(); itr++) {
        TEveElement *a = *itr;
        TString s(a->GetElementName());

        // MDI, beampipe, lumical can be enabled
        if (s.BeginsWith("MDI") || s.BeginsWith("Beampipe") || s.BeginsWith("LumiCal")) {
          a->SetRnrSelfChildren(true, true);
        }
        // turn off the envelopes of the vertex barrel
        else if (s.BeginsWith("Vertex barrel")) {
          a->SetRnrSelfChildren(false, true);
          for (TEveElement::List_i itrChild = a->BeginChildren(); itrChild != a->EndChildren(); itrChild++) {
            TEveElement *b = *itrChild;
            TString childName(b->GetElementName());
            if (childName.Contains("assembly")) {
              b->SetRnrSelfChildren(false, true);
            }
          }
        }
        // show vessel and gas for DCH
        else if (s.BeginsWith("Drift chamber")) {
          a->SetRnrSelfChildren(false, true);
          for (TEveElement::List_i itrChild = a->BeginChildren(); itrChild != a->EndChildren(); itrChild++) {
            TEveElement *b = *itrChild;
            b->SetRnrSelfChildren(true, false);
            if (overrideDetectorColors) b->SetMainColor(kGreen-5);
            ((TEveGeoShape *)b)->SetDrawFrame(false);
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
          TPRegexp re("ECalBarrel_vol_(\\w+)");
          TEveElement *projbarrel = a->FindChild(re);
          if (!projbarrel)
            continue;
          projbarrel->SetMainTransparency(useTransparencies ? 60 : 0);
          TPRegexp rebath("LAr_bath(\\w+)");
          TEveElement *projbath = projbarrel->FindChild(rebath);
          if (projbath) {
            projbath->SetRnrSelfChildren(false, true);
            for (TEveElement::List_i itr2 = projbath->BeginChildren(); itr2 != projbath->EndChildren(); itr2++) {
              TEveElement *el = *itr2;
              TString s2(el->GetElementName());
              if (s2.BeginsWith("passive_") || s2.BeginsWith("PCB_"))
                el->SetRnrSelfChildren(false, false);
            }
          }
        }
        // HCAL barrel - do not draw the envelope, rather show the layers
        else if (s.BeginsWith("HCal barrel"))
        {
          if (a->NumChildren() == 0)
            continue;
          TEveElement *envelope = *(a->BeginChildren());
          envelope->SetRnrSelfChildren(false, true);
          ((TEveGeoShape *)envelope)->SetDrawFrame(false);
          // loop over the children of the envelope
          for (TEveElement::List_i itr2 = envelope->BeginChildren(); itr2 != envelope->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            // std::cout << elname << std::endl;
            el->SetRnrSelfChildren(true, false);
            if (elname.Contains("HCalLayerVol")) {
              if (overrideDetectorColors) el->SetMainColor(kAzure-7);
              el->SetMainTransparency(0);
              // note: can do el->SetRnrSelfChildren(false, true); to see detailed row structure (tile/abs/..)
            }
          }
        }
        // HCAL endcap - do not draw the envelope
        else if (s.BeginsWith("HCal endcap"))
        {
          TPRegexp re("HCalThreePartsEndcap*");
          TEveElement *envelope = a->FindChild(re);
          if (!envelope)
            continue;
          envelope->SetRnrSelfChildren(false, true);
          ((TEveGeoShape *)envelope)->SetDrawFrame(false);
          for (TEveElement::List_i itr2 = envelope->BeginChildren(); itr2 != envelope->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            el->SetRnrSelfChildren(true, false);
            if (elname.Contains("HCalECLayerVol")) {
              if (overrideDetectorColors) el->SetMainColor(kAzure-7);
              el->SetMainTransparency(0);
              // note: can do el->SetRnrSelfChildren(false, true); to see detailed row structure (tile/abs/..)
            }
          }
        }
        // MUON barrel - do not draw the envelopes, rather show the layers
        else if (s.BeginsWith("Muon tagger barrel"))
        {
          if (a->NumChildren() == 0)
            continue;
          TEveElement *envelope = *(a->BeginChildren());
          envelope->SetRnrSelfChildren(false, true);
          ((TEveGeoShape *)envelope)->SetDrawFrame(false);
          for (TEveElement::List_i itr2 = envelope->BeginChildren(); itr2 != envelope->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            el->SetRnrSelfChildren(true, false);
            if (overrideDetectorColors) el->SetMainColor(kOrange);
            el->SetMainTransparency(0);
            ((TEveGeoShape *)el)->SetDrawFrame(true);
          }
        }
        // MUON endcap - do not draw the envelope, rather show the layers
        else if (s.BeginsWith("Muon tagger endcap"))
        {
          for (TEveElement::List_i itrenv = a->BeginChildren(); itrenv != a->EndChildren(); itrenv++) {
            TEveElement *envelope = *itrenv;
            envelope->SetRnrSelfChildren(false, true);
            ((TEveGeoShape *)envelope)->SetDrawFrame(false);
            for (TEveElement::List_i itr2 = envelope->BeginChildren(); itr2 != envelope->EndChildren(); itr2++)
            {
              TEveElement *el = *itr2;
              TString elname(el->GetElementName());
              el->SetRnrSelfChildren(true, false);
              if (overrideDetectorColors) el->SetMainColor(kOrange);
              el->SetMainTransparency(0);
              ((TEveGeoShape *)el)->SetDrawFrame(true);
            }
          }
        }
      }
    }
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

  gEve->DoRedraw3D();
  gEve->FullRedraw3D(true);


  // create the gui for event navigation
  makeGui();

  //
  // display the events
  //
  if (debug)
    std::cout << "evtFile: " << evtFile << std::endl;
  if (evtFile != "")
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
    eventReader = std::make_unique<EventReader>(eventFile.get(), doHCal);

    // setup the branches
    eventReader->setBranches();

    // decode the cellIDs with the CellIDEncoding stored in the file
    setCellIDEncodings();

    // print updated draw settings based on info found in ROOT file
    displayConfig.print();

    // read the number of events in the file
    nEvents = eventReader->nEvents;
    updateEventEntries();

    // load and display the requested event
    std::cout << std::endl;
    std::cout << "******************************************************************************" << std::endl;
    std::cout << "Reading the events" << std::endl;
    std::cout << "******************************************************************************" << std::endl
         << std::endl;
    drawEvent(initialEvent);
  }

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

  // the entry that emitted the signal (one per tab)
  auto it = std::find(eventEntries.begin(), eventEntries.end(), static_cast<TGNumberEntry *>(gTQSender));
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

void EventDisplay::updateEventEntries()
{
  for (TGNumberEntry *entry : eventEntries)
  {
    if (nEvents > 0)
      entry->SetLimits(TGNumberFormat::kNELLimitMinMax, 0, nEvents - 1);
    // do not emit ValueSet, which would reload the event
    entry->SetIntNumber(eventId, kFALSE);
  }
  for (TGLabel *label : eventCountLabels)
  {
    label->SetText(Form("of %d", nEvents));
    label->Resize(label->GetDefaultSize());
  }
  for (TGLabel *label : eventCountLabels)
    if (label->GetParent())
      ((TGCompositeFrame *) label->GetParent())->Layout();
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
    el->SetRnrSelf(showReadout);
  }
  gEve->FullRedraw3D(false);
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

    TGHorizontalFrame* hf = new TGHorizontalFrame(frm);
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
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "takeScreenshot()");
    b->SetToolTipText("Take a screenshot of the selected view");

    b = makeButton("GoBack.gif", "<");
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "bck()");
    b->SetToolTipText("Go to previous event");

    // event number entry: type the number and press Enter, or use its arrows
    TGNumberEntry *entry = new TGNumberEntry(hf, 0, 6, -1, TGNumberFormat::kNESInteger,
                                             TGNumberFormat::kNEANonNegative);
    entry->GetNumberEntry()->SetToolTipText("Event number: type it and press Enter to go to that event");
    entry->Connect("ValueSet(Long_t)", "EventDisplay", this, "goToEvent()");
    entry->SetState(evtFile != "");
    hf->AddFrame(entry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 10, 10));
    eventEntries.push_back(entry);
    TGLabel *countLabel = new TGLabel(hf, "of 0     ");
    hf->AddFrame(countLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 10, 10));
    eventCountLabels.push_back(countLabel);

    b = makeButton("GoForward.gif", ">");
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 10, 10, 10));
    b->SetToolTipText("Go to next event");
    b->Connect("Clicked()", "EventDisplay", this, "fwd()");

    b = makeButton("grid_icon.png", "Readout");
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 10, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "toggleReadout()");
    b->SetToolTipText("Hide/Show readouts");

    frm->AddFrame(hf);
    frm->Layout();
    frm->MapSubwindows();
  }

  // hide command tab
  gEve->GetBrowser()->HideBottomTab();
}

void EventDisplay::onTabSelected(Int_t tab)
{
  if (tab==0) {
    activeGLViewer = mainGLView;
  }
  else if (tab==1) {
    activeGLViewer = rhoPhiGLView;
    if (!initRhoPhiView) {
      rhoPhiView->Redraw(true);
      initRhoPhiView = true;
    }
  }
  else if (tab==2) {
    activeGLViewer = rhoZGLView;
    if (!initRhoZView) {
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

  // Create propagator
  TEveTrackPropagator* propagator = new TEveTrackPropagator();
  propagator->SetMagFieldObj(magField);
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
    propagator->InitTrack(start_pos, charge);
    TEveVectorD target_pos(g_vx, g_vy, g_vz);

    TEveVectorD propagated_mom(currentMomentum.X(), currentMomentum.Y(), currentMomentum.Z());
    Bool_t success = propagator->GoToVertex(target_pos, propagated_mom);
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
  TVector3 endpoint((*eventReader->genParticles_endpoint_x)[electronIndex]*mm,
                    (*eventReader->genParticles_endpoint_y)[electronIndex]*mm,
                    (*eventReader->genParticles_endpoint_z)[electronIndex]*mm);
  segments.push_back({
      TVector3(currentVertex.X(), currentVertex.Y(), currentVertex.Z()),
      TVector3(endpoint.X(), endpoint.Y(), endpoint.Z()),
      currentMomentum
    });
  return segments;
}

void EventDisplay::createElementList(TEveElementList*& list, const std::string& listName)
{
  if (list == nullptr)
  {
    list = new TEveElementList(listName.c_str());
    gEve->AddElement(list);
  }
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
  const std::vector<std::pair<std::string, DetectorGeometry::Readout>> readouts = {
    {"ecalBarrelCells", DetectorGeometry::Readout::ECalBarrel},
    {"ecalEndcapCells", DetectorGeometry::Readout::ECalEndcap},
    {"hcalBarrelCells", DetectorGeometry::Readout::HCalBarrel},
    {"hcalEndcapCells", DetectorGeometry::Readout::HCalEndcap}};
  for (const auto& [key, readoutType] : readouts)
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

void EventDisplay::addSimpleElement(TEveElement* el, SubDetector& det)
{
  if (! det.list) {
    det.list = new TEveElementList(det.name.c_str());
    geom->AddElement(det.list);
    if (! det.visibleByDefault)
      det.list->SetRnrSelfChildren(false, false);
  }
  if (debug) {
    std::cout << "Adding " << el->GetElementName() << " to " << det.name << std::endl;
  }
  det.list->AddElement(el);
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
  el->SetRnrSelfChildren(true, false);
  el->SetMainTransparency(useTransparencies ? 60 : 0);
  if (overrideDetectorColors) el->SetMainColor(kGreen-5);
  ((TEveGeoShape *)el)->SetNSegments(128);
  ((TEveGeoShape *)el)->SetDrawFrame(false);
}

void EventDisplay::addWrapper(TEveElement* el, SubDetector& det)
{
  // TODO
  // cannot use SimpleElement because it will add assembly to element list so I have to duplicate code. Need to refactor
  // addSimpleElement(el, det);
  if (! det.list) {
    det.list = new TEveElementList(det.name.c_str());
    geom->AddElement(det.list);
    if (! det.visibleByDefault)
      det.list->SetRnrSelfChildren(false, false);
  }
  if (debug) {
    std::cout << "Adding " << el->GetElementName() << " to " << det.name << std::endl;
  }

  // rather than the assembly we add directly the layers
  for (TEveElement::List_i itr2 = el->BeginChildren(); itr2 != el->EndChildren(); itr2++) {
    TEveElement *child = *itr2;
    det.list->AddElement(child);
    child->SetRnrSelfChildren(true, false);
    child->SetMainTransparency(useTransparencies ? 60 : 0);
    if (overrideDetectorColors) child->SetMainColor(kRed);
    ((TEveGeoShape *)child)->SetNSegments(128);
    ((TEveGeoShape *)child)->SetDrawFrame(false);
  }
}

void EventDisplay::addECalB(TEveElement* el, SubDetector& det)
{
  addSimpleElement(el, det);
  el->SetRnrSelfChildren(true, false);
  if (overrideDetectorColors) el->SetMainColor(kAzure+7);
  el->SetMainTransparency(useTransparencies ? 60 : 0);
  ((TEveGeoShape *)el)->SetDrawFrame(false);

  // daughter properties - essentially for projections...
  for (TEveElement::List_i itr2 = el->BeginChildren(); itr2 != el->EndChildren(); itr2++)
  {
    TEveElement *child = *itr2;
    if (overrideDetectorColors) child->SetMainColor(kAzure+7);
    TString elname(child->GetElementName());
    if (debug) std::cout << elname << std::endl;
    child->SetRnrSelfChildren(false, false);
    if (elname.BeginsWith("ECAL_Cryo_side"))
      child->SetRnrSelfChildren(false, false);
    else
      child->SetRnrSelfChildren(true, false);
    child->SetMainTransparency(useTransparencies ? 60 : 0);
    ((TEveGeoShape *)child)->SetNSegments(128);
    ((TEveGeoShape *)child)->SetDrawFrame(false);
  }
}


void EventDisplay::addECalEC(TEveElement* el, SubDetector& det)
{
  addSimpleElement(el, det);
  el->SetRnrSelfChildren(true, false);
  el->SetMainTransparency(useTransparencies ? 60 : 0);
  if (overrideDetectorColors) el->SetMainColor(kAzure+7);
  ((TEveGeoShape *)el)->SetDrawFrame(false);
}

void EventDisplay::addHCalB(TEveElement* el, SubDetector& det)
{
  // for the HCal barrel, we draw the envelope in the 3d model,
  // and the subvolumes in the 2D views

  // add element to list
  addSimpleElement(el, det);

  // for 3D view, draw envelope but not children
  // for projections will be overridden later
  el->SetRnrSelfChildren(true, false);
  // this if we wanted to not draw the envelope but not the layers.. but then too many details in 3D
  // el->SetRnrSelfChildren(false, true);
  el->SetMainTransparency(useTransparencies ? 60 : 0);
  if (overrideDetectorColors) el->SetMainColor(kAzure-7);
  ((TEveGeoShape *)el)->SetNSegments(128);
  ((TEveGeoShape *)el)->SetDrawFrame(false);

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
    // or not
    // ((TEveGeoShape *)child)->SetDrawFrame(false);
    ((TEveGeoShape *)child)->SetNSegments(128);
  }
}

void EventDisplay::addHCalEC(TEveElement* el, SubDetector& det)
{
  // for the HCal endcap, we draw the envelope in the 3d model,
  // and the subvolumes in the 2D views

  // add element to list
  addSimpleElement(el, det);

  // for 3D view, draw envelope but not children
  // for projections will be overridden later
  el->SetRnrSelfChildren(true, false);
  // this if we wanted to not draw the envelope but not the layers.. but then too many details in 3D
  // el->SetRnrSelfChildren(false, true);
  el->SetMainTransparency(useTransparencies ? 60 : 0);
  if (overrideDetectorColors) el->SetMainColor(kAzure-7);
  ((TEveGeoShape *)el)->SetNSegments(128);
  ((TEveGeoShape *)el)->SetDrawFrame(false);

  // setup the properties of the children of the envelope
  for (TEveElement::List_i itr = el->BeginChildren(); itr != el->EndChildren(); itr++)
  {
    TEveElement *child = *itr;
    TString elname(child->GetElementName());
    if (elname.Contains("ECLayer")) {
      // do not show the layers
      // child->SetRnrSelfChildren(true, false);
      // show the layers
      child->SetRnrSelfChildren(false, true);
      child->SetMainColor(kAzure-7);
      // child->SetMainTransparency(useTransparencies ? 60 : 0);
      child->SetMainTransparency(0);
      // set to false if we don't want to see the layer structure
      ((TEveGeoShape *)child)->SetDrawFrame(true);
      // ((TEveGeoShape *)child)->SetDrawFrame(false);
    }
  }

  // quite identical to HCALB only difference for B is ((TEveGeoShape *)child)->SetNSegments(128)... maybe we can merge code...
}

void EventDisplay::addMuonB(TEveElement* el, SubDetector& det)
{
  // add element to list
  addSimpleElement(el, det);

  // draw envelope but not children
  el->SetRnrSelfChildren(true, false);
  el->SetMainTransparency(useTransparencies ? 60 : 0);
  if (overrideDetectorColors) el->SetMainColor(kOrange);
  ((TEveGeoShape *)el)->SetNSegments(128);
  ((TEveGeoShape *)el)->SetDrawFrame(false);

  // setup the properties of the children of the envelope
  // for the projections (number of segments for rho-phi view)
  for (TEveElement::List_i itr = el->BeginChildren(); itr != el->EndChildren(); itr++)
  {
    TEveElement *child = *itr;
    // show the radial layers
    // child->SetRnrSelfChildren(false, true);
    ((TEveGeoShape *)child)->SetNSegments(128);
    // child->SetMainColor(kOrange);
    ((TEveGeoShape *)child)->SetDrawFrame(true);
  }
}

void EventDisplay::addMuonEC(TEveElement* el, SubDetector& det)
{
  // add element to list
  addSimpleElement(el, det);

  // draw envelope but not children
  el->SetRnrSelfChildren(true, false);
  el->SetMainTransparency(useTransparencies ? 60 : 0);
  if (overrideDetectorColors) el->SetMainColor(kOrange);
  ((TEveGeoShape *)el)->SetDrawFrame(false);

  // setup the properties of the children of the envelope
  // do not do it here, but directly in the projections..
  /*
  for (TEveElement::List_i itr = el->BeginChildren(); itr != el->EndChildren(); itr++)
  {
    TEveElement *child = *itr;
    // show the radial layers
    child->SetRnrSelfChildren(false, true);
    ((TEveGeoShape *)child)->SetNSegments(128);
    child->SetMainColor(kOrange);
    ((TEveGeoShape *)child)->SetDrawFrame(true);
  }
  */

  // identical to HCALB only difference for B is ((TEveGeoShape *)child)->SetNSegments(128)... maybe we can merge code...
}
