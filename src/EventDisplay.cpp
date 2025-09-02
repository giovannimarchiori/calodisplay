/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventDisplay: creates and populates the event display
//
/******************************************************************************/

#include "EventDisplay.h"
#include "DetectorGeometry.h"
#include "MagField.h"
#include "Globals.h"

#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPRegexp.h>

#include <TGButton.h>
#include <TGTab.h>

#include <TGeoTube.h>
#include <TGeoMatrix.h>

#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEveProjectionAxes.h>
#include <TTimeStamp.h>

#include <filesystem>
#include <unordered_map>
#include <iostream>

const bool debug = true;

// return the sign of a float
int sgn(float val)
{
  return (val > 0.) - (val < 0.);
}

// return solution of 2nd order poly (ax^2+bx+c = 0)
double solve_poly2(double a, double b, double c, int sol) {
  if (a==0.) return -c/b;
  double d = b*b-4*a*c;
  if (d<0) return -9999999999.;
  d = std::sqrt(d);
  if (sol>0)
    return (-b + d)/(2*a);
  else
    return (-b - d)/(2*a);
}

const char* partType(int pdgID) {
  if (pdgID == 111)
    return "pi0";
  else if (pdgID == 211)
    return "pi+";
  else if (pdgID == -211)
    return "pi-";
  else if (pdgID == 22)
    return "y";
  else if (pdgID == 11)
    return "e-";
  else if (pdgID == -11)
    return "e+";
  else if (pdgID == 13)
    return "mu-";
  else if (pdgID == -13)
    return "mu+";
  else
    return "unknown";
}

/******************************************************************************/
// HELPER FUNCTIONS related to the graphic system
/******************************************************************************/

// derived TGLAnnotation class that overrides MouseEnter method
// to avoid editing of annotation
TGLConstAnnotation::TGLConstAnnotation(TGLViewerBase *parent, const char *text, Float_t posx, Float_t posy) : TGLAnnotation(parent, text, posx, posy)
{
  ;
}

Bool_t TGLConstAnnotation::MouseEnter(TGLOvlSelectRecord & /*rec*/)
{
  fActive = kFALSE;
  return kTRUE;
}

void EventDisplay::FillClusters(std::string clusterType)
{
  std::vector<CaloCluster *> *clusterData = nullptr;
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
    for (CaloCluster *cluster : *clusterData)
    {
      delete cluster;
    }
    clusterData->resize(0);
  }

  // loop over clusters and fill the relevant info
  const double cm = geomReader->cm;
  const double mm = geomReader->mm;
  const int nECalLayers = geomReader->nLayers;
  const int nHCalLayers = geomReader->nLayersHCal;
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
    CaloCluster *cluster = new CaloCluster;
    clusterData->emplace_back(cluster);

    // set index and energy
    cluster->setIndex(i);
    cluster->setEnergy(E);

    // set barycenter
    if (debug)
    {
      std::cout << "  Setting cluster barycenter position .." << std::endl;
    }
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

    // retrieve index of first and last+1 cell in cell vectors
    if (debug)
    {
      std::cout << "  Retrieving cluster cells and calculating energy per layer .." << std::endl;
    }
    unsigned int first_hit, last_hit;
    unsigned int nCells;
    if (clusterType == "topo")
    {
      nCells = eventReader->CaloTopoClusterCells_energy->GetSize();
      first_hit = (*eventReader->CaloTopoClusters_hits_begin)[i];
      last_hit = (*eventReader->CaloTopoClusters_hits_end)[i];
    }
    else
    {
      nCells = eventReader->CaloClusterCells_energy->GetSize();
      first_hit = (*eventReader->CaloClusters_hits_begin)[i];
      last_hit = (*eventReader->CaloClusters_hits_end)[i];
    }

    // loop over cells to calculate energy per layer
    std::vector<float> energyVsECalLayer(nECalLayers, 0.);
    std::vector<float> energyVsHCalLayer(nHCalLayers, 0.);
    std::vector<float> energyVsMuonLayer(nMuonLayers, 0.);
    for (unsigned int iCell = 0; iCell < nCells; iCell++)
    {
      if (iCell < first_hit || iCell >= last_hit)
        continue;
      ULong_t cellID;
      float energy, x_center, y_center, z_center;
      if (clusterType == "topo")
      {
        cellID = (*eventReader->CaloTopoClusterCells_cellID)[iCell];
        energy = (*eventReader->CaloTopoClusterCells_energy)[iCell];
      }
      else
      {
        cellID = (*eventReader->CaloClusterCells_cellID)[iCell];
        energy = (*eventReader->CaloClusterCells_energy)[iCell];
      }
      ULong_t systemID = geomReader->SystemID(cellID);
      ULong_t layer;
      if (systemID == 4)
      {
        layer = geomReader->ECalBarrelLayer(cellID);
        energyVsECalLayer[layer] += energy;
      }
      else if (systemID == 8)
      {
        layer = geomReader->HCalBarrelLayer(cellID);
        energyVsHCalLayer[layer] += energy;
      }
      else if (systemID == 5)
      {
	// TODO: calculate properly ecal endcap layer
	// layer = geomReader->ECalEndcapLayer(cellID);
        // energyVsECalLayer[layer] += energy;
	energyVsECalLayer[0] += energy;
      }
      else if (systemID == 9)
      {
        // TODO: calculate properly hcal endcap layer
        // layer = geomReader->HCalBarrelLayer(cellID);
        //energyVsHCalLayer[layer] += energy;
	energyVsHCalLayer[0] += energy;
      }
      else if (systemID == 12)
      {
        // TODO: calculate properly muon barrel layer
	energyVsMuonLayer[0] += energy;
      }
      else if (systemID == 13)
      {
        // TODO: calculate properly muon endcap layer
	energyVsMuonLayer[0] += energy;
      }
      else
      {
        std::cout << "Unknown system " << systemID << std::endl;
        exit(1);
      }
    }
    cluster->setEnergyVsECalLayers(energyVsECalLayer);
    cluster->setEnergyVsHCalLayers(energyVsHCalLayer);
    cluster->setEnergyVsMuonLayers(energyVsMuonLayer);

    // loop again over cells to calculate barycenter per layer
    if (debug)
    {
      std::cout << "  Calculating cell barycenter per layer .." << std::endl;
    }
    std::vector<TVector3> barycenterVsECalLayer(nECalLayers);
    std::vector<float> sumWeightsVsECalLayer(nECalLayers, 0.0);
    std::vector<TVector3> barycenterVsHCalLayer(nHCalLayers);
    std::vector<float> sumWeightsVsHCalLayer(nHCalLayers, 0.0);
    std::vector<TVector3> barycenterVsMuonLayer(nMuonLayers);
    std::vector<float> sumWeightsVsMuonLayer(nMuonLayers, 0.0);
    for (unsigned int iCell = 0; iCell < nCells; iCell++)
    {
      if (iCell < first_hit || iCell >= last_hit)
        continue;
      ULong_t cellID;
      float energy, x, y, z;
      if (clusterType == "topo")
      {
        cellID = (*eventReader->CaloTopoClusterCells_cellID)[iCell];
        energy = (*eventReader->CaloTopoClusterCells_energy)[iCell];
        x = (*eventReader->CaloTopoClusterCells_position_x)[iCell] * mm;
        y = (*eventReader->CaloTopoClusterCells_position_y)[iCell] * mm;
        z = (*eventReader->CaloTopoClusterCells_position_z)[iCell] * mm;
      }
      else
      {
        cellID = (*eventReader->CaloClusterCells_cellID)[iCell];
        energy = (*eventReader->CaloClusterCells_energy)[iCell];
        x = (*eventReader->CaloClusterCells_position_x)[iCell] * mm;
        y = (*eventReader->CaloClusterCells_position_y)[iCell] * mm;
        z = (*eventReader->CaloClusterCells_position_z)[iCell] * mm;
      }
      ULong_t systemID = geomReader->SystemID(cellID);
      ULong_t layer;
      if (systemID == 4)
      {
        layer = geomReader->ECalBarrelLayer(cellID);
        if (energyVsECalLayer[layer] > 0)
        {
          float weight = energy / energyVsECalLayer[layer];
          barycenterVsECalLayer[layer].SetXYZ(
              barycenterVsECalLayer[layer].X() + x * weight,
              barycenterVsECalLayer[layer].Y() + y * weight,
              barycenterVsECalLayer[layer].Z() + z * weight);
          sumWeightsVsECalLayer[layer] += weight;
        }
      }
      else if (systemID == 8)
      {
        layer = geomReader->HCalBarrelLayer(cellID);
        if (energyVsHCalLayer[layer] > 0)
        {
          float weight = energy / energyVsHCalLayer[layer];
          barycenterVsHCalLayer[layer].SetXYZ(
              barycenterVsHCalLayer[layer].X() + x * weight,
              barycenterVsHCalLayer[layer].Y() + y * weight,
              barycenterVsHCalLayer[layer].Z() + z * weight);
          sumWeightsVsHCalLayer[layer] += weight;
        }
      }
    }
    for (unsigned int iLayer = 0; iLayer < nECalLayers; iLayer++)
    {
      if (energyVsECalLayer[iLayer] > 0)
      {
        barycenterVsECalLayer[iLayer].SetXYZ(
            barycenterVsECalLayer[iLayer].X() / sumWeightsVsECalLayer[iLayer],
            barycenterVsECalLayer[iLayer].Y() / sumWeightsVsECalLayer[iLayer],
            barycenterVsECalLayer[iLayer].Z() / sumWeightsVsECalLayer[iLayer]);
      }
    }
    cluster->setBarycenterVsECalLayers(barycenterVsECalLayer);

    for (unsigned int iLayer = 0; iLayer < nHCalLayers; iLayer++)
    {
      if (energyVsHCalLayer[iLayer] > 0.)
      {
        barycenterVsHCalLayer[iLayer].SetXYZ(
            barycenterVsHCalLayer[iLayer].X() / sumWeightsVsHCalLayer[iLayer],
            barycenterVsHCalLayer[iLayer].Y() / sumWeightsVsHCalLayer[iLayer],
            barycenterVsHCalLayer[iLayer].Z() / sumWeightsVsHCalLayer[iLayer]);
      }
    }
    cluster->setBarycenterVsHCalLayers(barycenterVsHCalLayer);

    cluster->setBarycenterVsMuonLayers(barycenterVsMuonLayer);

    cluster->print();
  }
}

void EventDisplay::DrawClusters(std::string clusterType)
{
  std::vector<CaloCluster *> *clusterData = nullptr;
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

  const double cm = geomReader->cm;
  const double mm = geomReader->mm;
  const double rMin = geomReader->rMin;
  const double alpha = geomReader->alpha;
  const double thetaGrid = geomReader->thetaGrid;
  const double gridPhi = geomReader->gridPhi;
  const double phiMin = geomReader->phiMin;
  const double thetaGridHCal = geomReader->thetaGridHCal;
  const double gridPhiHCal = geomReader->gridPhiHCal;
  const double phiMinHCal = geomReader->phiMinHCal;

  // create the graphic containers
  TEveRGBAPalette *pal = new TEveRGBAPalette(0, 1000);

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

    if (displayConfig.getBoolConfig("drawClusterDirection")) {
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

    if (displayConfig.getBoolConfig("drawClusterDirection")) {
      if (swclustersDirection == nullptr) {
	swclustersDirection = new TEveStraightLineSet("cluster directions");
	swclusters_3D->AddElement(swclustersDirection);
      }
      else
	swclustersDirection->DestroyElements();
    }
    clustersDirection = swclustersDirection;

  }
  if (displayConfig.getBoolConfig("drawClusterDirection")) {
    clustersCenter->SetMarkerStyle(4);
    clustersCenter->SetMarkerSize(5);
    clustersDirection->SetLineColor(kWhite);
    clustersDirection->SetLineWidth(5);
  }
  
  for (unsigned int i = 0; i < nClusters; i++)
  {
    float E = (clusterType == "topo") ? (*eventReader->CaloTopoClusters_energy)[i] : (*eventReader->CaloClusters_energy)[i];
    if (E < ClusterEnergyThreshold)
      continue;
    // topo cluster positions are in cm while calo clusters and hits/cells in mm ... ?
    if (clusterType == "topo")
    {
      float x0 = (*eventReader->CaloTopoClusters_position_x)[i] * mm;
      float y0 = (*eventReader->CaloTopoClusters_position_y)[i] * mm;
      float z0 = (*eventReader->CaloTopoClusters_position_z)[i] * mm;
      clustersCenter->SetNextPoint(x0, y0, z0);

      if (displayConfig.getBoolConfig("drawClusterDirection")) {
	float theta = (*eventReader->CaloTopoClusters_theta)[i];
	float phi = (*eventReader->CaloTopoClusters_phi)[i];
	float a = std::sin(theta)*std::sin(theta);
	float b = 2*std::sin(theta)*(x0*std::cos(phi) + y0*std::sin(phi));
	float c = x0*x0 + y0*y0 - rMin*rMin;
	float t1 = solve_poly2(a,b,c,1);
	double rMax = doHCal ? geomReader->rMaxHCal : geomReader->rMax;
	c = x0*x0 + y0*y0 - rMax*rMax;
	float t2 = solve_poly2(a,b,c,1);
	clustersDirection->AddLine(x0 + t1*std::sin(theta)*std::cos(phi),
				   y0 + t1*std::sin(theta)*std::sin(phi),
				   z0 + t1*std::cos(theta),
				   x0 + t2*std::sin(theta)*std::cos(phi),
				   y0 + t2*std::sin(theta)*std::sin(phi),
				   z0 + t2*std::cos(theta));
      }
    }
    else
    {
      float x0 = (*eventReader->CaloClusters_position_x)[i] * mm;
      float y0 = (*eventReader->CaloClusters_position_y)[i] * mm;
      float z0 = (*eventReader->CaloClusters_position_z)[i] * mm;
      clustersCenter->SetNextPoint(x0, y0, z0);

      if (displayConfig.getBoolConfig("drawClusterDirection")) {
	float theta = (*eventReader->CaloClusters_theta)[i];
	float phi = (*eventReader->CaloClusters_phi)[i];
	float a = std::sin(theta)*std::sin(theta);
	float b = 2*std::sin(theta)*(x0*std::cos(phi) + y0*std::sin(phi));
	float c = x0*x0 + y0*y0 - rMin*rMin;
	float t1 = solve_poly2(a,b,c,1);
	double rMax = doHCal ? geomReader->rMaxHCal : geomReader->rMax;
	c = x0*x0 + y0*y0 - rMax*rMax;
	float t2 = solve_poly2(a,b,c,1);
	clustersDirection->AddLine(x0 + t1*std::sin(theta)*std::cos(phi),
				   y0 + t1*std::sin(theta)*std::sin(phi),
				   z0 + t1*std::cos(theta),
				   x0 + t2*std::sin(theta)*std::cos(phi),
				   y0 + t2*std::sin(theta)*std::sin(phi),
				   z0 + t2*std::cos(theta));
      }
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
  int icl = 0;
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
        CaloCluster *clus = (*clusterData)[icl];
        for (unsigned int iLayer = 0; iLayer < clus->getNLayersECal(); iLayer++)
        {
          if (clus->getEnergyInECalLayer(iLayer) > 0)
          {
            layerBarycenters->SetNextPoint(
                                           clus->getBarycenterInECalLayer(iLayer).X() / cm,
                                           clus->getBarycenterInECalLayer(iLayer).Y() / cm,
                                           clus->getBarycenterInECalLayer(iLayer).Z() / cm);
          }
        }
        for (unsigned int iLayer = 0; iLayer < (*clusterData)[icl]->getNLayersHCal(); iLayer++)
        {
          if (clus->getEnergyInHCalLayer(iLayer) > 0)
          {
            layerBarycenters->SetNextPoint(
                                           clus->getBarycenterInHCalLayer(iLayer).X() / cm,
                                           clus->getBarycenterInHCalLayer(iLayer).Y() / cm,
                                           clus->getBarycenterInHCalLayer(iLayer).Z() / cm);
          }
        }
      }
      icl++;
    }
  }
  // loop over cells
  std::unordered_map<ULong64_t, double> cellEnergies_rhoz;
  std::unordered_map<ULong64_t, double> cellEnergies_rhophi;
  unsigned int nCells = (clusterType == "topo") ? eventReader->CaloTopoClusterCells_energy->GetSize() : eventReader->CaloClusterCells_energy->GetSize();
  for (unsigned int i = 0; i < nCells; i++)
  {
    int icl = -1;
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
      // TODO check again if < or <=
      if (i >= first_hit && i < last_hit)
      {
        icl = j;
        break;
      }
    }
    if (icl == -1)
      continue; // should never happen..
    ULong_t cellID;
    float energy, x_center, y_center, z_center;
    if (clusterType == "topo")
    {
      cellID = (*eventReader->CaloTopoClusterCells_cellID)[i];
      energy = (*eventReader->CaloTopoClusterCells_energy)[i];
      x_center = (*eventReader->CaloTopoClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->CaloTopoClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->CaloTopoClusterCells_position_z)[i] * mm;
    }
    else
    {
      cellID = (*eventReader->CaloClusterCells_cellID)[i];
      energy = (*eventReader->CaloClusterCells_energy)[i];
      x_center = (*eventReader->CaloClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->CaloClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->CaloClusterCells_position_z)[i] * mm;
    }
    float r_center = sqrt(x_center * x_center + y_center * y_center);
    int systemID = (int)DetectorGeometry::SystemID(cellID);
    // for 2D projections, assign unique ID
    // based on
    // ECAL: layer-theta for rho-z, layer-module for rho-phi
    // HCAL: layer-theta for rho-z, layer-phi for rho-phi
    // layer = 0..11 fits in 4 bits (0..15)
    // theta = 0..800 fits in 10 bits (0..1023)
    // module = 0..1535 fits in 11 bits (0..2043)
    // cluster = assume < 512 (9 bits)
    // so both IDs can fit in a 32-bit integer
    // but could even play it safer and just use cellID and remove unused fields
    float r_in, r_out;
    int layer, thetaID, moduleID, phiID;
    ULong64_t rhophiID(0), rhozID(0);

    // temporary
    // rhophiID = cellID;
    // rhozID = cellID;

    if (systemID == 4)
    {
      layer = (int)DetectorGeometry::ECalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::ECalBarrelThetaBin(cellID);
      moduleID = (int)DetectorGeometry::ECalBarrelModule(cellID);
      rhophiID = systemID + 16 * (icl + 512 * layer + 512 * 16 * moduleID);
      rhozID = systemID + 16 * (icl + 512 * layer + 512 * 16 * thetaID);
      r_in = geomReader->r[layer];
      r_out = geomReader->r[layer + 1];
    }
    else if (systemID == 8)
    {
      layer = (int)DetectorGeometry::HCalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::HCalBarrelThetaBin(cellID);
      phiID = (int)DetectorGeometry::HCalBarrelPhiBin(cellID);
      rhophiID = systemID + 16 * (icl + 512 * layer + 512 * 16 * phiID);
      rhozID = systemID + 16 * (icl + 512 * layer + 512 * 16 * thetaID);
      r_in = geomReader->rHCal[layer];
      r_out = geomReader->rHCal[layer + 1];
    }
    else if (systemID == 5)
    {
      // std::cout << "Not drawing cells for ECAL endcap: system " << systemID << std::endl;
      ;
    }
    else if (systemID == 9)
    {
      // std::cout << "Not drawing cells for HCAL endcap: system " << systemID << std::endl;
      ;
    }
    else if (systemID == 12)
    {
      // std::cout << "Not drawing cells for MUON barrel: system " << systemID << std::endl;
      ;
    }
    else if (systemID == 13)
    {
      // std::cout << "Not drawing cells for MUON barrel: system " << systemID << std::endl;
      ;
    }
    else
    {
      std::cout << "Unknown system " << systemID << std::endl;
      exit(1);
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
      if (systemID == 4)
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
      else if (systemID==8)
      {
        // HCAL
        float phiMin = phiMinHCal + gridPhiHCal * phiID;
        float phiMax = phiMinHCal + gridPhiHCal * (phiID + 1);
        float thetaMin = theta_center - thetaGridHCal / 2.;
        float thetaMax = theta_center + thetaGridHCal / 2.;

        verts3D[0] = r_in * cos(phiMin);
        verts3D[1] = r_in * sin(phiMin);
        verts3D[2] = r_in / tan(thetaMin);

        verts3D[3] = r_in * cos(phiMax);
        verts3D[4] = r_in * sin(phiMax);
        verts3D[5] = r_in / tan(thetaMin);

        verts3D[6] = r_in * cos(phiMax);
        verts3D[7] = r_in * sin(phiMax);
        verts3D[8] = r_in / tan(thetaMax);

        verts3D[9] = r_in * cos(phiMin);
        verts3D[10] = r_in * sin(phiMin);
        verts3D[11] = r_in / tan(thetaMax);

        verts3D[12] = r_out * cos(phiMin);
        verts3D[13] = r_out * sin(phiMin);
        verts3D[14] = r_out / tan(thetaMin);

        verts3D[15] = r_out * cos(phiMax);
        verts3D[16] = r_out * sin(phiMax);
        verts3D[17] = r_out / tan(thetaMin);

        verts3D[18] = r_out * cos(phiMax);
        verts3D[19] = r_out * sin(phiMax);
        verts3D[20] = r_out / tan(thetaMax);

        verts3D[21] = r_out * cos(phiMin);
        verts3D[22] = r_out * sin(phiMin);
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
      x_center = (*eventReader->CaloTopoClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->CaloTopoClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->CaloTopoClusterCells_position_z)[i] * mm;
    }
    else {
      x_center = (*eventReader->CaloClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->CaloClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->CaloClusterCells_position_z)[i] * mm;
    }
    float r_center = sqrt(x_center*x_center + y_center*y_center);

    float theta_center = atan2(r_center, z_center);
    //float phi_center = atan2(y_center, x_center);

    float verts[12];
    verts[0] = r_in / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
    verts[1] = r_in * sgn(y_center);
    verts[3] = r_in / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
    verts[4] = r_in * sgn(y_center);
    verts[6] = r_out / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
    verts[7] = r_out * sgn(y_center);
    verts[9] = r_out / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
    verts[10] = r_out * sgn(y_center);
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
    int icl = -1;
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
      // TODO check again if < or <=
      if (i >= first_hit && i < last_hit)
      {
        icl = j;
        break;
      }
    }
    if (icl == -1)
      continue; // should never happen..
    ULong_t cellID;
    float energy;
    if (clusterType == "topo")
    {
      cellID = (*eventReader->CaloTopoClusterCells_cellID)[i];
    }
    else
    {
      cellID = (*eventReader->CaloClusterCells_cellID)[i];
    }

    int system = (int)DetectorGeometry::SystemID(cellID);
    int layer, thetaID, moduleID;
    long phiID(0), rhozID(0);
    float r_in, r_out, thetaMin, thetaMax;
    if (system == 4)
    {
      layer = (int)DetectorGeometry::ECalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::ECalBarrelThetaBin(cellID);
      rhozID = system + 16 * (icl + 512 * layer + 512 * 16 * thetaID);
      r_in = geomReader->r[layer];
      r_out = geomReader->r[layer + 1];
    }
    else if (system == 8)
    {
      layer = (int)DetectorGeometry::HCalBarrelLayer(cellID);
      thetaID = (int)DetectorGeometry::HCalBarrelThetaBin(cellID);
      rhozID = system + 16 * (icl + 512 * layer + 512 * 16 * thetaID);
      r_in = geomReader->rHCal[layer];
      r_out = geomReader->rHCal[layer + 1];
    }
    if (cellEnergies_rhoz.count(rhozID) == 0)
      continue;
    energy = cellEnergies_rhoz[rhozID];
    cellEnergies_rhoz.erase(rhozID);
    float x_center, y_center, z_center;
    if (clusterType == "topo")
    {
      x_center = (*eventReader->CaloTopoClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->CaloTopoClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->CaloTopoClusterCells_position_z)[i] * mm;
    }
    else
    {
      x_center = (*eventReader->CaloClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->CaloClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->CaloClusterCells_position_z)[i] * mm;
    }
    float r_center = sqrt(x_center * x_center + y_center * y_center);
    float theta_center = atan2(r_center, z_center);
    if (system == 4)
    {
      thetaMin = theta_center - thetaGrid * geomReader->mergedCells_Theta[layer] / 2.;
      thetaMax = theta_center + thetaGrid * geomReader->mergedCells_Theta[layer] / 2.;
    }
    else if (system == 8)
    {
      thetaMin = theta_center - thetaGridHCal / 2.;
      thetaMax = theta_center + thetaGridHCal / 2.;
    }
    else
      continue;

    verts[0] = r_in / tan(thetaMin);
    verts[1] = r_in * sgn(y_center);
    verts[3] = r_in / tan(thetaMax);
    verts[4] = r_in * sgn(y_center);
    verts[6] = r_out / tan(thetaMax);
    verts[7] = r_out * sgn(y_center);
    verts[9] = r_out / tan(thetaMin);
    verts[10] = r_out * sgn(y_center);
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
    // continue;  // GM FIXME!!!! WHY SYSTEM IS 8????
    ULong64_t ikey = (ULong64_t)key;
    int system = ikey % 16;
    ikey = ikey / 16;
    int icl = ikey % 512;
    if (system == 4)
    {
      int layer = (ikey / 512) % 16;
      int moduleID = ikey / (512 * 16);
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
    else if (system==8)
    {
      int layer = (ikey / 512) % 16;
      int phiID = ikey / (512 * 16);
      // std::cout << "cluster layer module = " << icl << " " << layer << " " << moduleID << std::endl;
      float r_in = geomReader->rHCal[layer];
      float r_out = geomReader->rHCal[layer + 1];
      float phiMin = phiMinHCal + gridPhiHCal * phiID;
      float phiMax = phiMinHCal + gridPhiHCal * (phiID + 1);
      verts[0] = r_in * cos(phiMin);
      verts[1] = r_in * sin(phiMin);
      verts[3] = r_in * cos(phiMax);
      verts[4] = r_in * sin(phiMax);
      verts[6] = r_out * cos(phiMax);
      verts[7] = r_out * sin(phiMax);
      verts[9] = r_out * cos(phiMin);
      verts[10] = r_out * sin(phiMin);
      verts[2] = verts[5] = verts[8] = verts[11] = 0.;
      qs_rhophi[icl]->AddQuad(verts);
      qs_rhophi[icl]->QuadValue((int)(1000 * energy));
      qs_rhophi[icl]->QuadId(new TNamed(Form("%s %d cell L%d P%d", clusterType.c_str(), icl, layer, phiID), "Dong!"));
    }
  }
}

void EventDisplay::loadEvent(int event)
{

  if (event != -1)
    eventId = event;

  printf("Loading event %d ...\n", eventId);
  textEntry->SetTextColor(0xff0000);
  textEntry->SetText(Form("Loading event %d ...", eventId));

  eventReader->loadEvent(eventId);

  // find out pdgID and momentum of primary particle
  double pmax = 0.0;
  int ipmax = -1;
  for (unsigned int i = 0; i < eventReader->genParticles_generatorStatus->GetSize(); i++)
  {
    if ((*eventReader->genParticles_generatorStatus)[i] != 1)
      continue;
    float _px = (*eventReader->genParticles_momentum_x)[i];
    float _py = (*eventReader->genParticles_momentum_y)[i];
    float _pz = (*eventReader->genParticles_momentum_z)[i];
    float _p = sqrt(_px * _px + _py * _py + _pz * _pz);
    // cout << p << endl;
    if (_p == 0.)
      continue;
    if (_p > pmax)
    {
      pmax = _p;
      ipmax = i;
    }
  }
  int pdgID = (*eventReader->genParticles_PDG)[ipmax];

  const double cm = geomReader->cm;
  const double mm = geomReader->mm;
  double rMax = 5000.*mm;
  if (not showFullDetector) {
    rMax = doHCal ? geomReader->rMaxHCal : geomReader->rMax;
  }
  /*
  const double rMin = geomReader->rMin;
  const double alpha = geomReader->alpha;
  const double thetaGrid = geomReader->thetaGrid;
  const double gridPhi = geomReader->gridPhi;
  const double phiMin = geomReader->phiMin;
  */

  //
  // particles
  //
  
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
      //trkProp->SetMagField(-2.0); // tesla
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

    // handle differently the e/gamma vs pi0 particle guns
    // (for pi0, need to look for the two photons in secondary particle list
    // unfortunately cannot just use the latter to show particles since
    // info about endpoint is broken (identical to vertex)
    // TODO: UPDATE NOW ENDPOINT IS FILLED!

    for (unsigned int ip = 0; ip < eventReader->genParticles_generatorStatus->GetSize(); ip++)
    {
      float m = (*eventReader->genParticles_mass)[ip];
      float px = (*eventReader->genParticles_momentum_x)[ip];
      float py = (*eventReader->genParticles_momentum_y)[ip];
      float pz = (*eventReader->genParticles_momentum_z)[ip];
      float p = sqrt(px * px + py * py + pz * pz);
      if (p < ParticleEnergyThreshold)
	continue;
      float pT = sqrt(px * px + py * py);

      double t = (*eventReader->genParticles_time)[ip];
      double x1 = (*eventReader->genParticles_vertex_x)[ip] * mm;
      double y1 = (*eventReader->genParticles_vertex_y)[ip] * mm;
      double z1 = (*eventReader->genParticles_vertex_z)[ip] * mm;

      TEveMCTrack mct;
      mct.SetPdgCode((*eventReader->genParticles_PDG)[ip]);
      mct.SetMomentum(px, py, pz, sqrt(p * p + m * m));
      mct.SetProductionVertex(x1, y1, z1, t);
      TEveTrack *track = new TEveTrack(&mct, particles->GetPropagator());
      track->SetAttLineAttMarker(particles);
      track->SetElementTitle(Form("p = %.3f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
                                  p, acos(pz / p), atan2(py, px),
                                x1 / cm, y1 / cm, z1 / cm));
      TEveVectorF v;
      v[0] = (*eventReader->genParticles_endpoint_x)[ip] * mm;
      v[1] = (*eventReader->genParticles_endpoint_y)[ip] * mm;
      v[2] = (*eventReader->genParticles_endpoint_z)[ip] * mm;
      // cout << v[0] << " " << v[1] << " " << v[2] << endl;
      if ((*eventReader->genParticles_PDG)[ip]==22 ||
          (*eventReader->genParticles_PDG)[ip]==111
          ) {
	// don't do it for charged particles otherwise will
	// force the particle to pass through it and if there
	// are not enough other reference points along the
	// trajectory it will screw up the helix
        TEvePathMark mark(TEvePathMark::kDecay, v);
        track->AddPathMark(mark);
	track->SetRnrPoints(true);
      }
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
        // cout << "p = "<< p << endl;
        if (p < ParticleEnergyThreshold)
          continue;

        // cout << "PDG = "<< pdgID << endl;
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
        // cout << "x1 y1 z1 x2 y2 z2 = "
        //      << x1 << " " << y1 << " " << z1 << " "
        //      << x2 << " " << y2 << " " << z2 << endl;
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
    particles->MakeTracks();
  }

  //
  // hits (VTX/DCH/Wrapper)
  //
  if (displayConfig.getBoolConfig("drawVertexHits"))
  {
    if (hits == nullptr)
    {
      hits = new TEveElementList("hits");
      gEve->AddElement(hits);
    }
    // do we need to Reset() otherwise ?
    std::cout << "Creating vertex hits" << std::endl;
    if (vtxHits == nullptr)
    {
      vtxHits = new TEvePointSet();
      // vtxHits->SetName(Form("VTX hits (E>%.2f GeV)", HitEnergyThreshold));
      vtxHits->SetName("VTX hits");
      vtxHits->SetMarkerStyle(4);
      vtxHits->SetMarkerSize(1.6);
      vtxHits->SetMarkerColor(kRed);
      // gEve->AddElement(vtxHits);
      hits->AddElement(vtxHits);
    }
    else
      vtxHits->Reset();
    for (unsigned int i = 0; i < eventReader->VertexBarrelHits_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->VertexBarrelHits_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->VertexBarrelHits_cellID)[i];
      vtxHits->SetNextPoint(
                            (*eventReader->VertexBarrelHits_position_x)[i] * mm,
                            (*eventReader->VertexBarrelHits_position_y)[i] * mm,
                            (*eventReader->VertexBarrelHits_position_z)[i] * mm);
    }
    for (unsigned int i = 0; i < eventReader->VertexEndcapHits_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->VertexEndcapHits_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->VertexEndcapHits_cellID)[i];
      vtxHits->SetNextPoint(
                            (*eventReader->VertexEndcapHits_position_x)[i] * mm,
                            (*eventReader->VertexEndcapHits_position_y)[i] * mm,
                            (*eventReader->VertexEndcapHits_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawDriftChamberHits"))
  {
    if (hits == nullptr)
    {
      hits = new TEveElementList("hits");
      gEve->AddElement(hits);
    }
    // do we need to Reset() otherwise ?
    std::cout << "Creating drift chamber hits" << std::endl;
    if (dchHits == nullptr)
    {
      dchHits = new TEvePointSet();
      // dchHits->SetName(Form("DCH hits (E>%.2f GeV)", HitEnergyThreshold));
      dchHits->SetName("DCH hits");
      dchHits->SetMarkerStyle(4);
      dchHits->SetMarkerSize(1.6);
      dchHits->SetMarkerColor(kRed);
      // gEve->AddElement(dchHits);
      hits->AddElement(dchHits);
    }
    else
      dchHits->Reset();
    for (unsigned int i = 0; i < eventReader->DriftChamberHits_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->DriftChamberHits_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->DriftChamberHits_cellID)[i];
      dchHits->SetNextPoint(
                            (*eventReader->DriftChamberHits_position_x)[i] * mm,
                            (*eventReader->DriftChamberHits_position_y)[i] * mm,
                            (*eventReader->DriftChamberHits_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawSiWrapperHits"))
  {
    if (hits == nullptr)
    {
      hits = new TEveElementList("hits");
      gEve->AddElement(hits);
    }
    std::cout << "Creating wrapper hits" << std::endl;
    if (siwrHits == nullptr)
    {
      siwrHits = new TEvePointSet();
      // siwrHits->SetName(Form("VTX hits (E>%.2f GeV)", HitEnergyThreshold));
      siwrHits->SetName("SiWr hits");
      siwrHits->SetMarkerStyle(4);
      siwrHits->SetMarkerSize(1.6);
      siwrHits->SetMarkerColor(kRed);
      // gEve->AddElement(siwrHits);
      hits->AddElement(siwrHits);
    }
    else
      siwrHits->Reset();
    for (unsigned int i = 0; i < eventReader->SiWrapperBarrelHits_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->SiWrapperBarrelHits_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->SiWrapperBarrelHits_cellID)[i];
      siwrHits->SetNextPoint(
                             (*eventReader->SiWrapperBarrelHits_position_x)[i] * mm,
                             (*eventReader->SiWrapperBarrelHits_position_y)[i] * mm,
                             (*eventReader->SiWrapperBarrelHits_position_z)[i] * mm);
    }
    for (unsigned int i = 0; i < eventReader->SiWrapperEndcapHits_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->SiWrapperEndcapHits_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->SiWrapperEndcapHits_cellID)[i];
      siwrHits->SetNextPoint(
                             (*eventReader->SiWrapperEndcapHits_position_x)[i] * mm,
                             (*eventReader->SiWrapperEndcapHits_position_y)[i] * mm,
                             (*eventReader->SiWrapperEndcapHits_position_z)[i] * mm);
    }
  }

  //
  // hits (ECAL/HCAL)
  //
  if (displayConfig.getBoolConfig("drawECalBarrelHits"))
  {
    if (hits == nullptr)
    {
      hits = new TEveElementList("hits");
      gEve->AddElement(hits);
    }
    std::cout << "Creating ecal barrel hits" << std::endl;
    if (ecalHits == nullptr)
    {
      ecalHits = new TEvePointSet();
      ecalHits->SetName(Form("ECAL hits (E>%.2f GeV)", HitEnergyThreshold));
      ecalHits->SetMarkerStyle(4);
      ecalHits->SetMarkerSize(1);
      ecalHits->SetMarkerColor(kRed);
      // gEve->AddElement(ecalHits);
      hits->AddElement(ecalHits);
    }
    else
      ecalHits->Reset();

    for (unsigned int i = 0; i < eventReader->ECalBarrelHits_position_x->GetSize(); i++)
    {
      float E = (*eventReader->ECalBarrelHits_energy)[i];
      if (E < HitEnergyThreshold)
        continue;
      // ULong_t cellID = (*eventReader->ECalBarrelHits_cellID)[i];
      // ULong_t layer = DetectorGeometry:::Layer(cellID);
      ecalHits->SetNextPoint(
          (*eventReader->ECalBarrelHits_position_x)[i] * mm,
          (*eventReader->ECalBarrelHits_position_y)[i] * mm,
          (*eventReader->ECalBarrelHits_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawECalEndcapHits"))
  {
    if (hits == nullptr)
    {
      hits = new TEveElementList("hits");
      gEve->AddElement(hits);
    }
    std::cout << "Creating ecal endcap hits" << std::endl;
    if (ecalHits == nullptr)
    {
      ecalHits = new TEvePointSet();
      ecalHits->SetName(Form("ECAL hits (E>%.2f GeV)", HitEnergyThreshold));
      ecalHits->SetMarkerStyle(4);
      ecalHits->SetMarkerSize(1);
      ecalHits->SetMarkerColor(kRed);
      hits->AddElement(ecalHits);
    }
    else {
      // dont do it if ecal barrel hits are drawn - this was already done before
      if (not displayConfig.getBoolConfig("drawECalBarrelHits"))
        ecalHits->Reset();
    }
    for (unsigned int i = 0; i < eventReader->ECalEndcapHits_position_x->GetSize(); i++)
    {
      float E = (*eventReader->ECalEndcapHits_energy)[i];
      if (E < HitEnergyThreshold)
        continue;
      ecalHits->SetNextPoint(
          (*eventReader->ECalEndcapHits_position_x)[i] * mm,
          (*eventReader->ECalEndcapHits_position_y)[i] * mm,
          (*eventReader->ECalEndcapHits_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawHCalBarrelHits"))
  {
    if (hits == nullptr)
    {
      hits = new TEveElementList("hits");
      gEve->AddElement(hits);
    }
    std::cout << "Creating hcal barrel hits" << std::endl;
    if (hcalHits == nullptr)
    {
      hcalHits = new TEvePointSet();
      hcalHits->SetName(Form("HCAL hits (E>%.2f GeV)", HitEnergyThreshold));
      hcalHits->SetMarkerStyle(4);
      hcalHits->SetMarkerSize(1);
      hcalHits->SetMarkerColor(kRed);
      hits->AddElement(hcalHits);
    }
    else
      hcalHits->Reset();
    for (unsigned int i = 0; i < eventReader->HCalBarrelHits_position_x->GetSize(); i++)
    {
      float E = (*eventReader->HCalBarrelHits_energy)[i];
      if (E < HitEnergyThreshold)
        continue;
      hcalHits->SetNextPoint(
          (*eventReader->HCalBarrelHits_position_x)[i] * mm,
          (*eventReader->HCalBarrelHits_position_y)[i] * mm,
          (*eventReader->HCalBarrelHits_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawHCalEndcapHits"))
  {
    if (hits == nullptr)
    {
      hits = new TEveElementList("hits");
      gEve->AddElement(hits);
    }
    std::cout << "Creating hcal endcap hits" << std::endl;
    if (hcalHits == nullptr)
    {
      hcalHits = new TEvePointSet();
      hcalHits->SetName(Form("HCAL hits (E>%.2f GeV)", HitEnergyThreshold));
      hcalHits->SetMarkerStyle(4);
      hcalHits->SetMarkerSize(1);
      hcalHits->SetMarkerColor(kRed);
      hits->AddElement(hcalHits);
    }
    else {
      // dont do it if hcal barrel hits are drawn - this was already done before
      if (not displayConfig.getBoolConfig("drawHCalBarrelHits"))
        hcalHits->Reset();
    }
    for (unsigned int i = 0; i < eventReader->HCalEndcapHits_position_x->GetSize(); i++)
    {
      float E = (*eventReader->HCalEndcapHits_energy)[i];
      if (E < HitEnergyThreshold)
        continue;
      hcalHits->SetNextPoint(
          (*eventReader->HCalEndcapHits_position_x)[i] * mm,
          (*eventReader->HCalEndcapHits_position_y)[i] * mm,
          (*eventReader->HCalEndcapHits_position_z)[i] * mm);
    }
  }

  //
  // hits (muon)
  //
  if (displayConfig.getBoolConfig("drawMuonHits"))
  {
    if (hits == nullptr)
    {
      hits = new TEveElementList("hits");
      gEve->AddElement(hits);
    }
    // do we need to Reset() otherwise ?
    std::cout << "Creating muon hits" << std::endl;
    if (muonHits == nullptr)
    {
      muonHits = new TEvePointSet();
      // muonHits->SetName(Form("MUON hits (E>%.2f GeV)", HitEnergyThreshold));
      muonHits->SetName("MUON hits");
      muonHits->SetMarkerStyle(4);
      muonHits->SetMarkerSize(1.6);
      muonHits->SetMarkerColor(kRed);
      // gEve->AddElement(muonHits);
      hits->AddElement(muonHits);
    }
    else {
      muonHits->Reset();
    }
    for (unsigned int i = 0; i < eventReader->MuonBarrelHits_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->MuonBarrelHits_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->MuonBarrelHits_cellID)[i];
      muonHits->SetNextPoint(
                            (*eventReader->MuonBarrelHits_position_x)[i] * mm,
                            (*eventReader->MuonBarrelHits_position_y)[i] * mm,
                            (*eventReader->MuonBarrelHits_position_z)[i] * mm);
    }
    for (unsigned int i = 0; i < eventReader->MuonEndcapHits_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->MuonEndcapHits_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->MuonEndcapHits_cellID)[i];
      muonHits->SetNextPoint(
                            (*eventReader->MuonEndcapHits_position_x)[i] * mm,
                            (*eventReader->MuonEndcapHits_position_y)[i] * mm,
                            (*eventReader->MuonEndcapHits_position_z)[i] * mm);
    }
  }

  //
  // digis (VTX/DCH/Wrapper)
  //
  if (displayConfig.getBoolConfig("drawVertexDigis"))
  {
    if (digis == nullptr)
    {
      digis = new TEveElementList("digis");
      gEve->AddElement(digis);
    }
    // do we need to Reset() otherwise ?
    std::cout << "Creating vertex digis" << std::endl;
    if (vtxDigis == nullptr)
    {
      vtxDigis = new TEvePointSet();
      // vtxDigis->SetName(Form("VTX digis (E>%.2f GeV)", HitEnergyThreshold));
      vtxDigis->SetName("VTX digis");
      vtxDigis->SetMarkerStyle(4);
      vtxDigis->SetMarkerSize(2);
      vtxDigis->SetMarkerColor(kYellow);
      digis->AddElement(vtxDigis);
    }
    else
      vtxDigis->Reset();
    for (unsigned int i = 0; i < eventReader->VertexBarrelDigis_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->VertexBarrelDigis_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->VertexBarrelDigis_cellID)[i];
      vtxDigis->SetNextPoint(
			     (*eventReader->VertexBarrelDigis_position_x)[i] * mm,
			     (*eventReader->VertexBarrelDigis_position_y)[i] * mm,
			     (*eventReader->VertexBarrelDigis_position_z)[i] * mm);
    }
    for (unsigned int i = 0; i < eventReader->VertexEndcapDigis_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->VertexEndcapDigis_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->VertexEndcapDigis_cellID)[i];
      vtxDigis->SetNextPoint(
			     (*eventReader->VertexEndcapDigis_position_x)[i] * mm,
			     (*eventReader->VertexEndcapDigis_position_y)[i] * mm,
			     (*eventReader->VertexEndcapDigis_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawDriftChamberDigis"))
  {
    if (digis == nullptr)
    {
      digis = new TEveElementList("digis");
      gEve->AddElement(digis);
    }
    // do we need to Reset() otherwise ?
    std::cout << "Creating drift chamber digis" << std::endl;
    if (dchDigis == nullptr)
    {
      dchDigis = new TEvePointSet();
      // dchDigis->SetName(Form("DCH digis (E>%.2f GeV)", HitEnergyThreshold));
      dchDigis->SetName("DCH digis");
      dchDigis->SetMarkerStyle(4);
      dchDigis->SetMarkerSize(2);
      dchDigis->SetMarkerColor(kYellow);
      digis->AddElement(dchDigis);
    }
    else
      dchDigis->Reset();
    for (unsigned int i = 0; i < eventReader->DriftChamberDigis_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->DriftChamberDigis_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->DriftChamberDigis_cellID)[i];
      dchDigis->SetNextPoint(
                            (*eventReader->DriftChamberDigis_position_x)[i] * mm,
                            (*eventReader->DriftChamberDigis_position_y)[i] * mm,
                            (*eventReader->DriftChamberDigis_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawSiWrapperDigis"))
  {
    if (digis == nullptr)
    {
      digis = new TEveElementList("digis");
      gEve->AddElement(digis);
    }
    std::cout << "Creating wrapper digis" << std::endl;
    if (siwrDigis == nullptr)
    {
      siwrDigis = new TEvePointSet();
      // siwrDigis->SetName(Form("VTX digis (E>%.2f GeV)", HitEnergyThreshold));
      siwrDigis->SetName("SiWr digis");
      siwrDigis->SetMarkerStyle(4);
      siwrDigis->SetMarkerSize(2);
      siwrDigis->SetMarkerColor(kYellow);
      digis->AddElement(siwrDigis);
    }
    else
      siwrDigis->Reset();
    for (unsigned int i = 0; i < eventReader->SiWrapperBarrelDigis_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->SiWrapperBarrelDigis_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->SiWrapperBarrelDigis_cellID)[i];
      siwrDigis->SetNextPoint(
                             (*eventReader->SiWrapperBarrelDigis_position_x)[i] * mm,
                             (*eventReader->SiWrapperBarrelDigis_position_y)[i] * mm,
                             (*eventReader->SiWrapperBarrelDigis_position_z)[i] * mm);
    }
    for (unsigned int i = 0; i < eventReader->SiWrapperEndcapDigis_position_x->GetSize(); i++)
    {
      // float E = (*eventReader->SiWrapperEndcapDigis_energy)[i];
      // if (E < HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->SiWrapperEndcapDigis_cellID)[i];
      siwrDigis->SetNextPoint(
                             (*eventReader->SiWrapperEndcapDigis_position_x)[i] * mm,
                             (*eventReader->SiWrapperEndcapDigis_position_y)[i] * mm,
                             (*eventReader->SiWrapperEndcapDigis_position_z)[i] * mm);
    }
  }

  //
  // cells (ECAL/HCAL/muon)
  //
  if (displayConfig.getBoolConfig("drawECalBarrelCells"))
  {
    std::cout << "Creating ecal barrel cells" << std::endl;
    if (digis == nullptr)
    {
      digis = new TEveElementList("digis");
      gEve->AddElement(digis);
    }
    if (ecalCells == nullptr)
    {
      ecalCells = new TEvePointSet();
      ecalCells->SetName(Form("ECAL cells (E>%.2f GeV)", CellEnergyThreshold));
      ecalCells->SetMarkerStyle(4);
      ecalCells->SetMarkerSize(2);
      ecalCells->SetMarkerColor(kYellow);
      digis->AddElement(ecalCells);
    }
    else
      ecalCells->Reset();
    for (unsigned int i = 0; i < eventReader->ECalBarrelCells_position_x->GetSize(); i++)
    {
      float E = (*eventReader->ECalBarrelCells_energy)[i];
      if (E < CellEnergyThreshold)
        continue;
      ecalCells->SetNextPoint((*eventReader->ECalBarrelCells_position_x)[i] * mm,
                              (*eventReader->ECalBarrelCells_position_y)[i] * mm,
                              (*eventReader->ECalBarrelCells_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawECalEndcapCells"))
  {
    std::cout << "Creating ecal endcap cells" << std::endl;
    if (digis == nullptr)
    {
      digis = new TEveElementList("digis");
      gEve->AddElement(digis);
    }
    if (ecalCells == nullptr)
    {
      ecalCells = new TEvePointSet();
      ecalCells->SetName(Form("ECAL cells (E>%.2f GeV)", CellEnergyThreshold));
      ecalCells->SetMarkerStyle(4);
      ecalCells->SetMarkerSize(2);
      ecalCells->SetMarkerColor(kYellow);
      digis->AddElement(ecalCells);
    }
    else {
      // dont do it if ecal barrel cells are drawn - this was already done before
      if (not displayConfig.getBoolConfig("drawECalBarrelCells"))
        ecalCells->Reset();
    }
    for (unsigned int i = 0; i < eventReader->ECalEndcapCells_position_x->GetSize(); i++)
    {
      float E = (*eventReader->ECalEndcapCells_energy)[i];
      if (E < CellEnergyThreshold)
        continue;
      ecalCells->SetNextPoint((*eventReader->ECalEndcapCells_position_x)[i] * mm,
                              (*eventReader->ECalEndcapCells_position_y)[i] * mm,
                              (*eventReader->ECalEndcapCells_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawHCalBarrelCells"))
  {
    std::cout << "Creating hcal barrel cells" << std::endl;
    if (digis == nullptr)
    {
      digis = new TEveElementList("digis");
      gEve->AddElement(digis);
    }
    if (hcalCells == nullptr)
    {
      hcalCells = new TEvePointSet();
      hcalCells->SetName(Form("HCAL cells (E>%.2f GeV)", CellEnergyThreshold));
      hcalCells->SetMarkerStyle(4);
      hcalCells->SetMarkerSize(2);
      hcalCells->SetMarkerColor(kYellow);
      digis->AddElement(hcalCells);
    }
    else
      hcalCells->Reset();
    for (unsigned int i = 0; i < eventReader->HCalBarrelCells_position_x->GetSize(); i++)
    {
      float E = (*eventReader->HCalBarrelCells_energy)[i];
      if (E < CellEnergyThreshold)
        continue;
      hcalCells->SetNextPoint((*eventReader->HCalBarrelCells_position_x)[i] * mm,
                              (*eventReader->HCalBarrelCells_position_y)[i] * mm,
                              (*eventReader->HCalBarrelCells_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawHCalEndcapCells"))
  {
    std::cout << "Creating hcal endcap cells" << std::endl;
    if (digis == nullptr)
    {
      digis = new TEveElementList("digis");
      gEve->AddElement(digis);
    }
    if (hcalCells == nullptr)
    {
      hcalCells = new TEvePointSet();
      hcalCells->SetName(Form("HCAL cells (E>%.2f GeV)", CellEnergyThreshold));
      hcalCells->SetMarkerStyle(4);
      hcalCells->SetMarkerSize(2);
      hcalCells->SetMarkerColor(kYellow);
      digis->AddElement(hcalCells);
    }
    else {
      // dont do it if hcal barrel cells are drawn - this was already done before
      if (not displayConfig.getBoolConfig("drawHCalBarrelCells"))
        hcalCells->Reset();
    }
    for (unsigned int i = 0; i < eventReader->HCalEndcapCells_position_x->GetSize(); i++)
    {
      float E = (*eventReader->HCalEndcapCells_energy)[i];
      if (E < CellEnergyThreshold)
        continue;
      hcalCells->SetNextPoint((*eventReader->HCalEndcapCells_position_x)[i] * mm,
                              (*eventReader->HCalEndcapCells_position_y)[i] * mm,
                              (*eventReader->HCalEndcapCells_position_z)[i] * mm);
    }
  }

  if (displayConfig.getBoolConfig("drawMuonCells"))
  {
    std::cout << "Creating muon barrel cells" << std::endl;
    if (digis == nullptr)
    {
      digis = new TEveElementList("digis");
      gEve->AddElement(digis);
    }
    if (muonCells == nullptr)
    {
      muonCells = new TEvePointSet();
      muonCells->SetName(Form("MUON cells (E>%.2f GeV)", CellEnergyThreshold));
      muonCells->SetMarkerStyle(4);
      muonCells->SetMarkerSize(2);
      muonCells->SetMarkerColor(kYellow);
      digis->AddElement(muonCells);
    }
    else
      muonCells->Reset();
    for (unsigned int i = 0; i < eventReader->MuonBarrelCells_position_x->GetSize(); i++)
    {
      float E = (*eventReader->MuonBarrelCells_energy)[i];
      if (E < CellEnergyThreshold)
        continue;
      muonCells->SetNextPoint((*eventReader->MuonBarrelCells_position_x)[i] * mm,
                              (*eventReader->MuonBarrelCells_position_y)[i] * mm,
                              (*eventReader->MuonBarrelCells_position_z)[i] * mm);
    }
    std::cout << "Creating muon endcap cells" << std::endl;
    for (unsigned int i = 0; i < eventReader->MuonEndcapCells_position_x->GetSize(); i++)
    {
      float E = (*eventReader->MuonEndcapCells_energy)[i];
      if (E < CellEnergyThreshold)
        continue;
      muonCells->SetNextPoint((*eventReader->MuonEndcapCells_position_x)[i] * mm,
                              (*eventReader->MuonEndcapCells_position_y)[i] * mm,
                              (*eventReader->MuonEndcapCells_position_z)[i] * mm);
    }
  }

  //
  // cells merged (ECAL)
  //
  if (displayConfig.getBoolConfig("drawECalBarrelMergedCells"))
  {
    std::cout << "Creating merged ecal barrel cells" << std::endl;
    if (cells_merged == nullptr)
    {
      cells_merged = new TEvePointSet();
      cells_merged->SetName("cells_merged");
      cells_merged->SetMarkerStyle(4);
      cells_merged->SetMarkerSize(3);
      cells_merged->SetMarkerColor(kBlue);
      digis->AddElement(cells_merged);
    }
    else
      cells_merged->Reset();
    for (unsigned int i = 0; i < eventReader->ECalBarrelCells2_position_x->GetSize(); i++)
    {
      float E = (*eventReader->ECalBarrelCells2_energy)[i];
      // if (E<minCellE) continue;
      // ULong_t cellID = (*eventReader->ECalBarrelCells_cellID)[i];
      // ULong_t layer = DetectorGeometry::Layer(cellID);
      cells_merged->SetNextPoint((*eventReader->ECalBarrelCells2_position_x)[i] * mm,
                                 (*eventReader->ECalBarrelCells2_position_y)[i] * mm,
                                 (*eventReader->ECalBarrelCells2_position_z)[i] * mm);
    }
  }

  //
  // tracks
  //
  if (displayConfig.getBoolConfig("drawTracks"))
  {
    std::cout << "Creating tracks" << std::endl;
    if (tracks == nullptr)
    {
      tracks = new TEveTrackList("tracks");
      TEveTrackPropagator *trkProp = tracks->GetPropagator();
      // trkProp->SetMagField(-2.0); // tesla
      trkProp->SetMagFieldObj(magField);
      trkProp->SetMaxR(rMax);
      // trkProp->SetMaxZ(geomReader->zMax);
      trkProp->SetMaxZ(geomReader->zMaxEndCap);
      trkProp->SetRnrReferences(true);
      trkProp->SetFitReferences(false);
      tracks->SetMainColor(kCyan);
      tracks->SetLineWidth(2.);
      tracks->SetLineStyle(5);
      gEve->AddElement(tracks);
    }
    else
      tracks->DestroyElements();

    for (unsigned int ip = 0; ip < eventReader->TracksFromGenParticles_subdetectorHitNumbers_begin->GetSize(); ip++)
    {
      if (debug) std::cout << "Track " << ip << std::endl;
      unsigned int trackStates_begin = (*eventReader->TracksFromGenParticles_trackStates_begin)[ip];
      unsigned int trackStates_end = (*eventReader->TracksFromGenParticles_trackStates_end)[ip];
      if (debug) std::cout << "  trackStates begin , end = " << trackStates_begin << " , " << trackStates_end << std::endl;
      // store in vectors the track state positions and momenta
      unsigned int nTrackStates = 1 + trackStates_end - trackStates_begin;
      // std::vector<float> x(nTrackStates), y(nTrackStates), z(nTrackStates), omega(nTrackStates), tanLambda(nTrackStates), phi(nTrackStates);
      // for (unsigned int i=0; i<nTrackStates; i++) x[i]=-1e6; // set to some large value to check later if track state has been filled
      float x[5], y[5], z[5], omega[5], tanLambda[5], phi[5];
      for (unsigned int i=0; i<5; i++) x[i]=-1e6; // set to some large value to check later if track state has been filled
      for (unsigned int its = trackStates_begin; its < trackStates_end; its++) {
	int location = (*eventReader->_TracksFromGenParticles_trackStates_location)[its];
	if (debug) {
	  std::cout << "Trackstate " << its << std::endl;
	  std::cout << "  location = " << location << std::endl;
	}
	// if (location<1 or location>4) continue;
	if (location>4) continue;
	x[location] = (*eventReader->_TracksFromGenParticles_trackStates_referencePoint_x)[its];
	y[location] = (*eventReader->_TracksFromGenParticles_trackStates_referencePoint_y)[its];
	z[location] = (*eventReader->_TracksFromGenParticles_trackStates_referencePoint_z)[its];
	omega[location] = (*eventReader->_TracksFromGenParticles_trackStates_omega)[its];
	tanLambda[location] = (*eventReader->_TracksFromGenParticles_trackStates_tanLambda)[its];
	phi[location] = (*eventReader->_TracksFromGenParticles_trackStates_phi)[its];
	if (debug) {
	  std::cout << "  x = " << x[location-1] << std::endl;
	  std::cout << "  y = " << y[location-1] << std::endl;
	  std::cout << "  z = " << z[location-1] << std::endl;
	  std::cout << "  omega = " << omega[location-1] << std::endl;
	  std::cout << "  phi = " << phi[location-1] << std::endl;
	  std::cout << "  tanLambda = " << tanLambda[location-1] << std::endl;
	}
      }
      const int nSubDetectors(5);
      int nhits[nSubDetectors];
      unsigned int hits_begin = (*eventReader->TracksFromGenParticles_subdetectorHitNumbers_begin)[ip];
      unsigned int hits_end = (*eventReader->TracksFromGenParticles_subdetectorHitNumbers_end)[ip];
      if (debug) std::cout << "  subdetectorHitNumbers_begin , end = " << hits_begin << " , " << hits_end << std::endl;
      for (unsigned int ih = hits_begin; ih < hits_end; ih++) {
	nhits[ih-hits_begin] = (*eventReader->_TracksFromGenParticles_subdetectorHitNumbers)[ih];
      }
      
      const int tsOrig=2; // which track state to use for track origin: 0=at other, 1=at IP, 2=at first hit, 3=at last hit, 4=at ECAL
      if (std::fabs(x[tsOrig]>-5e5)) {
	TEveRecTrack t;
	float x1 = x[tsOrig]*mm;
	float y1 = y[tsOrig]*mm;
	float z1 = z[tsOrig]*mm;
	const float FCT = 2.99792458E-4f;
	const float bField = 2.0;
	int charge = sgn(omega[tsOrig]);
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
	for (int i=0; i<5; i++) {
	  if (std::fabs(x[i]<-5e5)) continue;
	  float radius = 1.f / std::fabs(omega[i]);
	  float pxy = FCT * bField * radius;
	  float px = pxy * std::cos(phi[i]);
	  float py = pxy * std::sin(phi[i]);
	  float pz = tanLambda[i] * pxy;
	  track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kReference,
					   TEveVectorD(x[i]*mm, y[i]*mm, z[i]*mm),
					   TEveVectorD(px, py, py)));
	  track->SetRnrPoints(true);
	  track->SetMarkerStyle(4);
	}
	// could also save decay ...	
	tracks->AddElement(track);
      }
    }
    tracks->MakeTracks();
  }

  
  //
  // clusters
  //
  if (displayConfig.getBoolConfig("drawTopoClusters"))
  {
    FillClusters("topo");
    DrawClusters("topo");
  }
  if (displayConfig.getBoolConfig("drawCaloClusters"))
  {
    FillClusters("sw");
    DrawClusters("sw");
  }

  //
  // event label
  //
  if (eventLabel == nullptr)
  {
    eventLabel = new TGLConstAnnotation(mainGLView,
                                        Form("%s, %.1f GeV\nEvent %d",
                                             partType(pdgID), pmax, eventId),
                                        0.1, 0.9);
    eventLabel->SetTextSize(0.05); // % of window diagonal
    eventLabel->SetAllowClose(false);
  }
  else
  {
    eventLabel->SetText(Form("%s, %.1f GeV\nEvent %d", partType(pdgID), pmax, eventId));
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

  textEntry->SetTextColor((Pixel_t)0x000000);
  textEntry->SetText(Form("Event %d loaded", eventId));
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
  if (geomFile.find("v03") != std::string::npos)
    detectorVersion = 3;
  else
    detectorVersion = 2;
  geomReader = new DetectorGeometry(detectorVersion);

  // create magnetic field
  magField = new MagField(geomReader->Bin, geomReader->Bout, geomReader->zMax, geomReader->rMax);

  std::cout << "******************************************************************************" << std::endl;
  std::cout << "Displaying the geometry" << std::endl;
  std::cout << "******************************************************************************" << std::endl
            << endl;

  // create the eve manageer
  TEveManager::Create();

  // Set title of main window
  gEve->GetBrowser()->SetWindowName("FCC-ee ALLEGRO detector event display");

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
  gEve->GetBrowser()->GetTabRight()->Connect("Selected(Int_t)", "EventDisplay", this, "onTabSelected(Int_t)");

  // Create the geometry and the readouts
  TEveElementList *geom = new TEveElementList("Geometry");
  TEveElementList *PCBs = new TEveElementList("PCBs");
  TEveElementList *actives = new TEveElementList("Active elements");
  TEveElementList *passives = new TEveElementList("Passive elements");
  TEveElementList *readout = new TEveElementList("Readout");
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

    if (showFullDetector)
    {
      // geom->AddElement(world);
      // world->SetRnrSelfChildren(false, true);

      // group together elements
      TPRegexp re_lc("LumiCal_*");
      TPRegexp re_dch("DCH_*");
      TPRegexp re_siwrapb("SiWrB_*");
      TPRegexp re_siwrapec("SiWrD_*");
      TPRegexp re_ecalb("ECalBarrel*");
      TPRegexp re_hcalb("HCalEnvelopeVolume*");
      TPRegexp re_hcalec("HCal(\\w+)Endcap(\\w+)");
      TPRegexp re_muonb("MuonTaggerBarrel*");
      TPRegexp re_muonec("MuonTaggerEndcap*");

      TEveElementList *mdi = new TEveElementList("MDI");
      geom->AddElement(mdi);
      // disable MDI rendering by default
      mdi->SetRnrSelfChildren(false, false);

      TEveElementList *beampipe = new TEveElementList("Beampipe");
      geom->AddElement(beampipe);
      // disable beampipe rendering by default
      beampipe->SetRnrSelfChildren(false, false);

      TEveElementList *lumical = new TEveElementList("LumiCal");
      geom->AddElement(lumical);
      // disable lumical rendering by default
      lumical->SetRnrSelfChildren(false, false);

      TEveElementList *vertexBarrel = new TEveElementList("Vertex barrel");
      geom->AddElement(vertexBarrel);

      TEveElementList *vertexEndcap = new TEveElementList("Vertex endcaps");
      geom->AddElement(vertexEndcap);

      TEveElementList *dch = new TEveElementList("Drift chamber");
      geom->AddElement(dch);

      TEveElementList *siwrapb = new TEveElementList("Si wrapper barrel");
      geom->AddElement(siwrapb);

      TEveElementList *siwrapec = new TEveElementList("Si wrapper endcap");
      geom->AddElement(siwrapec);

      TEveElementList *ecalb = new TEveElementList("ECal barrel");
      geom->AddElement(ecalb);

      TEveElementList *ecalec = new TEveElementList("ECal endcaps");
      geom->AddElement(ecalec);

      TEveElementList *hcalb = new TEveElementList("HCal barrel");
      geom->AddElement(hcalb);

      TEveElementList *hcalec = new TEveElementList("HCal endcaps");
      geom->AddElement(hcalec);

      TEveElementList *muontaggerb = new TEveElementList("Muon tagger barrel");
      geom->AddElement(muontaggerb);

      TEveElementList *muontaggerec = new TEveElementList("Muon tagger endcaps");
      geom->AddElement(muontaggerec);

      for (TEveElement::List_i itr = world->BeginChildren(); itr != world->EndChildren(); itr++)
      {
        TEveElement *a = *itr;
        TString s(a->GetElementName());

        if (s.Contains("ScreenSol") || s.Contains("CompSol")) {
          cout << "Adding " << s << " to MDI" << endl;
          mdi->AddElement(a);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
        }
        else if (s.Contains("BeamPipe") || s.Contains("Beampipe")) {
          cout << "Adding " << s << " to beampipe" << endl;
          beampipe->AddElement(a);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
        }
        else if (re_lc.MatchB(s)) {
          cout << "Adding " << s << " to lumical" << endl;
          lumical->AddElement(a);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
        }
        else if (s.Contains("Vertex"))
        {
          cout << "Found " << s << ", looping over its children" << endl;
          // for vertex, loop over children to add them either to barrel or endcap
          for (TEveElement::List_i itrVtx = a->BeginChildren(); itrVtx != a->EndChildren(); itrVtx++)
          {
            TEveElement *elVtx = *itrVtx;
            TString sVtx(elVtx->GetElementName());
            if (sVtx.Contains("VertexBarrel") or sVtx.Contains("VertexInnerBarrel")) {
              if (debug) cout << "Adding " << sVtx << " to vertex barrel" << endl;
              // vertexBarrel->AddElement(elVtx);
              // elVtx->SetMainColor(kRed);
              // rather than the assembly we add directly the layers
              for (TEveElement::List_i itr2 = elVtx->BeginChildren(); itr2 != elVtx->EndChildren(); itr2++) {
                TEveElement *el = *itr2;
                TString elName(el->GetElementName());
                vertexBarrel->AddElement(el);
                el->SetMainColor(kRed);
                ((TEveGeoShape *)el)->SetDrawFrame(false);
                if (elName.Contains("VertexInnerBarrel"))
                  ((TEveGeoShape *)el)->SetNSegments(128);
              }
            }
            else if (sVtx.Contains("VertexDisks")) {
              if (debug) cout << "Adding " << sVtx << " to vertex endcap" << endl;
              vertexEndcap->AddElement(elVtx);
              elVtx->SetMainColor(kRed);
              // rather than the assembly we add directly the layers
              //for (TEveElement::List_i itr2 = elVtx->BeginChildren(); itr2 != elVtx->EndChildren(); itr2++) {
              //TEveElement *el = *itr2;
              //vertexEndcap->AddElement(el);
              //el->SetMainColor(kRed);
              //((TEveGeoShape *)el)->SetDrawFrame(false);
              //}
            }
            else
              std::cout << "Unexpected volume: " << sVtx << std::endl;
          }
        }
        else if (re_dch.MatchB(s))
        {
          cout << "Adding " << s << " to drift chamber" << endl;
          dch->AddElement(a);
          a->SetRnrSelfChildren(true, false);
          a->SetMainTransparency(useTransparencies ? 60 : 0);
          //a->SetMainTransparency(0);
          a->SetMainColor(kGreen-5);
          ((TEveGeoShape *)a)->SetNSegments(128);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
        }
        else if (re_siwrapb.MatchB(s))
        {
          cout << "Adding " << s << " to Si-wrapper barrel" << endl;
          // rather than the assembly we add directly the layers
          for (TEveElement::List_i itr2 = a->BeginChildren(); itr2 != a->EndChildren(); itr2++) {
            TEveElement *el = *itr2;
            siwrapb->AddElement(el);
            el->SetRnrSelfChildren(true, false);
            el->SetMainTransparency(useTransparencies ? 60 : 0);
            el->SetMainColor(kRed);
            ((TEveGeoShape *)el)->SetNSegments(128);
            ((TEveGeoShape *)el)->SetDrawFrame(false);
          }
        }
        else if (re_siwrapec.MatchB(s))
        {
          cout << "Adding " << s << " to Si-wrapper endcap" << endl;
          // rather than the appendix we add directly the layers
          for (TEveElement::List_i itr2 = a->BeginChildren(); itr2 != a->EndChildren(); itr2++) {
            TEveElement *el = *itr2;
            siwrapec->AddElement(el);
            el->SetRnrSelfChildren(true, false);
            el->SetMainTransparency(useTransparencies ? 60 : 0);
            el->SetMainColor(kRed);
            ((TEveGeoShape *)el)->SetNSegments(128);
            ((TEveGeoShape *)el)->SetDrawFrame(false);
          }
        }
        else if (re_ecalb.MatchB(s))
        {
          cout << "Adding " << s << " to ECal barrel" << endl;
          // add the overall envelope
          // we draw the envelope in the 3d model,
          // and the volumes in the 2D views
          ecalb->AddElement(a);
          //a->SetRnrSelfChildren(false, true);
          a->SetRnrSelfChildren(true, false);
          a->SetMainColor(kAzure+7);
          a->SetMainTransparency(useTransparencies ? 60 : 0);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
          for (TEveElement::List_i itr2 = a->BeginChildren(); itr2 != a->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            el->SetMainColor(kAzure+7);
            TString elname(el->GetElementName());
            std::cout << elname << std::endl;
            el->SetRnrSelfChildren(false, false);
            if (elname.BeginsWith("ECAL_Cryo_side"))
              el->SetRnrSelfChildren(false, false);
            else
              el->SetRnrSelfChildren(true, false);
            el->SetMainTransparency(useTransparencies ? 60 : 0);
            ((TEveGeoShape *)el)->SetNSegments(128);
            ((TEveGeoShape *)el)->SetDrawFrame(false);
          }
        }
        else if (s.Contains("ECalEndcaps"))
        {
          cout << "Adding " << s << " to ECal endcap" << endl;
          // add the two envelopes in the two sides
          ecalec->AddElement(a);
          a->SetRnrSelfChildren(true, false);
          a->SetMainColor(kAzure+7);
          a->SetMainTransparency(useTransparencies ? 60 : 0);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
        }
        else if (re_hcalb.MatchB(s))
        {
          cout << "Adding " << s << " to HCal barrel" << endl;
          hcalb->AddElement(a);
          // for the HCal barrel, we draw the envelope in the 3d model,
          // and the subvolumes in the 2D views
          a->SetRnrSelfChildren(true, false);
          a->SetMainTransparency(useTransparencies ? 60 : 0);
          a->SetMainColor(kAzure-7);
          ((TEveGeoShape *)a)->SetNSegments(128);
          ((TEveGeoShape *)a)->SetDrawFrame(false);

          //a->SetRnrSelfChildren(false, true);
          for (TEveElement::List_i itr2 = a->BeginChildren(); itr2 != a->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            //el->SetRnrSelfChildren(true, false);
            el->SetMainColor(kAzure-7);
            el->SetMainTransparency(useTransparencies ? 60 : 0);
            ((TEveGeoShape *)el)->SetNSegments(128);
            ((TEveGeoShape *)el)->SetDrawFrame(false);
          }
        }
        else if (re_hcalec.MatchB(s))
        {
          cout << "Adding " << s << " to HCal endcap" << endl;
          hcalec->AddElement(a);
          // for the HCal endcap, we draw the envelope in the 3d model,
          // and the volumes in the 2D views
          a->SetRnrSelfChildren(true, false);
          //a->SetRnrSelfChildren(false, true);
          a->SetMainColor(kAzure-7);
          a->SetMainTransparency(useTransparencies ? 60 : 0);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
          for (TEveElement::List_i itr2 = a->BeginChildren(); itr2 != a->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            // el->SetRnrSelfChildren(true, false);
            el->SetMainColor(kAzure-7);
            el->SetMainTransparency(useTransparencies ? 60 : 0);
            ((TEveGeoShape *)el)->SetDrawFrame(false);
          }
        }
        else if (re_muonb.MatchB(s))
        {
          cout << "Adding " << s << " to Muon tagger barrel" << endl;
          muontaggerb->AddElement(a);
          a->SetRnrSelfChildren(true, false);
          a->SetMainTransparency(useTransparencies ? 60 : 0);
          a->SetMainColor(kOrange);
          ((TEveGeoShape *)a)->SetNSegments(128);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
        }
        else if (re_muonec.MatchB(s))
        {
          cout << "Adding " << s << " to Muon tagger endcap" << endl;
          muontaggerec->AddElement(a);
          a->SetRnrSelfChildren(true, false);
          a->SetMainColor(kOrange);
          a->SetMainTransparency(useTransparencies ? 60 : 0);
          ((TEveGeoShape *)a)->SetDrawFrame(false);
        }
        else
        {
          std::cout << "Unexpected volume: " << s << std::endl;
          geom->AddElement(a);
        }
      }
    }
    else
    {

      //
      // ECAL barrel
      //
      TPRegexp re("ECalBarrel*");
      TEveElement *ecalbarrel = world->FindChild(re);
      ecalbarrel->SetPickableRecursively(kTRUE);
      ((TEveGeoShape *)ecalbarrel)->SetNSegments(256);
      geom->AddElement(ecalbarrel);
      // hide the envelope, only show physical volumes in it
      ecalbarrel->SetRnrSelfChildren(false, true);

      // set transparency of the subvolumes of the bath
      re = TPRegexp("LAr_bath*");
      TEveElement *bath = ecalbarrel->FindChild(re);
      TEveElement::List_t matches;
      re = TPRegexp("ECAL_Cryo*");
      ecalbarrel->FindChildren(matches, re);
      for (TEveElement *a : matches)
      {
        a->SetMainTransparency(70);
        ((TEveGeoShape *)a)->SetNSegments(256);
      }
      re = TPRegexp("services*");
      ecalbarrel->FindChildren(matches, re);
      for (auto a : matches)
      {
        a->SetMainTransparency(70);
        ((TEveGeoShape *)a)->SetNSegments(256);
      }

      // turn off Cryo side in 3d, weird rendering
      re = TPRegexp("ECAL_Cryo_side_(\\w+)");
      TEveElement::List_t matches2;
      ecalbarrel->FindChildren(matches2, re);
      for (auto a : matches2)
      {
        a->SetRnrSelf(false);
      }
      // make lists of elements inside bath to turn on/off simultaneously
      if (bath)
      {
        TEveElementList *newbath = new TEveElementList("LAr_bath");
        ecalbarrel->AddElement(newbath);
        TEveElement::List_t matches;
        re = TPRegexp("PCB(\\w+)");
        bath->FindChildren(matches, re);
        for (auto a : matches)
          PCBs->AddElement(a);
        newbath->AddElement(PCBs);
        TEveElement::List_t matches2;
        re = TPRegexp("active*");
        bath->FindChildren(matches2, re);
        for (TEveElement *a : matches2)
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
        ecalbarrel->RemoveElement(bath);
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
      if (doHCal)
      {
        // add HCal envelope and subvolumes
        re = TPRegexp("HCalEnvelopeVolume*");
        TEveElement *hcalbarrel = world->FindChild(re);
        hcalbarrel->SetPickableRecursively(kTRUE);
        ((TEveGeoShape *)hcalbarrel)->SetNSegments(256);
        geom->AddElement(hcalbarrel);
        // do not draw envelope, only physical elements inside
        hcalbarrel->SetRnrSelfChildren(false, true);

        // re = TPRegexp("HCalLayerVol*");
        // hcalbarrel->FindChildren(matches, re);
        // for (auto a : matches) a->SetRnrSelfChildren(true, false);

        // set transparency and number of segments in 2D projection of plate faces and steel
        {
          re = TPRegexp("HCal(\\w+)PlateVol(\\w+)");
          TEveElement::List_t matches;
          hcalbarrel->FindChildren(matches, re);
          for (TEveElement *a : matches)
          {
            a->SetMainTransparency(70);
            ((TEveGeoShape *)a)->SetNSegments(256);
          }
        }

        re = TPRegexp("HCal*Steel*");
        hcalbarrel->FindChildren(matches, re);
        for (auto a : matches)
        {
          a->SetMainTransparency(70);
          ((TEveGeoShape *)a)->SetNSegments(256);
        }

        // group together the layers so that they can be turned on/off together
        // set number of segments of 2D projection of layer envelopes
        TEveElementList *hcalLayers = new TEveElementList("HCalLayers");
        hcalbarrel->AddElement(hcalLayers);
        re = TPRegexp("HCalLayerVol*");
        TEveElement::List_t matches4;
        hcalbarrel->FindChildren(matches4, re);
        for (TEveElement *a : matches4)
        {
          hcalLayers->AddElement(a);
          hcalbarrel->RemoveElement(a);
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
        TPRegexp re("ECAL_Cryo_side_(\\w+)");
        TEveElement *projcryoside = projbarrel->FindChild(re);
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
          for (TEveElement::List_i itr = a->BeginChildren(); itr != a->EndChildren(); itr++) {
            TEveElement *b = *itr;
            TString s(b->GetElementName());
            if (s.Contains("assembly")) {
              b->SetRnrSelfChildren(false, true);
            }
          }
        }

        // show active volumes in ecal barrel
        else if (s.BeginsWith("ECal barrel"))
        {
          TPRegexp re("ECalBarrel_vol_(\\w+)");
          TEveElement *projbarrel = a->FindChild(re);
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
          TPRegexp re("HCalEnvelopeVolume*");
          TEveElement *envelope = a->FindChild(re);
          envelope->SetRnrSelfChildren(false, true);
          // envelope->SetMainColor(kAzure-7);
          ((TEveGeoShape *)envelope)->SetDrawFrame(false);
          for (TEveElement::List_i itr2 = envelope->BeginChildren(); itr2 != envelope->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            el->SetRnrSelfChildren(true, false);
            // el->SetMainColor(kAzure-7);
            // el->SetMainTransparency(useTransparencies ? 60 : 0);
            // TEveGeoShape* gs = (TEveGeoShape *)el;
            // if (gs) {
            //   gs->SetDrawFrame(false);
            // }
          }
        }
      }
    }
  }

  // draw the merged ECAL readout segmentation in rho-phi
  TEveStraightLineSet *gridmod = new TEveStraightLineSet("ECAL phi readout merged");
  gridmod->SetLineColor(kViolet + 2);
  gridmod->SetLineWidth(5);
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
    hcalPhiReadout->SetLineColor(kViolet + 2);
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
  grid2->SetLineColor(kViolet + 2);
  grid2->SetLineWidth(5);
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
          for (TEveElement::List_i itr = a->BeginChildren(); itr != a->EndChildren(); itr++) {
            TEveElement *b = *itr;
            TString s(b->GetElementName());
            if (s.Contains("assembly")) {
              b->SetRnrSelfChildren(false, true);
            }
          }
        }
        // show vessel and gas for DCH
        else if (s.BeginsWith("Drift chamber")) {
          a->SetRnrSelfChildren(false, true);
          for (TEveElement::List_i itr = a->BeginChildren(); itr != a->EndChildren(); itr++) {
            TEveElement *b = *itr;
            b->SetRnrSelfChildren(true, false);
            b->SetMainColor(kGreen-5);
            ((TEveGeoShape *)b)->SetDrawFrame(false);
          }
        }
        // turn off the Si wrapper endcap envelope since it does not render well
        // instead show its daughter volume instead
        else if (s.BeginsWith("Si wrapper endcap")) {
          a->SetRnrSelfChildren(false, true);
          for (TEveElement::List_i itr = a->BeginChildren(); itr != a->EndChildren(); itr++) {
            TEveElement *b = *itr;
            TString s(b->GetElementName());
            if (s.BeginsWith("SiWrD_envelope")) {
              b->SetRnrSelfChildren(false, true);
              for (TEveElement::List_i itr = b->BeginChildren(); itr != b->EndChildren(); itr++) {
                TEveElement *c = *itr;
                c->SetRnrSelfChildren(true, false);
              }
            }
          }
        }
        // ECAL barrel
        else if (s.BeginsWith("ECal barrel")) {
          TPRegexp re("ECalBarrel_vol_(\\w+)");
          TEveElement *projbarrel = a->FindChild(re);
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
        // HCAL endcap - do not draw the envelope
        else if (s.BeginsWith("HCal endcap"))
        {
          TPRegexp re("HCalThreePartsEndcap*");
          TEveElement *envelope = a->FindChild(re);
          envelope->SetRnrSelfChildren(false, true);
          // envelope->SetMainColor(kAzure-7);
          ((TEveGeoShape *)envelope)->SetDrawFrame(false);
          for (TEveElement::List_i itr2 = envelope->BeginChildren(); itr2 != envelope->EndChildren(); itr2++)
          {
            TEveElement *el = *itr2;
            TString elname(el->GetElementName());
            el->SetRnrSelfChildren(true, false);
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

  // gEve->DoRedraw3D();
  // gEve->FullRedraw3D(true);


  // create the gui for event navigation
  makeGui();

  //
  // display the events
  //
  if (debug)
    std::cout << "evtFile: " << evtFile << std::endl;
  if (evtFile != "")
  {
    // create the gui for event navigation
    //makeGui();


    // setup the event reader
    std::cout << std::endl;
    std::cout << "******************************************************************************" << std::endl;
    std::cout << "Setting up the event reader" << std::endl;
    std::cout << "******************************************************************************" << std::endl
              << std::endl;

    std::cout << "Reading event data from file " << evtFile << std::endl
              << std::endl;

    // open the file with the events and create the reader
    TFile *f = TFile::Open(evtFile.c_str(), "READ");
    eventReader = new EventReader(f, doHCal);

    // setup the branches
    eventReader->SetBranches();

    // print updated draw settings based on info found in ROOT file
    displayConfig.Print();

    // read the number of events in the file
    nEvents = eventReader->nEvents;

    // load and display the requested event
    cout << endl;
    cout << "******************************************************************************" << endl;
    cout << "Reading the events" << endl;
    cout << "******************************************************************************" << endl
         << endl;
    loadEvent(initialEvent);
  }

  // Set the 3D view as the active tab and rotate the camera
  gEve->GetBrowser()->GetTabRight()->SetTab(0);

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
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText("No events loaded");
    printf("\nNo events loaded\n");
  }
  else if (eventId < nEvents - 1)
  {
    ++eventId;
    loadEvent(eventId);
  }
  else
  {
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText("Already at last event");
    printf("\nAlready at last event\n");
  }
}

void EventDisplay::bck()
{
  if (evtFile == "")
  {
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText("No events loaded");
    printf("\nNo events loaded\n");
  }
  else if (eventId > 0)
  {
    --eventId;
    loadEvent(eventId);
  }
  else
  {
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText("Already at first event");
    printf("\nAlready at first event\n");
  }
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
    cout << "Unknown view!!" << endl;

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

void EventDisplay::makeGui()
{
  // Create minimal GUI for event navigation.

  TEveBrowser *browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);

  TGMainFrame *frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("GUI");
  frmMain->SetCleanup(kDeepCleanup);

  TGHorizontalFrame *hf = new TGHorizontalFrame(frmMain);
  {
    TGPictureButton *b = 0;

    b = new TGPictureButton(hf, gClient->GetPicture("icons/TakeScreenshot.png"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "takeScreenshot()");
    b->SetToolTipText("Take a screenshot of the selected view");

    b = new TGPictureButton(hf, gClient->GetPicture("icons/GoBack.gif"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "bck()");
    b->SetToolTipText("Go to previous event");

    b = new TGPictureButton(hf, gClient->GetPicture("icons/GoForward.gif"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 10, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "fwd()");
    b->SetToolTipText("Go to next event");

    textEntry = new TGTextEntry(hf);
    textEntry->SetEnabled(kFALSE);
    hf->AddFrame(textEntry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY |
                                                  kLHintsExpandX,
                                              2, 10, 10, 10));
  }
  frmMain->AddFrame(hf);

  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();

  browser->StopEmbedding();
  browser->SetTabTitle("Controls", 0);
}

void EventDisplay::onTabSelected(Int_t tab)
{
  // std::cout << "tab selected: " << tab << std::endl;
  if (tab==0) activeGLViewer = mainGLView;
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
