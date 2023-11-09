/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventDisplay: creates and populates the event display
//
/******************************************************************************/

#include "EventDisplay.h"
#include "DetectorGeometry.h"

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
#include <TEveStraightLineSet.h>
#include <TEveTrackPropagator.h>
#include <TEveProjectionAxes.h>

#include <filesystem>
#include <unordered_map>

// return the sign of a float
int sgn(float val)
{
  return (val > 0.) - (val < 0.);
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

  if (clusterType == "topo")
    std::cout << "Creating topo clusters" << std::endl;
  else if (clusterType == "sw")
    std::cout << "Creating SW clusters" << std::endl;
  else
  {
    std::cout << "Unknown cluster type " << clusterType << std::endl;
    return;
  }

  TEvePointSet *clusters = nullptr;
  TEveElementList *clusters_3D = nullptr;
  TEveElementList *clusters_rhoz = nullptr;
  TEveElementList *clusters_rhophi = nullptr;

  const double cm = geomReader->cm;
  const double mm = geomReader->mm;
  const double rMin = geomReader->rMin;
  double rMax = doHCal ? geomReader->rMaxHCal : geomReader->rMax;
  const double alpha = geomReader->alpha;
  const double thetaGrid = geomReader->thetaGrid;
  const double gridPhi = geomReader->gridPhi;
  const double phiMin = geomReader->phiMin;

  if (clusterType == "topo")
  {
    clusters = topoclusters;
  }
  else
  {
    clusters = swclusters;
  }

  // centers of the clusters
  if (clusterType == "topo")
  {
    if (topoclusters == nullptr)
    {
      topoclusters = new TEvePointSet();
      topoclusters->SetName(Form("%s clusters (E>%.1f GeV)", clusterType.c_str(), ClusterEnergyThreshold));
      topoclusters->SetMarkerColor(kGreen);
      gEve->AddElement(topoclusters);
    }
    else
      topoclusters->Reset();
    clusters = topoclusters;
  }
  else
  {
    if (swclusters == nullptr)
    {
      swclusters = new TEvePointSet();
      swclusters->SetName(Form("%s clusters (E>%.1f GeV)", clusterType.c_str(), ClusterEnergyThreshold));
      swclusters->SetMarkerColor(kGray);
      gEve->AddElement(swclusters);
    }
    else
      swclusters->Reset();
    clusters = swclusters;
  }
  clusters->SetMarkerStyle(4);
  clusters->SetMarkerSize(6);

  unsigned int nClusters = (clusterType == "topo") ? eventReader->CorrectedCaloTopoClusters_position_x->GetSize() : eventReader->CaloClusters_position_x->GetSize();

  for (unsigned int i = 0; i < nClusters; i++)
  {
    float E = (clusterType == "topo") ? (*eventReader->CorrectedCaloTopoClusters_energy)[i] : (*eventReader->CaloClusters_energy)[i];
    if (E < ClusterEnergyThreshold)
      continue;
    // topo cluster positions are in cm while calo clusters and hits/cells in mm ... ?
    if (clusterType == "topo")
    {
      clusters->SetNextPoint((*eventReader->CorrectedCaloTopoClusters_position_x)[i] * cm,
                             (*eventReader->CorrectedCaloTopoClusters_position_y)[i] * cm,
                             (*eventReader->CorrectedCaloTopoClusters_position_z)[i] * cm);
    }
    else
    {
      clusters->SetNextPoint((*eventReader->CaloClusters_position_x)[i] * mm,
                             (*eventReader->CaloClusters_position_y)[i] * mm,
                             (*eventReader->CaloClusters_position_z)[i] * mm);
    }
  }

  // now create the visual representation of the clustered cells
  std::vector<TEveQuadSet *> qs_rhoz;
  std::vector<TEveQuadSet *> qs_rhophi;
  std::vector<TEveBoxSet *> bs;

  TEveRGBAPalette *pal = new TEveRGBAPalette(0, 1000);

  // first, create the containers
  // clusters in 3D
  if (clusterType == "topo")
  {
    if (topoclusters_3D == nullptr)
    {
      topoclusters_3D = new TEveElementList(Form("%s clusters in 3D (no E cut)", clusterType.c_str()));
      gEve->AddElement(topoclusters_3D);
    }
    else
      topoclusters_3D->DestroyElements();
    clusters_3D = topoclusters_3D;
  }
  else
  {
    if (swclusters_3D == nullptr)
    {
      swclusters_3D = new TEveElementList(Form("%s clusters in 3D (no E cut)", clusterType.c_str()));
      gEve->AddElement(swclusters_3D);
    }
    else
      swclusters_3D->DestroyElements();
    clusters_3D = swclusters_3D;
  }

  // clusters in 2D
  if (clusterType == "topo")
  {
    if (topoclusters_rhoz == nullptr)
    {
      topoclusters_rhoz = new TEveElementList(Form("%s clusters in rho-z (E>%.1f GeV)",
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
      topoclusters_rhophi = new TEveElementList(Form("%s clusters in rho-phi (E>%.1f GeV)",
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
      swclusters_rhoz = new TEveElementList(Form("%s clusters in rho-z (E>%.1f GeV)",
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
      swclusters_rhophi = new TEveElementList(Form("%s clusters in rho-phi (E>%.1f GeV)",
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

  // now loop over the clusters and fill the containers
  // in 3D this is done in this first loop
  // for 2D projection, in this loop integrate energy over projected dimension
  // and fill map of 2D cell ID vs energy
  // then the filling is done in a second loop later
  for (unsigned int i = 0; i < nClusters; i++)
  {
    float energy, xcl, ycl, zcl;
    if (clusterType == "topo")
    {
      energy = (*eventReader->CorrectedCaloTopoClusters_energy)[i];
      xcl = (*eventReader->CorrectedCaloTopoClusters_position_x)[i];
      ycl = (*eventReader->CorrectedCaloTopoClusters_position_y)[i];
      zcl = (*eventReader->CorrectedCaloTopoClusters_position_z)[i];
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
      aqs->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",
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
    }

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
  }

  // loop over cells
  std::unordered_map<int, double> cellEnergies_rhoz;
  std::unordered_map<int, double> cellEnergies_rhophi;
  unsigned int nCells = (clusterType == "topo") ? eventReader->PositionedCaloTopoClusterCells_energy->GetSize() : eventReader->PositionedCaloClusterCells_energy->GetSize();
  for (unsigned int i = 0; i < nCells; i++)
  {
    int icl = -1;
    for (unsigned int j = 0; j < nClusters; j++)
    {
      unsigned int first_hit, last_hit;
      if (clusterType == "topo")
      {
        first_hit = (*eventReader->CorrectedCaloTopoClusters_hits_begin)[j];
        last_hit = (*eventReader->CorrectedCaloTopoClusters_hits_end)[j];
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
      cellID = (*eventReader->PositionedCaloTopoClusterCells_cellID)[i];
      energy = (*eventReader->PositionedCaloTopoClusterCells_energy)[i];
      x_center = (*eventReader->PositionedCaloTopoClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->PositionedCaloTopoClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->PositionedCaloTopoClusterCells_position_z)[i] * mm;
    }
    else
    {
      cellID = (*eventReader->PositionedCaloClusterCells_cellID)[i];
      energy = (*eventReader->PositionedCaloClusterCells_energy)[i];
      x_center = (*eventReader->PositionedCaloClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->PositionedCaloClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->PositionedCaloClusterCells_position_z)[i] * mm;
    }
    float r_center = sqrt(x_center * x_center + y_center * y_center);
    // TODO: might need to add system (at least for topoclusters) to distinguish ecal/hcal when we will have clusters from both...
    int layer = (int)DetectorGeometry::ECalBarrelLayer(cellID);
    int thetaID = (int)DetectorGeometry::ECalBarrelThetaBin(cellID);
    int moduleID = (int)DetectorGeometry::ECalBarrelModule(cellID);
    // for 2D projections, assign unique ID based either on layer-theta (for rho-Z projection)
    // or layer-module (for rho-phi projection)
    // layer = 0..11 fits in 4 bits (0..15)
    // theta = 0..800 fits in 10 bits (0..1023)
    // module = 0..1535 fits in 11 bits (0..2043)
    // cluster = assume < 512 (9 bits)
    // so both IDs can fit in a 32-bit integer
    // but could even play it safer and just use cellID and remove unused fields
    // for clus
    int rhophiID = icl + 512 * layer + 512 * 16 * moduleID;
    int rhozID = icl + 512 * layer + 512 * 16 * thetaID;

    float r_in = geomReader->r[layer];
    float r_out = geomReader->r[layer + 1];
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
    float verts3D[24];
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
      x_center = (*eventReader->PositionedCaloTopoClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->PositionedCaloTopoClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->PositionedCaloTopoClusterCells_position_z)[i] * mm;
    }
    else {
      x_center = (*eventReader->PositionedCaloClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->PositionedCaloClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->PositionedCaloClusterCells_position_z)[i] * mm;
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

  for (unsigned int i = 0; i < nCells; i++)
  {
    int icl = -1;
    for (unsigned int j = 0; j < nClusters; j++)
    {
      unsigned int first_hit, last_hit;
      if (clusterType == "topo")
      {
        first_hit = (*eventReader->CorrectedCaloTopoClusters_hits_begin)[j];
        last_hit = (*eventReader->CorrectedCaloTopoClusters_hits_end)[j];
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
      cellID = (*eventReader->PositionedCaloTopoClusterCells_cellID)[i];
    }
    else
    {
      cellID = (*eventReader->PositionedCaloClusterCells_cellID)[i];
    }

    // TODO: update when we also have HCAL in topoclustering
    int layer = (int)DetectorGeometry::ECalBarrelLayer(cellID);
    int thetaID = (int)DetectorGeometry::ECalBarrelThetaBin(cellID);
    int rhozID = icl + 512 * layer + 512 * 16 * thetaID;
    if (cellEnergies_rhoz.count(rhozID) == 0)
      continue;
    energy = cellEnergies_rhoz[rhozID];
    cellEnergies_rhoz.erase(rhozID);
    float r_in = geomReader->r[layer];
    float r_out = geomReader->r[layer + 1];
    float x_center, y_center, z_center;
    if (clusterType == "topo")
    {
      x_center = (*eventReader->PositionedCaloTopoClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->PositionedCaloTopoClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->PositionedCaloTopoClusterCells_position_z)[i] * mm;
    }
    else
    {
      x_center = (*eventReader->PositionedCaloClusterCells_position_x)[i] * mm;
      y_center = (*eventReader->PositionedCaloClusterCells_position_y)[i] * mm;
      z_center = (*eventReader->PositionedCaloClusterCells_position_z)[i] * mm;
    }
    float r_center = sqrt(x_center * x_center + y_center * y_center);
    float theta_center = atan2(r_center, z_center);
    verts[0] = r_in / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);
    verts[1] = r_in * sgn(y_center);
    verts[3] = r_in / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);
    verts[4] = r_in * sgn(y_center);
    verts[6] = r_out / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);
    verts[7] = r_out * sgn(y_center);
    verts[9] = r_out / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer] / 2.);
    verts[10] = r_out * sgn(y_center);
    verts[2] = verts[5] = verts[8] = verts[11] = 0.;
    qs_rhoz[icl]->AddQuad(verts);
    qs_rhoz[icl]->QuadValue((int)(1000 * energy));
    qs_rhoz[icl]->QuadId(new TNamed(Form("%s %d cell L%d Th%d", clusterType.c_str(), icl, layer, thetaID), "Dong!"));
  }

  // cluster cells in rho-phi projection
  for (auto &[key, energy] : cellEnergies_rhophi)
  {
    int ikey = (int)key;
    int icl = ikey % 512;
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
}

void EventDisplay::loadEvent(int event)
{

  if (event != -1)
    eventId = event;

  printf("Loading event %d ...\n", eventId);
  textEntry->SetTextColor(0xff0000);
  textEntry->SetText(Form("Loading event %d ...", eventId));

  eventReader->loadEvent(eventId);

  TString partType;
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
  if (pdgID == 111)
    partType = "pi0";
  else if (pdgID == 211)
    partType = "pi+";
  else if (pdgID == 22)
    partType = "y";
  else if (pdgID == 11)
    partType = "e-";
  else if (pdgID == -11)
    partType = "e+";
  else
    partType = "unknown";

  const double cm = geomReader->cm;
  const double mm = geomReader->mm;
  double rMax = doHCal ? geomReader->rMaxHCal : geomReader->rMax;
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
  if (drawParticles)
  {
    cout << "Creating particles" << endl;
    if (particles == nullptr)
    {
      particles = new TEveTrackList("particles");
      TEveTrackPropagator *trkProp = particles->GetPropagator();
      trkProp->SetMagField(0.01);
      trkProp->SetMaxR(rMax);
      trkProp->SetMaxZ(geomReader->zMax);
      particles->SetMainColor(kWhite);
      particles->SetLineWidth(2);
      gEve->AddElement(particles);
    }
    else
      particles->DestroyElements();

    // handle differently the e/gamma vs pi0 particle guns
    // (for pi0, need to look for the two photons in secondary particle list
    // unfortunately cannot just use the latter to show particles since
    // info about endpoint is broken (identical to vertex)
    float m = (*eventReader->genParticles_mass)[ipmax];
    float px = (*eventReader->genParticles_momentum_x)[ipmax];
    float py = (*eventReader->genParticles_momentum_y)[ipmax];
    float pz = (*eventReader->genParticles_momentum_z)[ipmax];
    float p = sqrt(px * px + py * py + pz * pz);
    float pT = sqrt(px * px + py * py);

    double t = (*eventReader->genParticles_time)[ipmax];
    double x1 = (*eventReader->genParticles_vertex_x)[ipmax] * mm;
    double y1 = (*eventReader->genParticles_vertex_y)[ipmax] * mm;
    double z1 = (*eventReader->genParticles_vertex_z)[ipmax] * mm;

    // find when particle leaves detector volume (assuming it comes from 0,0,0)
    // for zMax use ECAL which is a bit larger than HCAL
    double tmax = rMax / pT;
    /*
    if (pz!=0.0) {
      if (geomReader->zMax/fabs(pz)<tmax)
  tmax = geomReader->zMax/fabs(pz);
    }
    */
    double x2 = x1 + px * tmax;
    double y2 = y1 + py * tmax;
    double z2 = z1 + pz * tmax;

    /*
    double r1 = sqrt(x1*x1+y1*y1);
    double r2 = sqrt(x2*x2+y2*y2);
    cout << "x1 y1 z1 x2 y2 z2 = "
   << x1 << " " << y1 << " " << z1 << " "
   << x2 << " " << y2 << " " << z2 << endl;
    cout << "r1 r2 = "
         << r1 << " " << r2 << endl;
    */

    TEveMCTrack mct;
    mct.SetPdgCode(pdgID);
    mct.SetMomentum(px, py, pz, sqrt(p * p + m * m));
    mct.SetProductionVertex(x1, y1, z1, t);
    TEveTrack *track = new TEveTrack(&mct, particles->GetPropagator());
    track->SetAttLineAttMarker(particles);
    track->SetElementTitle(Form("p = %.3f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
                                p, acos(pz / p), atan2(py, px),
                                x1 / cm, y1 / cm, z1 / cm));
    particles->AddElement(track);

    // if the particle is a pi0, also draw the two photons, and set the endpoint
    // of the pi0 track
    if (pdgID == 111)
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
    particles->MakeTracks();
  }

  //
  // hits (ECAL/HCAL)
  //
  if (drawHits)
  {
    cout << "Creating hits" << endl;
    if (ecalHits == nullptr)
    {
      ecalHits = new TEvePointSet();
      ecalHits->SetName(Form("ECAL hits (E>%.1f GeV)", HitEnergyThreshold));
      ecalHits->SetMarkerStyle(4);
      ecalHits->SetMarkerSize(1);
      ecalHits->SetMarkerColor(kRed);
      gEve->AddElement(ecalHits);
    }
    else
      ecalHits->Reset();

    for (unsigned int i = 0; i < eventReader->ECalBarrelPositionedHits_position_x->GetSize(); i++)
    {
      float E = (*eventReader->ECalBarrelPositionedHits_energy)[i];
      if (E < HitEnergyThreshold)
        continue;
      // ULong_t cellID = (*eventReader->ECalBarrelPositionedHits_cellID)[i];
      // ULong_t layer = DetectorGeometry:::Layer(cellID);
      ecalHits->SetNextPoint(
          (*eventReader->ECalBarrelPositionedHits_position_x)[i] * mm,
          (*eventReader->ECalBarrelPositionedHits_position_y)[i] * mm,
          (*eventReader->ECalBarrelPositionedHits_position_z)[i] * mm);
    }

    if (doHCal)
    {
      if (hcalHits == nullptr)
      {
        hcalHits = new TEvePointSet();
        hcalHits->SetName(Form("HCAL hits (E>%.1f GeV)", HitEnergyThreshold));
        hcalHits->SetMarkerStyle(4);
        hcalHits->SetMarkerSize(1);
        hcalHits->SetMarkerColor(kRed);
        gEve->AddElement(hcalHits);
      }
      else
        hcalHits->Reset();

      for (unsigned int i = 0; i < eventReader->HCalBarrelPositionedHits_position_x->GetSize(); i++)
      {
        float E = (*eventReader->HCalBarrelPositionedHits_energy)[i];
        if (E < HitEnergyThreshold)
          continue;
        hcalHits->SetNextPoint(
            (*eventReader->HCalBarrelPositionedHits_position_x)[i] * mm,
            (*eventReader->HCalBarrelPositionedHits_position_y)[i] * mm,
            (*eventReader->HCalBarrelPositionedHits_position_z)[i] * mm);
      }
    }
  }

  //
  // cells (ECAL/HCAL)
  //
  if (drawCells)
  {
    cout << "Creating cells" << endl;
    if (ecalCells == nullptr)
    {
      ecalCells = new TEvePointSet();
      ecalCells->SetName(Form("ECAL cells (E>%.1f GeV)", CellEnergyThreshold));
      ecalCells->SetMarkerStyle(4);
      ecalCells->SetMarkerSize(3);
      ecalCells->SetMarkerColor(kYellow);
      gEve->AddElement(ecalCells);
    }
    else
      ecalCells->Reset();
    for (unsigned int i = 0; i < eventReader->ECalBarrelPositionedCells_position_x->GetSize(); i++)
    {
      float E = (*eventReader->ECalBarrelPositionedCells_energy)[i];
      if (E < CellEnergyThreshold)
        continue;
      // ULong_t cellID = (*eventReader->ECalBarrelPositionedCells_cellID)[i];
      // ULong_t layer = DetectorGeometry::Layer(cellID);
      ecalCells->SetNextPoint((*eventReader->ECalBarrelPositionedCells_position_x)[i] * mm,
                              (*eventReader->ECalBarrelPositionedCells_position_y)[i] * mm,
                              (*eventReader->ECalBarrelPositionedCells_position_z)[i] * mm);
    }

    if (doHCal)
    {
      if (hcalCells == nullptr)
      {
        hcalCells = new TEvePointSet();
        hcalCells->SetName(Form("HCAL cells (E>%.1f GeV)", CellEnergyThreshold));
        hcalCells->SetMarkerStyle(4);
        hcalCells->SetMarkerSize(3);
        hcalCells->SetMarkerColor(kYellow);
        gEve->AddElement(hcalCells);
      }
      else
        hcalCells->Reset();
      for (unsigned int i = 0; i < eventReader->HCalBarrelPositionedCells_position_x->GetSize(); i++)
      {
        float E = (*eventReader->HCalBarrelPositionedCells_energy)[i];
        if (E < CellEnergyThreshold)
          continue;
        hcalCells->SetNextPoint((*eventReader->HCalBarrelPositionedCells_position_x)[i] * mm,
                                (*eventReader->HCalBarrelPositionedCells_position_y)[i] * mm,
                                (*eventReader->HCalBarrelPositionedCells_position_z)[i] * mm);
      }
    }
  }

  //
  // cells merged (ECAL)
  //
  if (drawMergedCells)
  {
    cout << "Creating merged cells" << endl;
    if (cells_merged == nullptr)
    {
      cells_merged = new TEvePointSet();
      cells_merged->SetName("cells_merged");
      cells_merged->SetMarkerStyle(4);
      cells_merged->SetMarkerSize(3);
      cells_merged->SetMarkerColor(kBlue);
      gEve->AddElement(cells_merged);
    }
    else
      cells_merged->Reset();
    for (unsigned int i = 0; i < eventReader->ECalBarrelPositionedCells2_position_x->GetSize(); i++)
    {
      float E = (*eventReader->ECalBarrelPositionedCells2_energy)[i];
      // if (E<minCellE) continue;
      // ULong_t cellID = (*eventReader->ECalBarrelPositionedCells_cellID)[i];
      // ULong_t layer = DetectorGeometry::Layer(cellID);
      cells_merged->SetNextPoint((*eventReader->ECalBarrelPositionedCells2_position_x)[i] * mm,
                                 (*eventReader->ECalBarrelPositionedCells2_position_y)[i] * mm,
                                 (*eventReader->ECalBarrelPositionedCells2_position_z)[i] * mm);
    }
  }

  //
  // clusters
  //
  if (drawTopoClusters)
    FillClusters("topo");
  if (drawSWClusters)
    FillClusters("sw");

  //
  // event label
  //
  if (eventLabel == nullptr)
  {
    eventLabel = new TGLConstAnnotation(gEve->GetDefaultGLViewer(),
                                        Form("%s, %.1f GeV\nEvent %d",
                                             partType.Data(), pmax, eventId),
                                        0.1, 0.9);
    eventLabel->SetTextSize(0.05); // % of window diagonal
    eventLabel->SetAllowClose(false);
  }
  else
  {
    eventLabel->SetText(Form("%s, %.1f GeV\nEvent %d", partType.Data(), pmax, eventId));
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

  cout << "Done" << endl
       << endl;

  textEntry->SetTextColor((Pixel_t)0x000000);
  textEntry->SetText(Form("Event %d loaded", eventId));
}

EventDisplay::EventDisplay()
{
}

void EventDisplay::startDisplay(int initialEvent)
{

  // calculate the geometry parameters
  geomReader = new DetectorGeometry;

  cout << "******************************************************************************" << endl;
  cout << "Displaying the geometry" << endl;
  cout << "******************************************************************************" << endl
       << endl;

  // create the eve manageer
  TEveManager::Create();

  // Set title of main window
  gEve->GetBrowser()->SetWindowName("Allegro calorimeter event display");

  // see palettes here: https://root.cern.ch/doc/master/classTColor.html
  // gStyle->SetPalette(kAvocado);
  gStyle->SetPalette(kSienna);

  // first tab
  gEve->GetDefaultGLViewer()->SetGuideState(TGLUtil::kAxesOrigin, false, false, 0);
  gEve->GetDefaultGLViewer()->DrawGuides();
  gEve->GetDefaultViewer()->SetElementName("3D view");
  gEve->GetDefaultGLViewer()->CurrentCamera().RotateRad(-.7, 0.5);

  // Create the geometry and the readout
  TEveElementList *geom = new TEveElementList("Geometry");
  TEveElementList *PCBs = new TEveElementList("PCBs");
  TEveElementList *actives = new TEveElementList("Active elements");
  TEveElementList *passives = new TEveElementList("Passive elements");
  TEveElementList *readout = new TEveElementList("Readout");
  if (useG4geom)
  {
    // auto fGeom = TFile::Open(geomFile.c_str(), "CACHEREAD");
    cout << "Reading Geant4 geometry from file " << geomFile << endl;
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
      geom->AddElement(world);
      world->SetRnrSelfChildren(false, true);
    }
    else
    {

      TPRegexp re;

      re = TPRegexp("ECalBarrel*");
      TEveElement *ecalbarrel = world->FindChild(re);
      ecalbarrel->SetPickableRecursively(kTRUE);
      geom->AddElement(ecalbarrel);

      // set transparency of the subvolumes of the bath
      re = TPRegexp("LAr_bath*");
      TEveElement *bath = ecalbarrel->FindChild(re);
      TEveElement::List_t matches;
      re = TPRegexp("ECAL_Cryo*");
      ecalbarrel->FindChildren(matches, re);
      for (auto a : matches)
        a->SetMainTransparency(70);
      re = TPRegexp("services*");
      ecalbarrel->FindChildren(matches, re);
      for (auto a : matches)
        a->SetMainTransparency(70);
      // make lists of elements inside bath to turn on/off simultaneously
      if (bath)
      {
        TEveElementList *newbath = new TEveElementList("LAr_bath");
        ecalbarrel->AddElement(newbath);
        TEveElement::List_t matches;
        re = TPRegexp("PCB*");
        bath->FindChildren(matches, re);
        for (auto a : matches)
          PCBs->AddElement(a);
        newbath->AddElement(PCBs);
        TEveElement::List_t matches2;
        re = TPRegexp("active*");
        bath->FindChildren(matches2, re);
        for (auto a : matches2)
          actives->AddElement(a);
        newbath->AddElement(actives);
        TEveElement::List_t matches3;
        re = TPRegexp("passive*");
        bath->FindChildren(matches3, re);
        for (auto a : matches3)
          passives->AddElement(a);
        newbath->AddElement(passives);
        ecalbarrel->RemoveElement(bath);
        // hide elements inside bath by default because they are slow in 3D
        newbath->SetRnrSelfChildren(true, false);
      }

      // HCAL
      if (doHCal)
      {
        // add HCal envelope and subvolumes
        re = TPRegexp("HCalEnvelopeVolume*");
        TEveElement *hcalbarrel = world->FindChild(re);
        hcalbarrel->SetPickableRecursively(kTRUE);
        geom->AddElement(hcalbarrel);

        // re = TPRegexp("HCalLayerVol*");
        // hcalbarrel->FindChildren(matches, re);
        // for (auto a : matches) a->SetRnrSelfChildren(true, false);

        // set transparency of plate faces and steel
        re = TPRegexp("HCal*PlateVol*");
        hcalbarrel->FindChildren(matches, re);
        for (auto a : matches)
          a->SetMainTransparency(70);

        re = TPRegexp("HCal*Steel*");
        hcalbarrel->FindChildren(matches, re);
        for (auto a : matches)
          a->SetMainTransparency(70);

        // group together the layers so that they can be turned on/off together
        TEveElementList *hcalLayers = new TEveElementList("HCalLayers");
        hcalbarrel->AddElement(hcalLayers);
        re = TPRegexp("HCalLayerVol*");
        TEveElement::List_t matches4;
        hcalbarrel->FindChildren(matches4, re);
        for (auto a : matches4)
        {
          hcalLayers->AddElement(a);
          hcalbarrel->RemoveElement(a);
        }
      }
    }
  }
  else
  {
    cout << "Creating simplified geometry based on calculated dimensions " << endl;

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
      double phi0 = iModule * geomReader->gridPhi - geomReader->gridPhi / 12.; // small extra shift is due to finite width of element (?)
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

  gEve->AddToListTree(readout, true);

  // create second tab (R-phi view)
  rhoPhiView = gEve->SpawnNewViewer("Projection Rho-Phi");
  // two scenes, for geometry and event
  rhoPhiScene = gEve->SpawnNewScene("Rho-Phi geometry",
                                    "Scene holding projected geometry data for the RhoPhi view.");
  rhoPhiView->AddScene(rhoPhiScene);
  if (evtFile != "")
  {
    rhoPhiEventScene = gEve->SpawnNewScene("RhoPhi Event Data",
                                           "Scene holding projected event-data for the RhoPhi view.");
    rhoPhiView->AddScene(rhoPhiEventScene);
  }

  rhoPhiEventSceneManual = gEve->SpawnNewScene("RhoPhi Event Data 2",
                                               "Scene holding hand-crafted event-data for the RhoPhi view.");
  rhoPhiView->AddScene(rhoPhiEventSceneManual);
  rhoPhiGLView = rhoPhiView->GetGLViewer();
  // set camera orientation
  rhoPhiGLView->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  // create 3D->2D projection manager for rho-phi
  rhoPhiProjManager = new TEveProjectionManager();
  rhoPhiProjManager->SetProjection(TEveProjection::kPT_RPhi);
  auto axes = new TEveProjectionAxes(rhoPhiProjManager);
  axes->SetElementName("Rho-Phi projection axes");
  rhoPhiScene->AddElement(axes);

  if (useG4geom)
    rhoPhiProjManager->ImportElements(geom, rhoPhiScene);
  else
  {
    rhoPhiProjManager->ImportElements(ecalbarrel, rhoPhiScene);
    if (doHCal) rhoPhiProjManager->ImportElements(hcalbarrel, rhoPhiScene);
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
      gridmod->AddLine(x1, y1, 0., x2, y2, 0.);
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
                                  "Scene holding projected geometry data for the RhoZ view.");
  rhoZView->AddScene(rhoZScene);
  if (evtFile != "")
  {
    rhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
                                         "Scene holding projected event-data for the RhoZ view.");
    rhoZView->AddScene(rhoZEventScene);
  }
  rhoZEventSceneManual = gEve->SpawnNewScene("RhoZ Event Data 2",
                                             "Scene holding hand-crafted event-data for the RhoZ view.");
  rhoZView->AddScene(rhoZEventSceneManual);
  rhoZGLView = rhoZView->GetGLViewer();
  rhoZGLView->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);

  rhoZProjManager = new TEveProjectionManager();
  rhoZProjManager->SetProjection(TEveProjection::kPT_RhoZ);
  auto axes2 = new TEveProjectionAxes(rhoZProjManager);
  axes2->SetElementName("Rho-Z projection axes");
  rhoZScene->AddElement(axes2);
  if (useG4geom)
    rhoZProjManager->ImportElements(geom, rhoZScene);
  else
  {
    rhoZProjManager->ImportElements(ecalbarrel, rhoZScene);
    if (doHCal) rhoZProjManager->ImportElements(hcalbarrel, rhoZScene);
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

  // draw the HCAL readout segmentation in eta or theta (rho-z view)
  if (doHCal)
  {
    //TEveStraightLineSet *hcalRhoZReadout = new TEveStraightLineSet("HCAL eta readout");
    TEveStraightLineSet *hcalRhoZReadout = new TEveStraightLineSet("HCAL theta readout");
    hcalRhoZReadout->SetLineColor(kViolet);
    hcalRhoZReadout->SetLineWidth(5);
    //for (int iEta = 0; iEta <= geomReader->nEtaBinsHCal; iEta++)
    //{
    //  double eta = geomReader->etaMinHCal + iEta * geomReader->etaGridHCal;
    //  double theta = 2 * TMath::ATan(TMath::Exp(-eta));
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

  gEve->Redraw3D(true);

  //
  // display the events
  //
  if (evtFile != "")
  {

    // create the gui for event navigation
    makeGui();

    // setup the event reader
    cout << endl;
    cout << "******************************************************************************" << endl;
    cout << "Setting up the event reader" << endl;
    cout << "******************************************************************************" << endl
         << endl;

    cout << "Reading event data from file " << evtFile << endl
         << endl;
    TFile *f = TFile::Open(evtFile.c_str(), "READ");
    eventReader = new EventReader(f, doHCal, drawSWClusters, drawTopoClusters);
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

  // Draw
  gEve->Redraw3D(true);
}

/******************************************************************************/
// GUI
/******************************************************************************/

void EventDisplay::fwd()
{
  if (eventId < nEvents - 1)
  {
    ++eventId;
    loadEvent(eventId);
  }
  else
  {
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText("Already at last event");
    printf("\nAlready at last event.\n");
  }
}

void EventDisplay::bck()
{
  if (eventId > 0)
  {
    --eventId;
    loadEvent(eventId);
  }
  else
  {
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText("Already at first event");
    printf("\nAlready at first event.\n");
  }
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

    TString icondir(Form("%s/icons/", gSystem->Getenv("ROOTSYS")));
    if (!std::filesystem::exists(icondir.Data()))
    {
      icondir = Form("%s/share/root/icons/", gSystem->Getenv("ROOTSYS"));
    }

    b = new TGPictureButton(hf, gClient->GetPicture(icondir + "GoBack.gif"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "bck()");

    b = new TGPictureButton(hf, gClient->GetPicture(icondir + "GoForward.gif"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 10, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "fwd()");

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
  browser->SetTabTitle("Event Control", 0);
}
