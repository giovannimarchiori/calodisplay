/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventDisplay: creates and populates the event display
//
/******************************************************************************/

#ifndef EVENTDISPLAY_H
#define EVENTDISPLAY_H

/******************************************************************************/
// dependencies
/******************************************************************************/
#include <TEveViewer.h>
#include <TEveGeoShape.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveElement.h>
#include <TEveTrack.h>
#include <TEvePointSet.h>
#include <TEveQuadSet.h>
#include <TEveBoxSet.h>
#include <TGLViewer.h>
#include <TGLAnnotation.h>
#include <TGTextEntry.h>

#include "DetectorGeometry.h"
#include "EventReader.h"
#include "CaloCluster.h"

using namespace std;

/******************************************************************************/
// auxiliary classes/functions
/******************************************************************************/

// helper function to return the sign of a float
int sgn(float val);

// derived TGLAnnotation class that overrides MouseEnter method
// to avoid editing of annotation
class TGLConstAnnotation : public TGLAnnotation
{
public:
  TGLConstAnnotation(TGLViewerBase *parent, const char *text, Float_t posx, Float_t posy);
  Bool_t MouseEnter(TGLOvlSelectRecord & /*rec*/);
};

/******************************************************************************/
// main class
/******************************************************************************/

class EventDisplay
{

public:
  /******************************************************************************/
  // SETTINGS FOR THE DISPLAY
  // - geometry file (and corresponding merging of theta cells/modules vs layer
  // - data file
  // - flags to draw or not elements of the event (particles, hits, ..)
  // - minimum energy of particles, hits, cells, clusters
  /******************************************************************************/

  bool doHCal = false;
  bool doEndcaps = false;
  bool showFullDetector = false;

  bool drawParticles = true;
  bool drawHits = true;
  bool drawCells = true;
  bool drawMergedCells = false;
  bool drawTopoClusters = true;
  bool drawSWClusters = false;
  bool drawClustersBarycenterVsLayer = true;

  // min particle energy (GeV)
  float ParticleEnergyThreshold = 1.0;
  // min hit and cell energy (GeV)
  float HitEnergyThreshold = 0.0;
  float CellEnergyThreshold = 0.0;
  // min cluster energy (in GeV)
  float ClusterEnergyThreshold = 1.;

  // G4 geometry file
  bool useG4geom = true;
  std::string geomFile = "data/allegro.root";

  // merging along theta and module directions and drawing options
  // std::string evtFile = "output_fullCalo_SimAndDigi_withTopoCluster_MagneticField_False_pMin_10000_MeV_ThetaMinMax_40_140_pdgId_11_pythiaFalse_NoiseFalse.root";
  // std::string evtFile = "test.root";
  // const std::vector<int> mergedCells_Theta = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  // const std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  // bool drawMergedCells = false;

  // std::string evtFile = "output_fullCalo_SimAndDigi_withTopoCluster_MagneticField_False_pMin_10000_MeV_ThetaMinMax_40_140_pdgId_11_pythiaFalse_NoiseFalse_testmerge.root";
  // const std::vector<int> mergedCells_Theta = {4, 8, 2, 1, 8, 4, 2, 1, 4, 2, 1, 8};
  // const std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};
  // bool drawMergedCells = true;

  std::string evtFile = "";

  Int_t nEvents = 0; // Number of events in file
  Int_t eventId = 0; // Current event id

  /******************************************************************************/
  // ELEMENTS OF THE EVENT DISPLAY
  /******************************************************************************/
  TEveTrackList *particles = nullptr;
  TEvePointSet *ecalHits = nullptr;
  TEvePointSet *hcalHits = nullptr;
  TEvePointSet *ecalCells = nullptr;
  TEvePointSet *hcalCells = nullptr;
  TEvePointSet *cells_merged = nullptr;
  TEvePointSet *topoclustersCenter = nullptr;
  TEvePointSet *swclustersCenter = nullptr;
  TEveElementList *topoclusters_rhoz = nullptr;
  TEveElementList *topoclusters_rhophi = nullptr;
  TEveElementList *swclusters_rhoz = nullptr;
  TEveElementList *swclusters_rhophi = nullptr;
  TEveElementList *topoclusters_3D = nullptr;
  TEveElementList *swclusters_3D = nullptr;
  TEveElementList *hits = nullptr;
  TEveElementList *digis = nullptr;
  TEveGeoShape *ecalbarrel = nullptr;
  TEveElementList *ecalendcap = nullptr;
  TEveGeoShape *hcalbarrel = nullptr;
  TEveElementList *hcalendcap = nullptr;
  TEveViewer *rhoPhiView = nullptr;
  TEveViewer *rhoZView = nullptr;
  TGLViewer *rhoPhiGLView = nullptr;
  TGLViewer *rhoZGLView = nullptr;
  TEveScene *rhoPhiScene = nullptr;
  TEveScene *rhoPhiEventScene = nullptr;
  TEveScene *rhoPhiEventSceneManual = nullptr;
  TEveScene *rhoZScene = nullptr;
  TEveScene *rhoZEventScene = nullptr;
  TEveScene *rhoZEventSceneManual = nullptr;
  TEveProjectionManager *rhoPhiProjManager = nullptr;
  TEveProjectionManager *rhoZProjManager = nullptr;

  TGLConstAnnotation *eventLabel = nullptr;
  TGTextEntry *textEntry = nullptr;

  DetectorGeometry *geomReader = nullptr;
  EventReader *eventReader = nullptr;
  std::vector<CaloCluster *> topoclusterData;
  std::vector<CaloCluster *> swclusterData;

  /******************************************************************************/
  // class methods
  /******************************************************************************/

  // constructor
  EventDisplay();

  // start event display
  void startDisplay(int initialEvent = 0);

  // fill cluster objects for given collection (sw or topo) and loaded event
  // invoked by loadEvent
  void FillClusters(std::string clusterType);

  // draw the clusters created by FillClusters
  void DrawClusters(std::string clusterType);

  // load event
  void loadEvent(int event);

  // move to next event
  void fwd();

  // move to previous event
  void bck();

  // create the buttons for navigating through the events
  void makeGui();
};

#endif
