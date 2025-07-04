/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
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
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TEvePointSet.h>
#include <TEveQuadSet.h>
#include <TEveBoxSet.h>
#include <TEveStraightLineSet.h>
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
  bool drawClustersBarycenterVsLayer = true;

  bool useTransparencies = false;
  
  // will read this from config file
  // tried to make them private but code was behaving weirdly..
  // min particle energy (GeV)
  float ParticleEnergyThreshold = 0.0;
  // min hit and cell energy (GeV)
  float HitEnergyThreshold = 0.0;
  float CellEnergyThreshold = 0.0;
  // min cluster energy (in GeV)
  float ClusterEnergyThreshold = 0.0;

  // G4 geometry file
  bool useG4geom = true;
  std::string geomFile = "data/ALLEGRO_o1_v03_noDCHcells.root";
  int detectorVersion = 3;

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
  // MC particles
  TEveTrackList *particles = nullptr;
  // Reconstructed tracks
  TEveTrackList *tracks = nullptr;
  // Magnetic field
  TEveMagField* magField = nullptr;
  
  // Sim hits
  TEvePointSet *vtxHits = nullptr;
  TEvePointSet *dchHits = nullptr;
  TEvePointSet *siwrHits = nullptr;
  TEvePointSet *ecalHits = nullptr;
  TEvePointSet *hcalHits = nullptr;
  TEvePointSet *muonHits = nullptr;
  // Digitised hits
  TEvePointSet *vtxDigis = nullptr;
  TEvePointSet *dchDigis = nullptr;
  TEvePointSet *siwrDigis = nullptr;
  TEvePointSet *ecalCells = nullptr;
  TEvePointSet *hcalCells = nullptr;
  TEvePointSet *muonCells = nullptr;
  TEvePointSet *cells_merged = nullptr;

  // calorimeter cluster barycenters
  TEvePointSet *topoclustersCenter = nullptr;
  TEvePointSet *swclustersCenter = nullptr;

  // cluster directions
  TEveStraightLineSet *swclustersDirection = nullptr;
  TEveStraightLineSet *topoclustersDirection = nullptr;

  // clusters (only works for ecal barrel, needs calculation of geometry for other detectors)
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
  TGLViewer *mainGLView = nullptr;
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

  TGLViewer* activeGLViewer = nullptr;

  bool initRhoPhiView = false;
  bool initRhoZView = false;

  /******************************************************************************/
  // class methods
  /******************************************************************************/

  // constructor - not needed, use compiler-generated one
  // EventDisplay();

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

  // keep track of active GL viewer
  // void onRhoPhiViewActivated();
  // void onRhoZViewActivated();
  // void on3DViewActivated();
  void onTabSelected(Int_t tab);

  void takeScreenshot();

};

#endif
