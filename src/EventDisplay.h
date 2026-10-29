/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
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
#include <TEveBrowser.h>
#include <TEveGeoShape.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveElement.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TEvePointSet.h>
#include <TEveQuadSet.h>
#include <TEveRGBAPalette.h>
#include <TEveBoxSet.h>
#include <TEveStraightLineSet.h>

#include <TGLViewer.h>

#include <TGNumberEntry.h>
#include <TGLabel.h>

#include <TVector3.h>
#include <TPRegexp.h>

#include <memory>
#include <set>

#include "DetectorGeometry.h"
#include "GeometryCache.h"
#include "EventReader.h"
#include "CaloCluster.h"
#include "GraphicUtils.h"
#include "Globals.h"


/******************************************************************************/
// auxiliary classes/functions
/******************************************************************************/

// struct to hold start/end/momentum info to calculate e- brem from MCParticles
struct ElectronSegment
{
  TVector3 start;
  TVector3 end;
  TVector3 momentum;
};

// derived TEvePointSet class that shows point info in console when selected
class MyPointSet : public TEvePointSet {
public:
  std::vector<std::string> fTooltips;  // per-point info
  bool fShowTooltips;                  // flag to enable/disable tooltips

  MyPointSet(const char* name = "MyPoints", bool showTooltips = true)
    : TEvePointSet(name), fShowTooltips(showTooltips) {
    if (showTooltips) {
        SetPickable(kTRUE);
        SetOwnIds(kTRUE);
      }
  }

  void Reset(Int_t n_points = 0, Int_t n_int_ids = 0) {
    // 1. call base class Reset
    TEvePointSet::Reset(n_points, n_int_ids);

    // 2. restore your desired behaviour
    if (fShowTooltips) {
      SetOwnIds(kTRUE);
      SetPickable(kTRUE);
    }

    // 3. resets the tooltips vector
    fTooltips.clear();
  }

  void PointSelected(Int_t idx) override {
    if (!fShowTooltips) return;
    std::cout << "tooltips size: " << fTooltips.size() << std::endl;
    if (idx < 0 || idx >= (Int_t)fTooltips.size()) return;
    std::cout << "Point selected: idx = " << idx << ", " << fTooltips[idx] << std::endl;
  }
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

  bool debug = false;
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
  // max radius for propagator
  float rMax = 0.0;

  // G4 geometry file
  bool useG4geom = true;
  std::string geomFile = "data/ALLEGRO_o1_v03_noDCHcells.root";
  int detectorVersion = 3;
  // read the G4 geometry without the volumes that are never drawn from a
  // cache (see GeometryCache.h), writing it the first time
  bool useGeometryCache = true;

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
  // Magnetic field: owned by the display and shared by the track propagators
  // (which must not own it, or it would be deleted by each of them)
  std::unique_ptr<TEveMagField> magField;

  // MC particles
  TEveTrackList *particles = nullptr;
  
  // Sim hits
  MyPointSet *vtxHits = nullptr;
  MyPointSet *mainTrkHits = nullptr; // drift chamber (DCH) or straw tube tracker (STT)
  MyPointSet *siwrHits = nullptr;
  MyPointSet *ecalHits = nullptr;
  MyPointSet *hcalHits = nullptr;
  MyPointSet *muonHits = nullptr;
  TEveElementList *hits = nullptr;

  // Digitised hits
  MyPointSet *vtxDigis = nullptr;
  MyPointSet *mainTrkDigis = nullptr;
  MyPointSet *siwrDigis = nullptr;
  MyPointSet *ecalCells = nullptr;
  MyPointSet *hcalCells = nullptr;
  MyPointSet *muonCells = nullptr;
  MyPointSet *cells_merged = nullptr;
  TEveElementList *digis = nullptr;

  // depth of the event in the projections, in front of the geometry (both are
  // otherwise at the same depth, and the filled volumes can cover the tracks,
  // hits and cluster cells, e.g. the drift chamber in the screenshots of the
  // rho-z view): of the projected event and of the cluster cells drawn by hand
  static constexpr float projectedEventDepth = 10.;

  // Clusters (a folder for each set of clusters, see ClusterSet)
  TEveElementList *clustersFolder = nullptr;

  // Reconstructed tracks
  TEveTrackList *tracks = nullptr;

  // clusters: the graphic elements of the clusters of a type (topo or
  // sliding window); the cells are drawn only for the ECAL and HCAL barrel
  // final window of the sliding-window clusters: nTheta x nPhi towers of size
  // deltaTheta x deltaPhi, centred on the tower containing the cluster
  // position, with tower indices (thetaMax - theta)/deltaTheta and
  // (phi - phiMin)/deltaPhi; read from the options of the job that produced
  // the clusters, stored in the event file (nTheta = 0 if not found)
  struct SlidingWindow
  {
    int nTheta = 0;
    int nPhi = 0;
    double deltaTheta = 0.;
    double deltaPhi = 0.;
    double thetaMax = 0.;
    double phiMin = 0.;
    // radial range of the barrel calorimeters used for the towers
    // (0 if none)
    double rMin = 0.;
    double rMax = 0.;
    // |z| range of the endcap calorimeters used for the towers (0 if none)
    double zMin = 0.;
    double zMax = 0.;
  };
  struct ClusterSet
  {
    std::string type;    // used in the names of the elements
    std::string label;   // used in the printouts
    Color_t centerColor; // colour of the barycenters
    bool drawWindow;     // draw the window of the sliding-window clusters
    SlidingWindow window;
    // the clusters in the event reader (set when the event file is opened)
    const EventReader::Clusters *clusters = nullptr;
    // the clusters above threshold, with their energy and barycenter per layer
    std::vector<std::unique_ptr<CaloCluster>> data;
    // the lists of the clusters in 3D and in the projections (drawn by hand,
    // added to the projected replicas of list3D: see attachClusterProjections)
    TEveElementList *list3D = nullptr;
    TEveElementList *listRhoZ = nullptr;
    TEveElementList *listRhoPhi = nullptr;
    // the barycenters and directions of the clusters in 3D
    TEvePointSet *center = nullptr;
    TEveStraightLineSet *direction = nullptr;
  };
  // one set for each collection of topo and sliding-window clusters
  // (set when the event file is opened)
  std::vector<ClusterSet> clusterSets;

  //
  TEveGeoShape *ecalbarrel = nullptr;
  TEveElementList *ecalendcap = nullptr;
  TEveGeoShape *hcalbarrel = nullptr;
  TEveElementList *hcalendcap = nullptr;

  // the various views/scenes/projection managers
  TEveViewer *rhoPhiView = nullptr;
  TEveViewer *rhoZView = nullptr;
  TGLViewer *rhoPhiGLView = nullptr;
  TGLViewer *rhoZGLView = nullptr;
  TGLViewer *mainGLView = nullptr;
  TEveScene *rhoPhiScene = nullptr;
  TEveScene *rhoPhiEventScene = nullptr;
  TEveScene *rhoZScene = nullptr;
  TEveScene *rhoZEventScene = nullptr;
  TEveProjectionManager *rhoPhiProjManager = nullptr;
  TEveProjectionManager *rhoZProjManager = nullptr;

  // the list of the geometry
  TEveElementList* geom = nullptr;
  // the full G4 geometry, kept until the display is set up to write the
  // geometry cache, and the elements created from it
  TEveGeoShapeExtract* geometryExtract = nullptr;
  GeometryCache::ExtractMap geometryExtracts;
  // identifies the display options that change which volumes are drawn
  std::string geometryCacheTag() const;

  // the list of the readout elements
  TEveElementList *readout = nullptr;
  // the readout elements not shown when the readout is switched on (e.g. the
  // unmerged ECAL theta grid), but which can be shown from the browser
  std::set<TEveElement *> readoutsHiddenByDefault;
  TGLTextLabel *processLabel = nullptr;
  // palettes for the energy of the cells in the clusters (values in MeV) in
  // each view, shared by all events, with a legend in the projections; their
  // range is set from the largest cell energy of the event in the view
  TEveRGBAPalette *cellPalette3D = nullptr;
  TEveRGBAPalette *cellPaletteRhoPhi = nullptr;
  TEveRGBAPalette *cellPaletteRhoZ = nullptr;
  int maxCellValue3D = 0;
  int maxCellValueRhoPhi = 0;
  int maxCellValueRhoZ = 0;
  // the legends of the palettes and their titles, that can be hidden
  std::vector<TGLOverlayElement *> paletteLegends;
  bool showPaletteLegends = true;
  TGLTextLabel *eventLabel = nullptr;
  // the event number entries and the labels with the range of the event
  // numbers in the navigation bar (one per tab)
  std::vector<TGNumberEntry *> eventEntries;
  std::vector<TGLabel *> eventRangeLabels;
  // the buttons to show/hide the readouts and the legends of the cell
  // energies (one per tab), shown pressed when these are shown
  std::vector<TGButton *> readoutButtons;
  std::vector<TGButton *> paletteLegendButtons;
  void updateToggleButtons();
  // the bars with the buttons, at the bottom of the tabs
  std::vector<TGCompositeFrame *> buttonBars;
  // when the main window is closed: delete the bars of the buttons before
  // TEve closes the tabs, which destroys the windows of the bars, but not the
  // bars (their pending redraws, e.g. of the labels in the hidden tabs, would
  // then fail and crash), and quit the application once TEve is terminated
  void onCloseWindow();

  std::unique_ptr<DetectorGeometry> geomReader;
  // the file must be declared before the reader, so that the reader is destroyed first
  std::unique_ptr<TFile> eventFile;
  std::unique_ptr<EventReader> eventReader;

  TGLViewer* activeGLViewer = nullptr;

  TEveBrowser* browser = nullptr;

  bool initRhoPhiView = false;
  bool initRhoZView = false;
  // the geometry and the readout have been projected in the views
  bool rhoPhiProjectionBuilt = false;
  bool rhoZProjectionBuilt = false;

  bool showReadout = true;

  /******************************************************************************/
  // class methods
  /******************************************************************************/

  // constructor - not needed, use compiler-generated one
  // EventDisplay();

  // start event display
  void startDisplay(int initialEvent = 0);

  // fill cluster objects for given collection (sw or topo) and loaded event
  // invoked by drawEvent
  void fillClusters(ClusterSet &set);

  // draw the MC particles
  void drawParticles();

  // draw the hits
  void drawHits();

  // draw the digis
  void drawDigis();

  // draw the reconstructed tracks
  void drawTracks();
  
  // draw the clusters created by fillClusters
  void drawClusters(ClusterSet &set);
  // add the projections of the clusters to the projected replicas of their
  // 3D lists, after the event has been projected
  void attachClusterProjections();

  // draw the event
  void drawEvent(int event);

  // save the views (3d, rhophi, rhoz) of the events [firstEvent, firstEvent+nEvents)
  // to <prefix>_ev<event>_<view>.png, with the given width in pixels
  void saveScreenshots(const std::string &prefix, int firstEvent, int nEvents,
                       const std::vector<std::string> &views, int width);

  // framing of the views in the screenshots: automatic for the projections,
  // on the most energetic cluster ("cluster"), charged particle ("track") or
  // both ("track+cluster"),
  // with the framed box scaled by cameraFocusScale (0 for the default of the
  // focus) filling the image, or for some views set explicitly
  // (cameraSettings, see --camera)
  std::string cameraFocus = "all";
  double cameraFocusScale = 0.;
  std::map<std::string, std::vector<double>> cameraSettings;
  // set the camera of a view (3d, rhophi, rhoz) for the current event
  void setCamera(const std::string &view);
  // print the camera of the selected view, as a --camera option
  void printCamera();

  // move to next event
  void fwd();

  // move to previous event
  void bck();

  // move to the event typed in the event number entry
  void goToEvent();

  // move to the next/previous event with the up/down arrows of the event
  // number entry (val: the argument of TGNumberEntry::ValueChanged)
  void stepEvent(Long_t val);

  // show the current event and the range of the event numbers in the navigation bar
  void updateEventEntries();

  // show/hide the readouts
  void toggleReadout();

  // show/hide the legends of the cluster cell energies
  void togglePaletteLegends();

  // the bounding box of the object to frame (see cameraFocus) in the
  // coordinates of a view; returns false if there is none
  bool getFocusBox(const std::string &view, double boxMin[3], double boxMax[3]);
  // set an orthographic (2D) camera to show the region centred on (x, y) with
  // the given zoom, or a perspective (3D) camera at eye, looking at center,
  // with the given vertical direction and field of view
  void setOrthoCamera(TGLViewer *viewer, double x, double y, double zoom);
  void setPerspectiveCamera(TGLViewer *viewer, const double center[3], const double eye[3],
                            const double up[3], double fov);

  // create the buttons for navigating through the events
  void makeGui();

  // keep track of active GL viewer
  // void onRhoPhiViewActivated();
  // void onRhoZViewActivated();
  // void on3DViewActivated();
  void onTabSelected(Int_t tab);

  // take a screenshot
  void takeScreenshot();

  // find path marks for electron brem to make more realistic rendering of electron trajectories
  std::vector<TEvePathMarkD> getElectronBremsstrahlungMarks(
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
    const TTreeReaderArray<UInt_t>& daughters_end);

  std::vector<ElectronSegment> getElectronBremsstrahlungSegments(
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
    const TTreeReaderArray<UInt_t>& daughters_end);

  void createElementList(TEveElementList*& list, const std::string& listName);

  // draw a collection of hits, digis or cells with energy above Ethr
  // (negative for no threshold in the name of the collection)
  template <typename T> void drawHitCollection(MyPointSet*& hitSet, const std::string& name,
                                               const EventReader::HitCollection<T>& hitCollection,
                                               const float Ethr,
                                               bool digi = false,    // if true, use style for digis instead of hits
                                               bool addTooltips = false,  // if true, add tooltips (do it only for debug!)
                                               bool append = false)  // if append is true, and the pointer is not null, the collection is not reset
  {
    if (hitSet == nullptr)
    {
      hitSet = new MyPointSet("MyPointSet", addTooltips && digi && debug);
      if (Ethr>=0.0)
        hitSet->SetName(Form("%s (E>%.2f GeV)", name.c_str(), Ethr));
      else
        hitSet->SetName(name.c_str());

      if (digi) {
        hitSet->SetMarkerStyle(digiMarkerStyle);
        hitSet->SetMarkerSize(digiMarkerSize);
        hitSet->SetMarkerColor(digiMarkerColor);
        digis->AddElement(hitSet);
      }
      else {
        hitSet->SetMarkerStyle(hitsMarkerStyle);
        hitSet->SetMarkerSize(hitsMarkerSize);
        hitSet->SetMarkerColor(hitsMarkerColor);
        hits->AddElement(hitSet);
      }
    }
    else {
      if (!append)
        hitSet->Reset();
    }

    for (unsigned int i = 0; i < hitCollection.size(); i++)
    {
      // skip hit if E<threshold
      if (hitCollection.energy && (*hitCollection.energy)[i] < Ethr)
        continue;
      hitSet->SetNextPoint((*hitCollection.x)[i] * Units::mm,
                           (*hitCollection.y)[i] * Units::mm,
                           (*hitCollection.z)[i] * Units::mm);
    }
  }

private:
  struct SubDetector
  {
    std::string name;
    std::string regex;  // regular expression for the volume name (empty: use contains)
    std::vector<std::string> contains;
    TEveElementList* list = nullptr;
    bool visibleByDefault = true;
    std::function<void(TEveElement*, SubDetector&)> handler;
    // optionally
    // DetectorType type;
  };
  std::vector<SubDetector> m_subDetectors;

  SubDetector* findSubDetector(TString const& s);
  // the list of the elements of a subdetector, created at the first call
  TEveElementList *subDetectorList(SubDetector &det);
  // style of a volume in 3D: drawn without its children, with the transparency
  // of the detectors, the colour of the subdetector (if overrideDetectorColors),
  // nSegments segments (if > 0) and without frame
  void setVolumeStyle(TEveElement *el, Color_t color, int nSegments);
  // in the projections: draw the content of an envelope instead of the
  // envelope, with setChild called for each of its children
  void showEnvelopeContent(TEveElement *envelope,
                           const std::function<void(TEveElement *, const TString &)> &setChild);
  // in the projections: the layers of the vertex barrel, without their envelopes
  void showVertexBarrelLayers(TEveElement *vertexBarrel);
  // in the projections: the active elements in the LAr bath of the ECAL barrel
  // (calorimeters only mode), optionally without the cryostat side
  void showECalBathCells(TEveElement *projectedGeometry, bool hideCryoSide);
  // in the projections: the content of the ECAL barrel (full detector), with
  // optionally transparent active layers
  void showECalBarrelContent(TEveElement *ecalBarrel, bool transparentActiveLayers);
  // in the projected ECAL endcap, show the cryostat and, in the LAr bath, the
  // support tubes and the active volumes (without their content)
  void showECalEndcapContent(TEveElement *ecalEndcap);
  // in the projections: the layers of an envelope of the muon tagger
  void showMuonLayers(TEveElement *envelope);

  // the steps of startDisplay
  // read the energy thresholds and the detector version
  void readConfiguration();
  // create the magnetic field (ideal solenoid or field map, see bField)
  void createMagneticField();
  // create the eve manager, the main window and the 3D view
  void createMainWindow();
  // create the geometry (returns false if it cannot be read)
  bool buildGeometry();
  // the G4 geometry, from the geometry cache or the geometry file
  bool buildG4Geometry();
  // the full detector (--fulldet), with the handlers of the subdetectors
  void buildFullDetector(TEveElement *world);
  // the calorimeters only (returns false if the ECAL barrel is not found)
  bool buildCalorimeters(TEveElement *world);
  // the simplified geometry, from the calculated dimensions (--nog4)
  void buildSimplifiedGeometry();
  // the rho-phi and rho-z views: the projection of the geometry and the readout
  void createRhoPhiView();
  // project the geometry and the readout in the views (once, when they are
  // first shown), showing the progress in the status bar; the projected event
  // (eventScene) is added to the view after the geometry
  void buildRhoPhiProjection();
  void buildRhoZProjection();
  void buildProjection(const char *view, TEveViewer *viewer, TEveScene *scene, TEveScene *eventScene,
                       const std::function<void()> &build);
  void projectGeometryRhoPhi();
  void drawRhoPhiReadout();
  void createRhoZView();
  void projectGeometryRhoZ();
  void drawRhoZReadout();
  // write the geometry cache (when the geometry was read from the geometry file)
  void writeGeometryCache();
  // open the event file and set up the event reader
  void openEventFile();
  // legends of the cluster cell energies in the projections
  void addPaletteLegends();
  // set the 3D view as the active tab, rotate the camera and set the clipping
  void setInitialCamera();

  // decode the cellIDs with the CellIDEncoding stored in the event file
  void setCellIDEncodings();

  // set the window of a set of sliding-window clusters
  void setSlidingWindow(ClusterSet &set);
  // theta and phi range of the final window of a cluster
  void getSlidingWindowRange(const SlidingWindow &window, float xcl, float ycl, float zcl,
                             double &thetaLow, double &thetaHigh,
                             double &phiLow, double &phiHigh) const;
  // draw the final window of a sliding-window cluster in the projections
  void drawSlidingWindow(const SlidingWindow &window, float xcl, float ycl, float zcl, const char *name,
                         TEveElement *rhoz, TEveElement *rhophi);

  void addSimpleElement(TEveElement* element, SubDetector&);
  void addVertex(TEveElement* element, SubDetector&);
  void addDCH(TEveElement* element, SubDetector&);
  void addWrapper(TEveElement* element, SubDetector&);
  void addECalB(TEveElement* element, SubDetector&);
  void addECalEC(TEveElement* element, SubDetector&);
  void addHCalB(TEveElement* element, SubDetector&);
  void addHCalEC(TEveElement* element, SubDetector&);
  void addMuonB(TEveElement* element, SubDetector&);
  void addMuonEC(TEveElement* element, SubDetector&);

};

#endif
