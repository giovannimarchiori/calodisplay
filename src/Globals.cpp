/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// Globals.h/cpp: define here global variables/objects
//
/******************************************************************************/

#include "Globals.h"

DisplayConfig displayConfig;

Float_t hitsMarkerSize = 1.2;
Style_t hitsMarkerStyle = 4;
Color_t hitsMarkerColor = kRed;

Float_t digiMarkerSize = 2.0;
Style_t digiMarkerStyle = 4;
Color_t digiMarkerColor = kYellow;

// size of the markers in the detector (cm), and their min and max size (in
// marker size units, 5 pixels each), with the zoom (see MarkerScaler): the
// max size is the one without scaling
const MarkerScaling hitsMarkerScaling = {0.6, 0.6, 1.2};
const MarkerScaling digiMarkerScaling = {1.0, 1.0, 2.0};
const MarkerScaling mergedCellsMarkerScaling = {1.5, 1.5, 3.0};
const MarkerScaling clusterCenterMarkerScaling = {4.0, 2.0, 5.0};
const MarkerScaling layerBarycenterMarkerScaling = {2.0, 1.2, 3.0};

bool overrideDetectorColors = true;
