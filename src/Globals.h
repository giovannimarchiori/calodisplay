/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// Globals.h/cpp: define here global variables/objects
//
/******************************************************************************/

#ifndef GLOBALS_H
#define GLOBALS_H

#include "DisplayConfig.h"
#include "GraphicUtils.h"

#include "TStyle.h"
#include "TColor.h"

extern DisplayConfig displayConfig;

extern Float_t hitsMarkerSize;
extern Style_t hitsMarkerStyle;
extern Color_t hitsMarkerColor;

extern Float_t digiMarkerSize;
extern Style_t digiMarkerStyle;
extern Color_t digiMarkerColor;

// size of the markers in the detector, with the zoom (see MarkerScaler)
extern const MarkerScaling hitsMarkerScaling;
extern const MarkerScaling digiMarkerScaling;
extern const MarkerScaling mergedCellsMarkerScaling;
extern const MarkerScaling clusterCenterMarkerScaling;
extern const MarkerScaling layerBarycenterMarkerScaling;

extern bool overrideDetectorColors;

#endif
