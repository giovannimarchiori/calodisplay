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

bool overrideDetectorColors = true;
