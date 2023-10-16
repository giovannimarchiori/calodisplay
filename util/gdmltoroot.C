/******************************************************************************/
// Convert GDML to ROOT (EVE)
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// Run with
// root
// .L gdmltoroot.C+
// gdmltoroot(inFile.gdml, outFile.gdml, volname, false)
//
/******************************************************************************/


/******************************************************************************/
// dependencies
/******************************************************************************/
#include <TFile.h>
#include <TGLViewer.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TEveManager.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveGeoNode.h>
#include <iostream>
using namespace std;

/******************************************************************************/
// DISPLAY DETECTOR STORED IN GDML AND CONVERT TO ROOT EVE FORMAT)
/******************************************************************************/

void gdmltoroot(std::string inFile = "GeantDetector.gdml",
		std::string outFile = "barrel.root",
		std::string topVolume = "ECalBarrel_vol",
		bool verbose = false,
		int dim=3)
{
  TEveManager::Create();
  
  TFile::SetCacheFileDir(".");
  gGeoManager = TGeoManager::Import(inFile.c_str());
  gGeoManager->DefaultColors();
  gGeoManager->SetNsegments(200); // has no effect in tube rendering in OpenGL (overridden by OGL)
  // makes the code crash..
  //gGeoManager->SetDefaultUnits(TGeoManager::kRootUnits);
  //gGeoManager->SetDefaultUnits(TGeoManager::kG4Units);

  TEveGeoTopNode* tn = nullptr;
  if (topVolume == "world") {
    // full detector
    tn = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
    tn->SetVisLevel(3);
    gEve->AddGlobalElement(tn);
    gEve->AddToListTree(tn, true);
    tn->ExpandIntoListTreesRecursively();
    tn->SaveExtract(outFile.c_str(), topVolume.c_str(), false);
  }
  else {
    TGeoVolume* vol = gGeoManager->GetVolume(topVolume.c_str());
    if (vol==nullptr) {
      std::cerr << "Volume " << topVolume << " not found" << std::endl;
      return;
    }
    vol->SetFillColorAlpha(kBlue, 0.05);
    vol->SetTransparency(90);
    TObjArray* nodes = gGeoManager->GetTopVolume()->GetNodes();
    TObjArrayIter nodeIter(nodes);
    TGeoNode* barrel = nullptr;
    while (TGeoNode* node = (TGeoNode*) nodeIter.Next()) {
      TString name(node->GetName());
      if (verbose) std::cout << name << std::endl;
      if (name.BeginsWith(topVolume.c_str())) {
	barrel = node;
	break;
      }
    }
    if (barrel == nullptr) {
      cout << "Could not find " << topVolume << " volume in GDML" << endl;
      return;
    }
    tn = new TEveGeoTopNode(gGeoManager, barrel);
    tn->SetVisOption(1);
    tn->SetVisLevel(3);
    gEve->AddGlobalElement(tn);
    gEve->AddToListTree(tn, true);
    tn->ExpandIntoListTreesRecursively();
    tn->SetCSGExportNSeg(2000);
    tn->SaveExtract(outFile.c_str(), topVolume.c_str(), false);    
  }
  gEve->FullRedraw3D(kTRUE);
  
  // EClipType not exported to CINT (see TGLUtil.h):
  // 0 - no clip, 1 - clip plane, 2 - clip box
  auto v = gEve->GetDefaultGLViewer();
  
  if (dim==3) {
    // 3d view
    v->GetClipSet()->SetClipType(TGLClip::EType(1));
    v->CurrentCamera().RotateRad(-.7, 0.5);
  }
  else {
    // 2d view
    v->GetClipSet()->SetClipType(TGLClip::EType(0));
    v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  }
  v->DoDraw();
}
