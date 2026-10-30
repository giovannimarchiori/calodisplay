/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class TGLConstAnnotation: GL annotation that cannot be edited with the mouse
// class TGLTextLabel: GL text label with unicode text
// class MarkerScaler: sets the size of the markers of a viewer from the zoom
//
/******************************************************************************/

#ifndef __GRAPHICUTILS__
#define __GRAPHICUTILS__

#include <TGLAnnotation.h>
#include <TGLFontManager.h>
#include <TGLOverlay.h>

#include <string>

class TEveElement;
class TEveScene;
class TGLViewer;

// derived TGLAnnotation class that overrides MouseEnter method
// to avoid editing of annotation
class TGLConstAnnotation : public TGLAnnotation
{
public:
  TGLConstAnnotation(TGLViewerBase *parent, const char *text, Float_t posx, Float_t posy);
  Bool_t MouseEnter(TGLOvlSelectRecord & /*rec*/);
};

// text label in a GL viewer, with unicode text (e.g. greek letters and
// superscripts, which TGLAnnotation cannot draw), in a font that has them.
// It is drawn like a TGLAnnotation without box: posx, posy ([0, 1]) is the
// top-left corner of the text, size its height as a fraction of the viewer
// height
class TGLTextLabel : public TGLOverlayElement
{
public:
  TGLTextLabel(TGLViewerBase *parent, const std::wstring &text, Float_t posx, Float_t posy, Float_t size);
  ~TGLTextLabel() override;
  void SetText(const std::wstring &text) { fText = text; }
  void Render(TGLRnrCtx &rnrCtx) override;

private:
  TGLViewerBase *fParent;
  std::wstring fText;
  Float_t fPosX;
  Float_t fPosY;
  Float_t fTextSize;
  TGLFont fFont;
};

// size of the markers of a point set: its size in the detector, converted to
// pixels with the zoom of each viewer, between a minimum and a maximum (in
// marker size units, 5 pixels each, as in TAttMarker::SetMarkerSize)
struct MarkerScaling
{
  Float_t sizeCm;
  Float_t minSize;
  Float_t maxSize;
};

// invisible overlay element of a GL viewer that sets the size of the markers
// of the point sets in the scene of the viewer (e.g. the event) that have a
// MarkerScaling (see SetScaling), so that they shrink when zooming out and
// grow when zooming in. The size is computed when the viewer is drawn, from
// its camera, and applied (when it changes by more than 10%) at the next draw,
// without drawing again (except the first time): during a zoom, the markers
// follow with one step of delay. Apply sets it at once with the last zoom,
// e.g. for a new event
class MarkerScaler : public TGLOverlayElement
{
public:
  MarkerScaler(TGLViewer *viewer, TEveScene *scene) : fViewer(viewer), fScene(scene) {}
  void Render(TGLRnrCtx &rnrCtx) override;
  // returns true if the size of some markers was changed
  bool Apply();
  // set the scaling of the markers of a point set (stored in its user data,
  // as a pointer to the given scaling, which must outlive it)
  static void SetScaling(TEveElement *pointSet, const MarkerScaling *scaling);

private:
  TGLViewer *fViewer;
  TEveScene *fScene;
  // pixels per cm at the center of the view, for the screen (0 before the
  // first draw)
  Double_t fPixelsPerCm = 0.;
  bool apply(TEveElement *element);
};

#endif
