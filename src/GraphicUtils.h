/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class TGLConstAnnotation: GL annotation that cannot be edited with the mouse
// class TGLTextLabel: GL text label with unicode text
//
/******************************************************************************/

#ifndef __GRAPHICUTILS__
#define __GRAPHICUTILS__

#include <TGLAnnotation.h>
#include <TGLFontManager.h>
#include <TGLOverlay.h>

#include <string>

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

#endif
