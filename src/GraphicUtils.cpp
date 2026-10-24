/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class TGLConstAnnotation: GL annotation that cannot be edited with the mouse
//
/******************************************************************************/

#include "GraphicUtils.h"

#include <TGLRnrCtx.h>
#include <TGLCamera.h>
#include <TGLUtil.h>
#include <TGLViewerBase.h>
#include <TMath.h>
#include <TVirtualX.h>
#ifdef __APPLE__
// OpenGL is deprecated on macOS (but still available, and used by ROOT)
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

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

TGLTextLabel::TGLTextLabel(TGLViewerBase *parent, const std::wstring &text, Float_t posx, Float_t posy, Float_t size)
    : fParent(parent), fText(text), fPosX(posx), fPosY(posy), fTextSize(size)
{
  parent->AddOverlayElement(this);
}

TGLTextLabel::~TGLTextLabel()
{
  fParent->RemoveOverlayElement(this);
}

void TGLTextLabel::Render(TGLRnrCtx &rnrCtx)
{
  if (rnrCtx.Selection() || fText.empty())
    return;
  const TGLRect &vp = rnrCtx.RefCamera().RefViewport();
  if (vp.Width() == 0 && vp.Height() == 0)
    return;

  // as in TGLAnnotation::Render
  Float_t oldDepthRange[2];
  glGetFloatv(GL_DEPTH_RANGE, oldDepthRange);
  glDepthRange(0, 0.001);
  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  TGLCapabilitySwitch lightsOff(GL_LIGHTING, kFALSE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDisable(GL_CULL_FACE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // [0, 1] window units
  rnrCtx.ProjectionMatrixPushIdentity();
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(-1.0f, -1.0f, 0.0f);
  glScalef(2.0f, 2.0f, 1.0f);
  glTranslatef(fPosX, fPosY, 0.0f);

  // the font of the annotations (arial) has no greek letters: use FreeSans
  Int_t fs = TGLFontManager::GetFontSize(TMath::Nint(vp.Height() * fTextSize), 12, 64);
  rnrCtx.RegisterFontNoScale(fs, "FreeSans.otf", TGLFont::kTexture, fFont);
  const Float_t descent = fFont.GetDescent();
  const Float_t lineHeight = fFont.GetLineHeight();
  const Float_t sy = fTextSize / (lineHeight + descent);
  glScalef(sy / vp.Aspect(), sy, 1.0f);

  TGLUtil::Color(rnrCtx.ColorSet().Markup().GetColorIndex());
  fFont.PreRender();
  // TGLFont::Render with unicode text centres the text on the given position
  // (with the centred text alignment): shift it so that the text starts at
  // (descent, -lineHeight), as in TGLAnnotation
  Float_t llx, lly, llz, urx, ury, urz;
  fFont.BBox(fText.c_str(), llx, lly, llz, urx, ury, urz);
  const Short_t align = gVirtualX->GetTextAlign();
  gVirtualX->SetTextAlign(22);
  fFont.Render(fText.c_str(), descent + 0.5 * (urx - llx), -lineHeight + 0.5 * (ury - lly), 0., 1.);
  gVirtualX->SetTextAlign(align);
  fFont.PostRender();

  glPopMatrix();
  rnrCtx.ProjectionMatrixPop();
  glPopAttrib();
  glDepthRange(oldDepthRange[0], oldDepthRange[1]);
}
