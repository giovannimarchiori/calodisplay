#include "GraphicUtils.h"

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
