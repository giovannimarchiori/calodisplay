/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class TGLConstAnnotation: GL annotation that cannot be edited with the mouse
//
/******************************************************************************/

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
