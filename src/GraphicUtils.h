#ifndef __GRAPHICUTILS__
#define __GRAPHICUTILS__

#include <TGLAnnotation.h>

// derived TGLAnnotation class that overrides MouseEnter method
// to avoid editing of annotation
class TGLConstAnnotation : public TGLAnnotation
{
public:
  TGLConstAnnotation(TGLViewerBase *parent, const char *text, Float_t posx, Float_t posy);
  Bool_t MouseEnter(TGLOvlSelectRecord & /*rec*/);
};

#endif
