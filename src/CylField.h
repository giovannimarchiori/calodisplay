#ifndef __CYL_FIELD_H__
#define __CYL_FIELD_H__

#include "TEveTrackPropagator.h"
#include <TH2D.h>
#include <string>

class CylField : public TEveMagField
{
private:
  TH2D* fBr;  // Br(r,z)
  TH2D* fBz;  // Bz(r,z)

  double fRmin, fRmax;
  double fZmin, fZmax;

  double fMaxFieldMag;

public:
  CylField(const std::string& filename,
           const std::string& treename = "tree");
  virtual ~CylField();

  // Override double-precision interface
  virtual Double_t GetMaxFieldMagD() const override;
  virtual TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t z) const override;

  void PrintFieldMap() const;
  void DrawFieldMap(const std::string& canvasName = "cField",
                    int nArrowsR = 20,
                    int nArrowsZ = 20,
                    double scale = 1.0) const;
};

#endif
