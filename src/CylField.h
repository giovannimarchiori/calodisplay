#ifndef __CYL_FIELD_H__
#define __CYL_FIELD_H__

#include "TEveTrackPropagator.h"
#include <TH2D.h>
#include <string>
#include <vector>

class CylField : public TEveMagField
{
private:
  bool m_magnetIsOn;
  bool m_reverse;

  TH2D* fBr;  // Br(r,z), bin centres on the grid points (for printing/drawing)
  TH2D* fBz;  // Bz(r,z), bin centres on the grid points (for printing/drawing)

  // regular grid used for the field lookup, indexed as [ir*fNz + iz]
  std::vector<double> fBrGrid;
  std::vector<double> fBzGrid;
  int fNr, fNz;
  double fDr, fDz;

  double fRmin, fRmax;
  double fZmin, fZmax;

  // if the map only covers z>=0, use Br(r,-z) = -Br(r,z), Bz(r,-z) = Bz(r,z)
  bool fMirrorZ;

  double fMaxFieldMag;

public:
  CylField(const std::string& filename,
           const std::string& treename = "tree");
  virtual ~CylField();

  void setMagnetState(bool state);
  bool isMagnetOn() const;
  void setReverseState(bool state);
  bool isReverse() const;

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
