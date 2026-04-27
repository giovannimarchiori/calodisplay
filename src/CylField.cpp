#include "CylField.h"

#include "Units.h"
using Units::m;

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

CylField::CylField(const std::string& filename,
                   const std::string& treename)
  : fBr(nullptr), fBz(nullptr)
{
  TFile file(filename.c_str(), "READ");
  if (file.IsZombie()) {
    throw std::runtime_error("Cannot open ROOT file");
  }

  TTree* tree = dynamic_cast<TTree*>(file.Get(treename.c_str()));
  if (!tree) {
    throw std::runtime_error("Cannot find TTree");
  }

  float r, z, Br, Bz;
  tree->SetBranchAddress("r",  &r);
  tree->SetBranchAddress("z",  &z);
  tree->SetBranchAddress("Br", &Br);
  tree->SetBranchAddress("Bz", &Bz);

  // Collect unique sorted grid coordinates
  std::vector<double> rvals, zvals;
  Long64_t n = tree->GetEntries();
  rvals.reserve(n);
  zvals.reserve(n);

  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    rvals.push_back((double) r*m); // the ROOT file contains r, z in meters
    zvals.push_back((double) z*m); // the ROOT file contains r, z in meters
  }

  std::sort(rvals.begin(), rvals.end());
  rvals.erase(std::unique(rvals.begin(), rvals.end()), rvals.end());
  std::sort(zvals.begin(), zvals.end());
  zvals.erase(std::unique(zvals.begin(), zvals.end()), zvals.end());

  const int nr = rvals.size();
  const int nz = zvals.size();

  if (nr < 2 || nz < 2) {
    throw std::runtime_error("Grid must have at least 2x2 points");
  }

  fRmin = rvals.front();
  fRmax = rvals.back();
  fZmin = zvals.front();
  fZmax = zvals.back();

  // Create histograms (regular grid assumed)
  fBr = new TH2D("hBr", "Br", nr, fRmin, fRmax, nz, fZmin, fZmax);
  fBz = new TH2D("hBz", "Bz", nr, fRmin, fRmax, nz, fZmin, fZmax);

  // Fill histograms and calculate max field
  fMaxFieldMag = 0.0;
  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    int binx = fBr->GetXaxis()->FindBin(r);
    int biny = fBr->GetYaxis()->FindBin(z);
    fBr->SetBinContent(binx, biny, (double) Br);
    fBz->SetBinContent(binx, biny, (double) Bz);
    const double mag = std::sqrt(Br*Br + Bz*Bz);
    if (mag > fMaxFieldMag) {
      fMaxFieldMag = mag;
    }
  }

  file.Close();
}

CylField::~CylField()
{
  delete fBr;
  delete fBz;
}

Double_t CylField::GetMaxFieldMagD() const
{
  return fMaxFieldMag;
}

TEveVectorD CylField::GetFieldD(Double_t x, Double_t y, Double_t z) const
{
  const double r = std::sqrt(x*x + y*y);

  // Outside map return zero field
  if (r < fRmin || r > fRmax || z < fZmin || z > fZmax) {
    return TEveVectorD(0.0, 0.0, 0.0);
  }

  const double Br = fBr->Interpolate(r, z);
  const double Bz = fBz->Interpolate(r, z);

  double Bx = 0.0;
  double By = 0.0;

  // Avoid division by zero near axis

  if (r > 1e-12) {
    const double invr = 1.0 / r;
    Bx = Br * x * invr;
    By = Br * y * invr;
  }

  return TEveVectorD(Bx, By, Bz);
}
