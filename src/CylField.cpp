#include "CylField.h"

#include "Units.h"
using Units::m;

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TH2D.h>
#include <TStyle.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

CylField::CylField(const std::string& filename,
                   const std::string& treename)
  : fBr(nullptr), fBz(nullptr)
{
  std::cout << "Reading magnetic field map from file " << filename << std::endl;
  
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
  std::cout << "Tree " << treename << " has " << n << " entries" << std::endl;
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
  std::cout << "r: min, max, npoints = " << fRmin << " " << fRmax << " " << nr << std::endl;
  std::cout << "z: min, max, npoints = " << fZmin << " " << fZmax << " " << nz << std::endl;
  
  // Create histograms (regular grid assumed)
  fBr = new TH2D("hBr", "Br", nr, fRmin, fRmax, nz, fZmin, fZmax);
  fBz = new TH2D("hBz", "Bz", nr, fRmin, fRmax, nz, fZmin, fZmax);

  // Fill histograms and calculate max field
  fMaxFieldMag = 0.0;
  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    int binx = fBr->GetXaxis()->FindBin(r*m);
    int biny = fBr->GetYaxis()->FindBin(z*m);
    fBr->SetBinContent(binx, biny, (double) Br);
    fBz->SetBinContent(binx, biny, (double) Bz);
    const double mag = std::sqrt(Br*Br + Bz*Bz);
    if (mag > fMaxFieldMag) {
      fMaxFieldMag = mag;
    }
  }

  fBr->SetDirectory(0);
  fBz->SetDirectory(0);
  file.Close();

  std::cout << "DONE" << std::endl;

  PrintFieldMap();
  DrawFieldMap("field.pdf", 20, 20, 20.0);
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

void CylField::PrintFieldMap() const
{
  const int nx = fBr->GetNbinsX();
  const int ny = fBr->GetNbinsY();
  std::cout << std::setw(10) << "r"
            << std::setw(10) << "z"
            << std::setw(15) << "Br"
            << std::setw(15) << "Bz"
            << "\n";
  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      double r  = fBr->GetXaxis()->GetBinCenter(ix);
      double z  = fBr->GetYaxis()->GetBinCenter(iy);
      double Br = fBr->GetBinContent(ix, iy);
      double Bz = fBz->GetBinContent(ix, iy);
      std::cout << std::setw(10) << r
                << std::setw(10) << z
                << std::setw(15) << Br
                << std::setw(15) << Bz
                << "\n";
    }
  }
}

void CylField::DrawFieldMap(const std::string& canvasName,
                            int nArrowsR,
                            int nArrowsZ,
                            double scale) const
{
    // Create magnitude histogram with swapped axes: (z, r)
  TH2D* hMag = new TH2D("hMag", "|B|;z;r",
                        fBz->GetNbinsY(), fZmin, fZmax,   // x = z
                        fBr->GetNbinsX(), fRmin, fRmax);  // y = r
  const int nx = fBr->GetNbinsX();
  const int ny = fBr->GetNbinsY();
  
  // Fill magnitude map
  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      double r  = fBr->GetXaxis()->GetBinCenter(ix);
      double z  = fBr->GetYaxis()->GetBinCenter(iy);
      double Br = fBr->GetBinContent(ix, iy);
      double Bz = fBz->GetBinContent(ix, iy);
      double mag = std::sqrt(Br*Br + Bz*Bz);
      hMag->Fill(z, r, mag);  // swapped!
    }
  }
  
  // Draw heatmap
  TCanvas* c = new TCanvas("cField", "Field map", 800, 700);
  gStyle->SetOptStat(0);
  hMag->Draw("COLZ");
  
  // Overlay arrows
  for (int ir = 0; ir < nArrowsR; ++ir) {
    for (int iz = 0; iz < nArrowsZ; ++iz) {
      double r = fRmin + (fRmax - fRmin) * (ir + 0.5) / nArrowsR;
      double z = fZmin + (fZmax - fZmin) * (iz + 0.5) / nArrowsZ;
      double Br = fBr->Interpolate(r, z);
      double Bz = fBz->Interpolate(r, z);
      double mag = std::sqrt(Br*Br + Bz*Bz);
      if (mag < 1e-12) continue;
      
      // Direction in (z,r) plane:
      // x-axis = z → component = Bz
      // y-axis = r → component = Br
      double dz = Bz;
      double dr = Br;
      
      // Normalize direction
      double norm = std::sqrt(dz*dz + dr*dr);
      dz /= norm;
      dr /= norm;
      
      // Length proportional to magnitude
      double len = scale * mag;
      double z2 = z + dz * len;
      double r2 = r + dr * len;
      TArrow* arrow = new TArrow(z, r, z2, r2, 0.015, "|>");
      arrow->SetLineColor(kBlack);
      arrow->SetFillColor(kBlack);
      arrow->Draw();
    }
  }
  c->Update();
  c->SaveAs(canvasName.c_str());
  c->Close();
}
