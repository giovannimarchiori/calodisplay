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
  : m_magnetIsOn(true), m_reverse(false), fBr(nullptr), fBz(nullptr)
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
  fNr = nr;
  fNz = nz;
  fRmin = rvals.front();
  fRmax = rvals.back();
  fZmin = zvals.front();
  fZmax = zvals.back();
  fDr = (fRmax - fRmin) / (nr - 1);
  fDz = (fZmax - fZmin) / (nz - 1);
  fMirrorZ = (fZmin >= 0.);
  std::cout << "r: min, max, npoints = " << fRmin << " " << fRmax << " " << nr << std::endl;
  std::cout << "z: min, max, npoints = " << fZmin << " " << fZmax << " " << nz << std::endl;
  if (fMirrorZ)
    std::cout << "Map covers only z>=0: field at z<0 obtained by symmetry" << std::endl;

  // Create histograms with bin centres on the grid points (regular grid assumed)
  fBr = new TH2D("hBr", "Br", nr, fRmin - fDr/2., fRmax + fDr/2., nz, fZmin - fDz/2., fZmax + fDz/2.);
  fBz = new TH2D("hBz", "Bz", nr, fRmin - fDr/2., fRmax + fDr/2., nz, fZmin - fDz/2., fZmax + fDz/2.);

  // Fill grid and histograms and calculate max field
  fBrGrid.assign(nr * nz, 0.);
  fBzGrid.assign(nr * nz, 0.);
  std::vector<bool> filled(nr * nz, false);
  fMaxFieldMag = 0.0;
  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    const double rr = (double) r*m;
    const double zz = (double) z*m;
    const long ir = std::lround((rr - fRmin) / fDr);
    const long iz = std::lround((zz - fZmin) / fDz);
    if (std::fabs(rr - (fRmin + ir*fDr)) > 1e-3*fDr || std::fabs(zz - (fZmin + iz*fDz)) > 1e-3*fDz) {
      throw std::runtime_error("Field map points are not on a regular grid");
    }
    fBrGrid[ir*nz + iz] = Br;
    fBzGrid[ir*nz + iz] = Bz;
    filled[ir*nz + iz] = true;
    fBr->SetBinContent(ir+1, iz+1, (double) Br);
    fBz->SetBinContent(ir+1, iz+1, (double) Bz);
    const double mag = std::sqrt(Br*Br + Bz*Bz);
    if (mag > fMaxFieldMag) {
      fMaxFieldMag = mag;
    }
  }
  if (std::find(filled.begin(), filled.end(), false) != filled.end()) {
    throw std::runtime_error("Field map grid is incomplete");
  }

  fBr->SetDirectory(0);
  fBz->SetDirectory(0);
  file.Close();

  std::cout << "DONE" << std::endl;
}

CylField::~CylField()
{
  delete fBr;
  delete fBz;
}

void CylField::setMagnetState(bool state)
{
  if (state != m_magnetIsOn)
    {
      if ( state )
	std::cout << "Magnet state is changed to ON" << std::endl;
      else
	std::cout << "Magnet state is changed to OFF" << std::endl;
    }
  m_magnetIsOn = state;
}

bool CylField::isMagnetOn() const
{
  return m_magnetIsOn;
}

void CylField::setReverseState(bool state)
{
  m_reverse = state;
}

bool CylField::isReverse() const
{
  return m_reverse;
}

Double_t CylField::GetMaxFieldMagD() const
{
  return fMaxFieldMag;
}

TEveVectorD CylField::GetFieldD(Double_t x, Double_t y, Double_t z) const
{
  if (!m_magnetIsOn)
    return TEveVectorD(0,0,0);

  const double r = std::sqrt(x*x + y*y);

  // use symmetry for z<0 if the map only covers z>=0
  double zMap = z;
  double signBr = 1.0;
  if (fMirrorZ && z < 0.) {
    zMap = -z;
    signBr = -1.0;
  }

  // Outside map return zero field
  if (r < fRmin || r > fRmax || zMap < fZmin || zMap > fZmax) {
    return TEveVectorD(0.0, 0.0, 0.0);
  }

  // bilinear interpolation between the 4 surrounding grid points
  const double ur = (r - fRmin) / fDr;
  const double uz = (zMap - fZmin) / fDz;
  const int ir = std::min((int) ur, fNr - 2);
  const int iz = std::min((int) uz, fNz - 2);
  const double tr = ur - ir;
  const double tz = uz - iz;
  const int i00 = ir*fNz + iz;
  const int i10 = i00 + fNz;
  const double w00 = (1-tr)*(1-tz), w10 = tr*(1-tz), w01 = (1-tr)*tz, w11 = tr*tz;
  const double Br = signBr * (w00*fBrGrid[i00] + w10*fBrGrid[i10] + w01*fBrGrid[i00+1] + w11*fBrGrid[i10+1]);
  const double Bz = w00*fBzGrid[i00] + w10*fBzGrid[i10] + w01*fBzGrid[i00+1] + w11*fBzGrid[i10+1];

  double Bx = 0.0;
  double By = 0.0;

  // Avoid division by zero near axis

  if (r > 1e-12) {
    const double invr = 1.0 / r;
    Bx = Br * x * invr;
    By = Br * y * invr;
  }

  if (m_reverse)
    return TEveVectorD(-Bx, -By, -Bz);
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
                        fBz->GetNbinsY(), fZmin - fDz/2., fZmax + fDz/2.,   // x = z
                        fBr->GetNbinsX(), fRmin - fDr/2., fRmax + fDr/2.);  // y = r
  // owned by the canvas (deleted with it), not by the current directory
  hMag->SetDirectory(nullptr);
  hMag->SetBit(kCanDelete);
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
      arrow->SetBit(kCanDelete);
      arrow->Draw();
    }
  }
  c->Update();
  c->SaveAs(canvasName.c_str());
  // deletes also the histogram and the arrows (kCanDelete)
  delete c;
}
