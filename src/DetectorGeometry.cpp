/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class DetectorGeometry: contains detector geometry information
//
/******************************************************************************/

#include "DetectorGeometry.h"

#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
using std::cout;
using std::endl;
using Units::units;
using Units::cm;
using Units::mm;

DetectorGeometry::DetectorGeometry(int version) {
  calcGeom(version);
}


/******************************************************************************/
// GEOMETRY HELPER FUNCTIONS AND DERIVED GEOMETRY PARAMETERS
/******************************************************************************/

// calculate length along electrode at r=r_out for electrode
// starting at r=r_in and inclined in phi by angle
// for r_in>r_out, when the solution exists, it's always > 0
double DetectorGeometry::getL(double angle, double r_in, double r_out) {
  if (r_out > r_in)
    return sqrt(r_out*r_out - r_in*r_in*sin(angle)*sin(angle)) - r_in*cos(angle);
  else
    return r_in*cos(angle) - sqrt(r_out*r_out - r_in*r_in*sin(angle)*sin(angle));
}

// calculate phi offset between point at r=r_out with respect 
// to point at r=r_in if they are on the same electrode 
// inclined in phi by angle
double DetectorGeometry::getPhi(double angle, double r_in, double r_out) {
  double L = getL(angle, r_in, r_out);
  if (r_out > r_in)
    return TMath::ASin(L/r_out * sin(angle));
  else
    return -TMath::ASin(L/r_out * sin(angle));
}

// calculate radius of point on electrode starting at 
// r=r_in, inclined in phi by angle, and at distance L 
// from beginning of electrode
double DetectorGeometry::getR(double angle, double r_in, double L) {
  return sqrt( (r_in+L*cos(angle)) * (r_in+L*cos(angle))
	       + (L*sin(angle)) * (L*sin(angle)) );
}


// calculate the derived geometry parameters of the detector
void DetectorGeometry::calcGeom(int version) {

  cout << "******************************************************************************" << endl;
  cout << "Calculating the geometry parameters" << endl;
  cout << "******************************************************************************" << endl << endl;

  cout << "ECAL" << endl << endl;

  drNom = std::vector<double>({
      15.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 
      35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm
    });
  if (version==3)
  {
    // adjust alpha
    alpha = 50.18 * TMath::Pi()/180.;
    // reduce number of layers by 1
    mergedCells_Theta.pop_back();
    mergedModules.pop_back();
    // adjust values of radial thicknesses of each layer
    drNom = std::vector<double>({
      1.65*cm, 3.36*cm, 3.46*cm, 3.56*cm, 3.68*cm, 3.79*cm,
      3.92*cm, 4.05*cm, 4.20*cm, 4.35*cm, 4.52*cm
    });
  }
  phiMin = -alpha - gridPhi/2.0;

  nLayers = drNom.size();  

  // print initial information
  cout << "r(min) = " << rMin << " " << units << endl;
  cout << "r(max) = " << rMax << " " << units << endl;
  cout << "total thickness = " << rMax-rMin << " " << units << endl;

  // calculate total length
  Ltot = getL(alpha, rMin, rMax);
  cout << "electrode length = " << Ltot << " " << units << endl << endl;
 
  // calculate length along electrode of each layer  
  cout << "n(layers) = " << nLayers << endl << endl;
  std::vector<double> L;
  double sumL(0.0);
  dL.clear();
  for (int i=0; i<nLayers; i++) {
    dL.push_back(drNom[i]/cos(alpha));
    sumL += dL[i];
  }
  cout << "electrode length/layer: " << endl;
  L.push_back(0.0);
  for (int i=0; i<nLayers; i++) {
    dL[i] = dL[i] * Ltot / sumL;
    cout <<  i  <<  " " << dL[i] << " " << units << endl;
    L.push_back(L[i] + dL[i]);
  }
  cout << endl;
  
  // calculate r and dr of each layer
  cout << "radial position of each layer: " << endl;
  r.resize(nLayers+1);
  dr.resize(nLayers);
  r[0] = rMin;
  for (int iLayer=1; iLayer<=nLayers; iLayer++) {
    r[iLayer] = getR(alpha, r[0], L[iLayer]);
    dr[iLayer-1] = r[iLayer]-r[iLayer-1];
    cout <<  iLayer-1  <<  " " << r[iLayer-1] << " - " << r[iLayer] << " " << units << endl;
  }
  cout << endl;

  // print calculated thickness of each layer
  cout << "radial thickness of each layer: " << endl;
  for (int iLayer=0; iLayer<nLayers; iLayer++) {
    cout <<  iLayer  <<  " " << dr[iLayer] << " " << units << endl;
  }
  cout << endl;

  // calculate radial position of points at middle length of
  // electrode in each layer
  cout << "r of center: " << endl;
  std::vector<double> rMed;
  for (int i=0; i<nLayers; i++) {
    rMed.push_back(getR(alpha, r[0], (L[i]+L[i+1])/2.0));
    cout <<  i  <<  " " << rMed[i] << " " << units << endl;
  }
  cout << endl;

  rAvg = getR(alpha, rMin, Ltot/2.0);
  dPhiAvg = getPhi(alpha, rMin, rAvg);
  cout << "r(avg) = " << rAvg << " " << units << endl;
  cout << "phi grid : " << gridPhi << endl;
  cout << "phi offset of electrode mid point wrt initial point : " << dPhiAvg << endl;
  cout << "initial phi of edge of module 0 : " << phiMin << endl << endl;

  
  cout << "HCAL" << endl << endl;
  
  cout << "r(min) = " << rMinHCal << " " << units << endl;
  cout << "r(max) = " << rMaxHCal << " " << units << endl;
  cout << "total thickness = " << rMaxHCal-rMinHCal << " " << units << endl << endl;
  cout << "n(layers) = " << nLayersHCal << endl << endl;
  cout << "radial position of each layer: " << endl;
  for (int iLayer=0; iLayer<nLayersHCal; iLayer++) {
    cout <<  iLayer  <<  " " << rHCal[iLayer] << " - " << rHCal[iLayer+1] << " " << units << endl;
  }
  cout << endl;
  cout << "radial thickness of each layer: " << endl;
  for (int iLayer=0; iLayer<nLayersHCal; iLayer++) {
    cout <<  iLayer  <<  " " << rHCal[iLayer+1]-rHCal[iLayer] << " " << units << endl;
  }
  cout << endl;
  cout << "z range of envelope = " << zMinHCal << " - " << zMaxHCal << " " << units << endl;
  cout << "z range of front endplate = "
       << zMaxHCal + dzHCalEndPlate - thicknessHCalEndPlate << " - "
       << zMaxHCal + dzHCalEndPlate << " " << units << endl;
  cout << "z range of back endplate = "
       << zMinHCal - dzHCalEndPlate << " - "
       << zMinHCal - dzHCalEndPlate + thicknessHCalEndPlate << " " << units << endl;
  cout << "z range of detector = "
       << zMinHCal-dzHCalEndPlate+thicknessHCalEndPlate << " - "
       << zMaxHCal+dzHCalEndPlate-thicknessHCalEndPlate << " " << units << endl;
  cout << endl;
  cout << "theta coverage of each layer: " << endl;
  for (int iLayer=0; iLayer<nLayersHCal; iLayer++) {
    float thetaMinLayer = TMath::ATan2(rHCal[iLayer], zMaxHCal+dzHCalEndPlate-thicknessHCalEndPlate);
    float thetaMaxLayer = TMath::ATan2(rHCal[iLayer], zMinHCal-dzHCalEndPlate+thicknessHCalEndPlate);
    cout <<  iLayer  <<  " " << thetaMinLayer << " - " << thetaMaxLayer << endl;
  }
  cout << endl << endl;
  
}

bool DetectorGeometry::isInECalBarrel(double rho, double z) {
  if ( (rho>=rMin) && (rho<=rMax) &&
       (z>=zMin) && (z<=zMax) )
    return true;
  else
    return false;
}

bool DetectorGeometry::isInECalEndCap(double rho, double z) {
  if ( (rho>=rMinEndCap) && (rho<=rMaxEndCap) &&
       ( (z>=zMinEndCap && z<=zMaxEndCap) ||
         (z<=-zMinEndCap && z>=-zMaxEndCap) ) )
    return true;
  else
    return false;
}


/******************************************************************************/
// HELPER FUNCTIONS related to the readout
/******************************************************************************/

// default cellID encodings, used if the event file has no CellIDEncoding
// metadata for the readout
// ECal barrel: system:0:4,cryo:4:1,type:5:3,subtype:8:3,layer:11:8,module:19:11,theta:30:10
// HCal barrel: system:0:4,layer:4:5,theta:9:9,phi:18:10
// ECal endcap: system:0:4,cryo:4:1,type:5:3,subtype:8:3,side:11:-2,wheel:13:3,layer:16:12,module:28:11,rho:39:8,z:47:8
// HCal endcap: system:0:4,type:4:3,layer:7:6,row:13:11,theta:24:11,phi:35:10
DetectorGeometry::CellIDField DetectorGeometry::fSystem = {0, 4};
DetectorGeometry::CellIDField DetectorGeometry::fECalBarrelLayer = {11, 8};
DetectorGeometry::CellIDField DetectorGeometry::fECalBarrelModule = {19, 11};
DetectorGeometry::CellIDField DetectorGeometry::fECalBarrelTheta = {30, 10};
DetectorGeometry::CellIDField DetectorGeometry::fHCalBarrelLayer = {4, 5};
DetectorGeometry::CellIDField DetectorGeometry::fHCalBarrelTheta = {9, 9};
DetectorGeometry::CellIDField DetectorGeometry::fHCalBarrelPhi = {18, 10};
DetectorGeometry::CellIDField DetectorGeometry::fECalEndCapSide = {11, -2};
DetectorGeometry::CellIDField DetectorGeometry::fECalEndCapWheel = {13, 3};
DetectorGeometry::CellIDField DetectorGeometry::fECalEndCapLayer = {16, 12};
DetectorGeometry::CellIDField DetectorGeometry::fECalEndCapModule = {28, 11};
DetectorGeometry::CellIDField DetectorGeometry::fECalEndCapRho = {39, 8};
DetectorGeometry::CellIDField DetectorGeometry::fECalEndCapZ = {47, 8};
DetectorGeometry::CellIDField DetectorGeometry::fHCalEndCapLayer = {7, 6};
DetectorGeometry::CellIDField DetectorGeometry::fHCalEndCapTheta = {24, 11};
DetectorGeometry::CellIDField DetectorGeometry::fHCalEndCapPhi = {35, 10};
DetectorGeometry::CellIDField DetectorGeometry::fMuonBarrelLayer = {5, 5};
DetectorGeometry::CellIDField DetectorGeometry::fMuonEndcapLayer = {5, 5};

bool DetectorGeometry::SetCellIDEncoding(Readout readout, const std::string& encoding)
{
  // parse the fields "name:offset:width" (or "name:width", with the offset
  // following the previous field), separated by commas
  std::map<std::string, CellIDField> fields;
  std::stringstream ss(encoding);
  std::string item;
  int nextOffset = 0;
  while (std::getline(ss, item, ','))
  {
    std::vector<std::string> parts;
    std::stringstream itemss(item);
    std::string part;
    while (std::getline(itemss, part, ':'))
      parts.push_back(part);
    try
    {
      CellIDField field;
      if (parts.size() == 3)
        field = {std::stoi(parts[1]), std::stoi(parts[2])};
      else if (parts.size() == 2)
        field = {nextOffset, std::stoi(parts[1])};
      else
        throw std::invalid_argument(item);
      fields[parts[0]] = field;
      nextOffset = field.offset + std::abs(field.width);
    }
    catch (const std::exception&)
    {
      cout << "WARNING: cannot parse field '" << item << "' of CellIDEncoding " << encoding << ", ignoring the encoding" << endl;
      return false;
    }
  }

  // the fields needed for each readout, and where to store them
  std::vector<std::pair<std::string, CellIDField*>> needed = {{"system", &fSystem}};
  std::string readoutName;
  switch (readout)
  {
  case Readout::ECalBarrel:
    readoutName = "ECal barrel";
    needed.insert(needed.end(), {{"layer", &fECalBarrelLayer}, {"module", &fECalBarrelModule}, {"theta", &fECalBarrelTheta}});
    break;
  case Readout::HCalBarrel:
    readoutName = "HCal barrel";
    needed.insert(needed.end(), {{"layer", &fHCalBarrelLayer}, {"theta", &fHCalBarrelTheta}, {"phi", &fHCalBarrelPhi}});
    break;
  case Readout::ECalEndcap:
    readoutName = "ECal endcap";
    needed.insert(needed.end(), {{"side", &fECalEndCapSide}, {"wheel", &fECalEndCapWheel}, {"layer", &fECalEndCapLayer},
                                 {"module", &fECalEndCapModule}, {"rho", &fECalEndCapRho}, {"z", &fECalEndCapZ}});
    break;
  case Readout::HCalEndcap:
    readoutName = "HCal endcap";
    needed.insert(needed.end(), {{"layer", &fHCalEndCapLayer}, {"theta", &fHCalEndCapTheta}, {"phi", &fHCalEndCapPhi}});
    break;
  case Readout::MuonBarrel:
    readoutName = "muon tagger barrel";
    needed.insert(needed.end(), {{"layer", &fMuonBarrelLayer}});
    break;
  case Readout::MuonEndcap:
    readoutName = "muon tagger endcap";
    needed.insert(needed.end(), {{"layer", &fMuonEndcapLayer}});
    break;
  }
  for (const auto& [name, field] : needed)
  {
    if (fields.find(name) == fields.end())
    {
      cout << "WARNING: field " << name << " missing in CellIDEncoding " << encoding
           << " of the " << readoutName << ", using the default encoding" << endl;
      return false;
    }
  }
  // the system field is used to find the readout of a cell, so it must be the
  // same for all the readouts: keep the first one and warn if they differ
  static bool systemSet = false;
  const CellIDField& system = fields["system"];
  if (systemSet && (system.offset != fSystem.offset || system.width != fSystem.width))
    cout << "WARNING: the system field (" << system.offset << ":" << system.width << ") in the CellIDEncoding of the "
         << readoutName << " differs from the one of the other readouts (" << fSystem.offset << ":" << fSystem.width
         << "), keeping the latter" << endl;
  else
  {
    fSystem = system;
    systemSet = true;
  }
  for (const auto& [name, field] : needed)
    if (field != &fSystem)
      *field = fields[name];
  cout << "Using CellIDEncoding " << encoding << " for the " << readoutName << endl;
  return true;
}

Long64_t DetectorGeometry::ReadField(const CellIDField& field, ULong_t cellID)
{
  const int width = std::abs(field.width);
  const ULong64_t mask = (width < 64) ? ((1ULL << width) - 1) : ~0ULL;
  const ULong64_t value = ((ULong64_t) cellID >> field.offset) & mask;
  // signed field: extend the sign
  if (field.width < 0 && ((value >> (width - 1)) & 1))
    return (Long64_t) value - (Long64_t) (1ULL << width);
  return (Long64_t) value;
}

ULong_t DetectorGeometry::ReadNbitsAtPositionMFromCellID(int n, int m, ULong_t cellID) {
  return ReadField({m, n}, cellID);
}

ULong_t DetectorGeometry::SystemID(ULong_t cellID) {
  return ReadField(fSystem, cellID);
}

ULong_t DetectorGeometry::ECalBarrelLayer(ULong_t cellID) {
  return ReadField(fECalBarrelLayer, cellID);
}

ULong_t DetectorGeometry::ECalBarrelModule(ULong_t cellID) {
  return ReadField(fECalBarrelModule, cellID);
}

ULong_t DetectorGeometry::ECalBarrelThetaBin(ULong_t cellID) {
  return ReadField(fECalBarrelTheta, cellID);
}

ULong_t DetectorGeometry::HCalBarrelLayer(ULong_t cellID) {
  return ReadField(fHCalBarrelLayer, cellID);
}

ULong_t DetectorGeometry::HCalBarrelThetaBin(ULong_t cellID) {
  return ReadField(fHCalBarrelTheta, cellID);
}

ULong_t DetectorGeometry::HCalBarrelPhiBin(ULong_t cellID) {
  return ReadField(fHCalBarrelPhi, cellID);
}

Long_t DetectorGeometry::ECalEndCapSide(ULong_t cellID) {
  return ReadField(fECalEndCapSide, cellID);
}

ULong_t DetectorGeometry::ECalEndCapWheel(ULong_t cellID) {
  return ReadField(fECalEndCapWheel, cellID);
}

ULong_t DetectorGeometry::ECalEndCapLayer(ULong_t cellID) {
  return ReadField(fECalEndCapLayer, cellID);
}

ULong_t DetectorGeometry::ECalEndCapModule(ULong_t cellID) {
  return ReadField(fECalEndCapModule, cellID);
}

ULong_t DetectorGeometry::ECalEndCapRhoBin(ULong_t cellID) {
  return ReadField(fECalEndCapRho, cellID);
}

ULong_t DetectorGeometry::ECalEndCapZBin(ULong_t cellID) {
  return ReadField(fECalEndCapZ, cellID);
}

bool DetectorGeometry::HCalEndCapCellBounds(ULong_t cellID, double &rhoMin, double &rhoMax, double &zMin, double &zMax,
                                            double &phiMin, double &phiMax) const {
  // the section and the radial range of the layer
  const int layer = (int)HCalEndCapLayer(cellID);
  int first = 0, section = -1;
  for (size_t s = 0; s < numLayersHCalEndCap.size() && section < 0; s++) {
    double r = rMinHCalEndCapSections[s];
    for (size_t i = 0; i < dRLayersHCalEndCap.size() && section < 0; i++) {
      for (int k = 0; k < numLayersHCalEndCap[s][i]; k++, first++) {
        if (first == layer) {
          section = s;
          rhoMin = r;
          rhoMax = r + dRLayersHCalEndCap[i];
          break;
        }
        r += dRLayersHCalEndCap[i];
      }
    }
  }
  if (section < 0)
    return false;

  // the range in z of the theta bin: half way to the neighbouring bins in the
  // layer (at its centre), or the limit of the section (on the side of the cell)
  const double radius = 0.5 * (rhoMin + rhoMax);
  const int thetaBin = (int)HCalEndCapThetaBin(cellID);
  auto binToZ = [&](int bin) { return radius / std::tan(thetaOffsetHCalEndCap + bin * thetaGridHCalEndCap); };
  const double z = binToZ(thetaBin);
  const double sign = z > 0 ? 1. : -1.;
  const double low = sign * zMinHCalEndCapSections[section], high = sign * zMaxHCalEndCapSections[section];
  if (std::fabs(z) < zMinHCalEndCapSections[section] || std::fabs(z) > zMaxHCalEndCapSections[section])
    return false;
  auto edge = [&](int neighbour) {
    const double zn = binToZ(neighbour);
    if (std::fabs(zn) <= zMinHCalEndCapSections[section])
      return low;
    if (std::fabs(zn) >= zMaxHCalEndCapSections[section])
      return high;
    return 0.5 * (z + zn);
  };
  const double edge1 = edge(thetaBin - 1), edge2 = edge(thetaBin + 1);
  zMin = std::min(edge1, edge2);
  zMax = std::max(edge1, edge2);

  // the phi bin, centred at -pi + (bin + 1/2) * 2pi/nbins
  const double gridPhi = TMath::TwoPi() / nPhiBinsHCalEndCap;
  phiMin = -TMath::Pi() + HCalEndCapPhiBin(cellID) * gridPhi;
  phiMax = phiMin + gridPhi;
  return true;
}

bool DetectorGeometry::ECalEndCapCellRhoZ(ULong_t cellID, double &rho, double &z, double &dRho, double &dZ) const {
  ULong_t wheel = ECalEndCapWheel(cellID);
  if (wheel >= (ULong_t)nWheelsECalEndCap)
    return false;
  ULong_t rhoBin = ECalEndCapRhoBin(cellID);
  ULong_t zBin = ECalEndCapZBin(cellID);
  if (rhoBin >= (ULong_t)nRhoBinsECalEndCap[wheel] || zBin >= (ULong_t)nZBinsECalEndCap[wheel])
    return false;
  dRho = gridSizeRhoECalEndCap[wheel];
  dZ = gridSizeZECalEndCap[wheel];
  rho = offsetRhoECalEndCap[wheel] + (rhoBin + 0.5) * dRho;
  z = (ECalEndCapSide(cellID) > 0 ? 1. : -1.) * (offsetZECalEndCap[wheel] + (zBin + 0.5) * dZ);
  return true;
}

ULong_t DetectorGeometry::HCalEndCapLayer(ULong_t cellID) {
  return ReadField(fHCalEndCapLayer, cellID);
}

ULong_t DetectorGeometry::MuonBarrelLayer(ULong_t cellID) {
  return ReadField(fMuonBarrelLayer, cellID);
}

ULong_t DetectorGeometry::MuonEndcapLayer(ULong_t cellID) {
  return ReadField(fMuonEndcapLayer, cellID);
}

ULong_t DetectorGeometry::HCalEndCapThetaBin(ULong_t cellID) {
  return ReadField(fHCalEndCapTheta, cellID);
}

ULong_t DetectorGeometry::HCalEndCapPhiBin(ULong_t cellID) {
  return ReadField(fHCalEndCapPhi, cellID);
}
