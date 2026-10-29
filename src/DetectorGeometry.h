/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class DetectorGeometry: contains detector geometry information
//
/******************************************************************************/


#ifndef DETECTORGEOMETRY_H
#define DETECTORGEOMETRY_H

/******************************************************************************/
// dependencies
/******************************************************************************/
#include <TMath.h>
#include <string>
#include <vector>
#include "Units.h"

class DetectorGeometry {
public:

  /******************************************************************************/
  // GEOMETRICAL CONSTANTS
  // They should match the geometry used in the geometry file and the readout
  // used to produce the event file
  /******************************************************************************/
  
  // magnetic field (only in barrel: outside it is assumed to be zero)
  // zmin, zmax = ECAL barrel zmin, zmax (EMBarrel_dz = BarECal_dz-CryoBarrelSide)
  // Rsol = BarCryoECal_rmax-CryoBarrelBackWarm-CryoBarrelBackCold
  const double Bin = 2.;  // in Tesla, magnetic field for zmin<z<zmax, R<Rsol (default, can be changed with bField)
  const double Bout = -1.;  // in Tesla, magnetic field for zmin<z<zmax, R>Rsol

  // ECAL barrel

  // z extension
  const double zMin = -3100.*Units::mm;
  const double zMax =  3100.*Units::mm;
  
  // radial extension
  const double rMin = 217.28*Units::cm;
  const double rMax = 257.83*Units::cm;

  // nominal radial thickness of layers in barrel
  std::vector<double> drNom;
  
  // number of layers in barrel
  int nLayers;
  
  // number of electrodes ECAL barrel
  const int nModules = 1536;
  
  // inclination angle of electrodes in ECAL barrel
  double alpha = 50*TMath::Pi()/180.;
  
  // grid in theta
  // - size of cell
  const double thetaGrid = 0.00981748/4;

  // - theta of edge of first cell
  const double thetaMin = 0.59027850 - thetaGrid/2.;

  // - n(bins)
  const int nThetaBins = 800;

  // merged cells in theta and modules direction
  std::vector<int> mergedCells_Theta = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  

  // ECAL endcap

  // z extent
  const double rMinEndCap = 350.0*Units::mm;
  const double rMaxEndCap = 2900.0*Units::mm;
  const double zMinEndCap = 3200.*Units::mm;
  const double zMaxEndCap = 3850.*Units::mm;

  // number of layers
  const int nLayersECalEndCap = 10;

  // readout segmentation (FCCSWEndcapTurbine_k4geo): in each wheel, the cells
  // are bins of size gridSizeRho in rho starting from offsetRho, and of size
  // gridSizeZ in |z| starting from offsetZ (from ECalEndcaps_Turbine_o1_v03.xml:
  // offset_rho = EMECRminWheel<i>+EMECSupportTubeThickness, grid_size_rho =
  // EMECGridSizeRho<i>, offset_z = EMEC_z1, grid_size_z = EMECGridSizeZ<i>)
  const int nWheelsECalEndCap = 3;
  const std::vector<int> nRhoBinsECalEndCap = {10, 14, 34};
  const std::vector<int> nZBinsECalEndCap = {10, 10, 10};
  const std::vector<double> offsetRhoECalEndCap = {410.0*Units::mm, 774.2605*Units::mm, 1470.2352*Units::mm};
  const std::vector<double> gridSizeRhoECalEndCap = {35.4260*Units::mm, 48.9982*Units::mm, 38.8166*Units::mm};
  const std::vector<double> offsetZECalEndCap = {3260.0*Units::mm, 3260.0*Units::mm, 3260.0*Units::mm};
  const std::vector<double> gridSizeZECalEndCap = {45.0*Units::mm, 45.0*Units::mm, 45.0*Units::mm};


  // HCAL barrel

  // radial positions of layer boundaries
  const std::vector<double> rHCal = {
    281.05*Units::cm, 286.05*Units::cm, 291.05*Units::cm, 296.05*Units::cm, 301.05*Units::cm,
    311.05*Units::cm, 321.05*Units::cm, 331.05*Units::cm, 341.05*Units::cm, 351.05*Units::cm,
    361.05*Units::cm, 381.05*Units::cm, 401.05*Units::cm, 421.05*Units::cm};

  // z: envelope up to +-280*cm
  //    end plates between +-279.05 and +-279.55*cm
  //    detector between -279.05cm and 279.05cm
  const double zMaxHCal =  280.0*Units::cm;
  const double zMinHCal = -280.0*Units::cm;
  const double dzHCalEndPlate = -0.45*Units::cm; // distance of end plate from envelope
  const double thicknessHCalEndPlate = 0.5*Units::cm; // thickness of endplate

  // might have to shift offset by half grid size.. will check
  const double thetaGridHCal = 0.02218;
  const double thetaOffsetHCal = 0.783406;
  const int nThetaBinsHCal = 72;
  const int nPhiBinsHCal = 256;

  // derived
  const int nLayersHCal = rHCal.size()-1;
  const double rMinHCal = rHCal[0];
  const double rMaxHCal = rHCal[nLayersHCal];
  const double thetaMinHCal = thetaOffsetHCal - thetaGridHCal/2.0;
  const double thetaMaxHCal = thetaOffsetHCal + nThetaBinsHCal*thetaGridHCal;
  const double gridPhiHCal = TMath::TwoPi()/nPhiBinsHCal;
  const double phiMinHCal = -TMath::Pi();

  // HCAL endcap

  // readout segmentation (FCCSWHCalPhiTheta_k4geo, from
  // HCalEndcaps_ThreeParts_TileCal_v03.xml): three sections in z, each with
  // layers in rho, numbered consecutively over the sections (in each section,
  // numLayersHCalEndCap[section][i] layers of thickness dRLayersHCalEndCap[i],
  // starting from rMinHCalEndCapSections), and bins in theta and phi. In each
  // layer, the cell of a theta bin extends in z around r/tan(theta) up to half
  // way to the neighbouring bins, within the section
  const int nLayersHCalEndCap = 37;
  const std::vector<double> zMinHCalEndCapSections = {2900.*Units::mm, 3400.*Units::mm, 3900.*Units::mm};
  const std::vector<double> zMaxHCalEndCapSections = {3400.*Units::mm, 3900.*Units::mm, 5450.*Units::mm};
  // (inner radius of the section + face plate thickness + plate space)
  const std::vector<double> rMinHCalEndCapSections = {3610.5*Units::mm, 2910.5*Units::mm, 360.5*Units::mm};
  const std::vector<std::vector<int>> numLayersHCalEndCap = {{5, 1, 0}, {4, 3, 2}, {4, 10, 8}};
  const std::vector<double> dRLayersHCalEndCap = {100.*Units::mm, 150.*Units::mm, 250.*Units::mm};
  const double thetaGridHCalEndCap = 0.022180;
  const double thetaOffsetHCalEndCap = 0.007106;
  const int nPhiBinsHCalEndCap = 256;

  // muon tagger: radial layers in the barrel, layers in z in the endcaps
  // (MuonTaggerBarrelLayers and MuonTaggerEndcapLayers in MuonTaggerPhiTheta.xml)
  const int nLayersMuonBarrel = 2;
  const int nLayersMuonEndcap = 2;
  // extent of the barrel and of the endcaps (DectDimensions.xml)
  const double rMinMuonBarrel = 4550.*Units::mm;
  const double rMaxMuonBarrel = 5000.*Units::mm;
  const double zMaxMuonBarrel = 6000.*Units::mm;
  const double rMinMuonEndcap = 1100.*Units::mm;
  const double rMaxMuonEndcap = 4500.*Units::mm;
  const double zMinMuonEndcap = 5550.*Units::mm;
  const double zMaxMuonEndcap = 6000.*Units::mm;
  // readout segmentation, the same in the barrel and the endcaps
  // (FCCSWGridPhiTheta_k4geo): bins in theta and phi, centred at
  // offset + i * grid size
  const double thetaGridMuon = 0.0087266;
  const double thetaOffsetMuon = 0.1047277;
  const int nPhiBinsMuon = 704;
  const double phiOffsetMuon = -TMath::Pi() + TMath::Pi() / nPhiBinsMuon;
  
  /******************************************************************************/
  // GEOMETRY HELPER FUNCTIONS AND DERIVED GEOMETRY PARAMETERS
  /******************************************************************************/
  
  
  // calculate length along electrode at r=r_out for electrode
  // starting at r=r_in and inclined in phi by angle
  // for r_in>r_out, when the solution exists, it's always > 0
  double getL(double angle, double r_in, double r_out);

  // calculate phi offset between point at r=r_out with respect 
  // to point at r=r_in if they are on the same electrode 
  // inclined in phi by angle
  double getPhi(double angle, double r_in, double r_out);

  // calculate radius of point on electrode starting at 
  // r=r_in, inclined in phi by angle, and at distance L 
  // from beginning of electrode
  double getR(double angle, double r_in, double L);

  // length of electrodes
  double Ltot;

  // phi separation of modules
  const double gridPhi = TMath::TwoPi()/nModules;

  // r for electrode length = L/2
  double rAvg;

  // delta phi of point at L/2
  double dPhiAvg;

  // phi edge of module 0
  double phiMin;

  // other quantities, calculated by calcGeom
  // radial position of each layer
  std::vector<double> r;
  std::vector<double> dr;
  // length of electrode along each layer (dRnom/cos(alpha))
  std::vector<double> dL;

  // constructor
  DetectorGeometry(int version=2);
  
  // calculate derived parameters of geometry depending on the main ones
  // (such as radial/electrode length of each layer)
  // and print them to screen
  void calcGeom(int version);

  bool isInECalBarrel(double rho, double z);
  bool isInECalEndCap(double rho, double z);

  // the centre in rho and z (signed) and the size of an ECAL endcap cell in the
  // readout segmentation; returns false if its wheel or bins are out of range
  bool ECalEndCapCellRhoZ(ULong_t cellID, double &rho, double &z, double &dRho, double &dZ) const;

  // the boundaries in rho, z (signed) and phi of an HCAL endcap cell in the
  // readout segmentation; returns false if its layer or theta bin are out of
  // range
  bool HCalEndCapCellBounds(ULong_t cellID, double &rhoMin, double &rhoMax, double &zMin, double &zMax,
                            double &phiMin, double &phiMax) const;

  /******************************************************************************/
  // HELPER FUNCTIONS related to the readout
  /******************************************************************************/

  // a field of the cellID: bit offset and width (negative width: signed
  // field), as in the DD4hep CellIDEncoding strings
  struct CellIDField
  {
    int offset;
    int width;
  };

  // the values of the system field of the cellID for each subdetector
  enum System : ULong_t
  {
    kECalBarrel = 4,
    kECalEndcap = 5,
    kHCalBarrel = 8,
    kHCalEndcap = 9,
    kMuonBarrel = 12,
    kMuonEndcap = 13
  };

  // the calorimeter readouts for which the cellID is decoded
  enum class Readout { ECalBarrel, ECalEndcap, HCalBarrel, HCalEndcap, MuonBarrel, MuonEndcap };

  // set the cellID fields of a calorimeter readout from its CellIDEncoding
  // string (e.g. "system:0:4,layer:4:5,..."), as stored in the metadata of
  // the event file. Until this is called, default encodings are used.
  // Returns false (and keeps the current fields) if a needed field is missing
  static bool SetCellIDEncoding(Readout readout, const std::string& encoding);

  // extract the value of a field from cellID
  static Long64_t ReadField(const CellIDField& field, ULong_t cellID);

  // extract number encoded in N bits starting at bit M from cellID
  static ULong_t ReadNbitsAtPositionMFromCellID(int n, int m, ULong_t cellID);

  // extract system number from cellID
  static ULong_t SystemID(ULong_t cellID);

  // ECal barrel: extract layer, module, theta bin from cellID
  static ULong_t ECalBarrelLayer(ULong_t cellID);
  static ULong_t ECalBarrelModule(ULong_t cellID);
  static ULong_t ECalBarrelThetaBin(ULong_t cellID);

  // HCal barrel: extract layer, theta, phi bin from cellID
  static ULong_t HCalBarrelLayer(ULong_t cellID);
  static ULong_t HCalBarrelThetaBin(ULong_t cellID);
  static ULong_t HCalBarrelPhiBin(ULong_t cellID);

  // ECal endcap: extract side, wheel, layer, module, rho, z
  static Long_t ECalEndCapSide(ULong_t cellID);
  static ULong_t ECalEndCapWheel(ULong_t cellID);
  static ULong_t ECalEndCapLayer(ULong_t cellID);
  static ULong_t ECalEndCapModule(ULong_t cellID);
  static ULong_t ECalEndCapRhoBin(ULong_t cellID);
  static ULong_t ECalEndCapZBin(ULong_t cellID);

  // HCal endcap: extract layer, theta, phi bin from cellID
  static ULong_t HCalEndCapLayer(ULong_t cellID);
  static ULong_t HCalEndCapThetaBin(ULong_t cellID);
  static ULong_t HCalEndCapPhiBin(ULong_t cellID);

  // muon tagger: extract layer from cellID
  static ULong_t MuonBarrelLayer(ULong_t cellID);
  static ULong_t MuonEndcapLayer(ULong_t cellID);

private:
  // the cellID fields used by the functions above
  static CellIDField fSystem;
  static CellIDField fECalBarrelLayer, fECalBarrelModule, fECalBarrelTheta;
  static CellIDField fHCalBarrelLayer, fHCalBarrelTheta, fHCalBarrelPhi;
  static CellIDField fECalEndCapSide, fECalEndCapWheel, fECalEndCapLayer,
                     fECalEndCapModule, fECalEndCapRho, fECalEndCapZ;
  static CellIDField fHCalEndCapLayer, fHCalEndCapTheta, fHCalEndCapPhi;
  static CellIDField fMuonBarrelLayer, fMuonEndcapLayer;

  //ClassDef(DetectorGeometry, 0)
};

#endif
