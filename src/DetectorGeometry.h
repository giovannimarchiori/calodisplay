/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
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

class DetectorGeometry {
public:

  /******************************************************************************/
  // GEOMETRICAL CONSTANTS
  // They should match the geometry used in the geometry file and the readout
  // used to produce the event file
  /******************************************************************************/
  
  std::vector<int> mergedCells_Theta = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  
  // units
  // G4 uses mm but Root uses cm..
  const double cm = 1.0;
  const double mm = 0.1;
  const std::string units = "cm";

  // magnetic field (only in barrel: outside it is assumed to be zero)
  // zmin, zmax = ECAL barrel zmin, zmax (EMBarrel_dz = BarECal_dz-CryoBarrelSide)
  // Rsol = BarCryoECal_rmax-CryoBarrelBackWarm-CryoBarrelBackCold
  const double Bin = -2.;  // in Tesla, magnetic field for zmin<z<zmax, R<Rsol
  const double Bout = 1.;  // in Tesla, magnetic field for zmin<z<zmax, R>Rsol

  
  // z extension of barrel
  const double zMin = -3100.*mm;
  const double zMax =  3100.*mm;
  
  // radial extension of barrel
  const double rMin = 217.28*cm;
  const double rMax = 257.83*cm;

  // endcap
  const double zMinEndCap = 3200.*mm;
  const double zMaxEndCap = 3850.*mm;
  
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
  

  // HCAL
  const std::vector<double> rHCal = {
    281.05*cm, 286.05*cm, 291.05*cm, 296.05*cm, 301.05*cm,
    311.05*cm, 321.05*cm, 331.05*cm, 341.05*cm, 351.05*cm,
    361.05*cm, 381.05*cm, 401.05*cm, 421.05*cm};

  // z: envelope up to +-280*cm
  //    end plates between +-279.05 and +-279.55*cm
  //    detector between -279.05cm and 279.05cm
  const double zMaxHCal =  280.0*cm;
  const double zMinHCal = -280.0*cm;
  const double dzHCalEndPlate = -0.45*cm; // distance of end plate from envelope
  const double thicknessHCalEndPlate = 0.5*cm; // thickness of endplate

  //const double etaGridHCal = 0.025;
  //const double etaOffsetHCal = -0.9;
  // might have to shift offset by half grid size.. will check
  const double thetaGridHCal = 0.02218;
  //const double thetaOffsetHCal = 0.772316;
  const double thetaOffsetHCal = 0.783406;
  const int nThetaBinsHCal = 72;
  const int nPhiBinsHCal = 256;

  // derived
  const int nLayersHCal = rHCal.size()-1;
  const double rMinHCal = rHCal[0];
  const double rMaxHCal = rHCal[nLayersHCal];
  //const double etaMinHCal = etaOffsetHCal - etaGridHCal/2.0;
  //const double etaMaxHCal = -etaOffsetHCal + etaGridHCal/2.0;
  //const int nEtaBinsHCal = (etaMaxHCal - etaMinHCal)/etaGridHCal;
  const double thetaMinHCal = thetaOffsetHCal - thetaGridHCal/2.0;
  const double thetaMaxHCal = thetaOffsetHCal + nThetaBinsHCal*thetaGridHCal;
  const double gridPhiHCal = TMath::TwoPi()/nPhiBinsHCal;
  const double phiMinHCal = -TMath::Pi();

  // muon
  const int nLayersMuon = 2;
  
  /******************************************************************************/
  // GEOMETRY HELPER FUNCTIONS AND DERIVED GEOMETRY PARAMETERS
  /******************************************************************************/
  
  
  // calculate length along electrode at r=r_out for electrode
  // starting at r=r_in and inclined in phi by alpha
  // for r_in>r_out, when the solution exists, it's always > 0
  double getL(double alpha, double r_in, double r_out);

  // calculate phi offset between point at r=r_out with respect 
  // to point at r=r_in if they are on the same electrode 
  // inclined in phi by alpha
  double getPhi(double alpha, double r_in, double r_out);

  // calculate radius of point on electrode starting at 
  // r=r_in, inclined in phi by alpha, and at distance L 
  // from beginning of electrode
  double getR(double alpha, double r_in, double L);

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


  /******************************************************************************/
  // HELPER FUNCTIONS related to the readout
  /******************************************************************************/

  // extract number encoded in N bits starting at bit M from cellID  
  static ULong_t ReadNbitsAtPositionMFromCellID(int n, int m, ULong_t cellID);

  // extract system number from cellID
  static ULong_t SystemID(ULong_t cellID);

  // extract layer number from cellID
  static ULong_t ECalBarrelLayer(ULong_t cellID);
  
  // extract module number from cellID
  static ULong_t ECalBarrelModule(ULong_t cellID);
  
  // extract theta bin from cellID
  static ULong_t ECalBarrelThetaBin(ULong_t cellID);

  // extract layer number from cellID
  static ULong_t HCalBarrelLayer(ULong_t cellID);

  // extract row number from cellID
  // static ULong_t HCalBarrelRow(ULong_t cellID);

  // extract theta bin from cellID
  static ULong_t HCalBarrelThetaBin(ULong_t cellID);

  // extract phi bin from cellID
  static ULong_t HCalBarrelPhiBin(ULong_t cellID);
 
  //ClassDef(DetectorGeometry, 0)
};

#endif
