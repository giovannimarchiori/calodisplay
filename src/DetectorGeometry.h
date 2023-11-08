/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
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
  
  const std::vector<int> mergedCells_Theta = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  const std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  
  // units
  // G4 uses mm but Root uses cm..
  const double cm = 1.0;
  const double mm = 0.1;
  const std::string units = "cm";
  
  // z extension of barrel
  const double zMin = -3100.*mm;
  const double zMax =  3100.*mm;
  
  // radial extension of barrel
  const double rMin = 217.28*cm;
  const double rMax = 257.33*cm;

  // nominal radial thickness of layers
  std::vector<double> drNom;
  
  // number of layers
  int nLayers;
  
  // number of electrodes ECAL
  const int nModules = 1536;
  
  // inclination angle of electrodes
  const double alpha = 50*TMath::Pi()/180.;
  
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

  const double etaGridHCal = 0.025;
  const double etaOffsetHCal = -0.9;
  const int nPhiBinsHCal = 256;
  
  // derived
  const int nLayersHCal = rHCal.size()-1;
  const double rMinHCal = rHCal[0];
  const double rMaxHCal = rHCal[nLayersHCal];
  const double etaMinHCal = etaOffsetHCal - etaGridHCal/2.0;
  const double etaMaxHCal = -etaOffsetHCal + etaGridHCal/2.0;
  const int nEtaBinsHCal = (etaMaxHCal - etaMinHCal)/etaGridHCal;
  const double gridPhiHCal = TMath::TwoPi()/nPhiBinsHCal;
  const double phiMinHCal = -TMath::Pi();

  
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
  const double phiMin = -alpha - gridPhi/2.0;

  // other quantities, calculated by calcGeom
  // radial position of each layer
  std::vector<double> r;
  std::vector<double> dr;
  // length of electrode along each layer (dRnom/cos(alpha))
  std::vector<double> dL;

  // constructor
  DetectorGeometry();
  
  // calculate derived parameters of geometry depending on the main ones
  // (such as radial/electrode length of each layer)
  // and print them to screen
  void calcGeom();


  /******************************************************************************/
  // HELPER FUNCTIONS related to the readout
  /******************************************************************************/

  // extract number encoded in N bits starting at bit M from cellID  
  static ULong_t ReadNbitsAtPositionMFromCellID(int n, int m, ULong_t cellID);

  // extract system number from cellID
  static ULong_t SystemID(ULong_t cellID);

  // extract layer number from cellID
  static ULong_t Layer(ULong_t cellID);
  
  // extract module number from cellID
  static ULong_t Module(ULong_t cellID);
  
  // extract theta bin from cellID
  static ULong_t ThetaBin(ULong_t cellID);

  //ClassDef(DetectorGeometry, 0)
};

#endif
