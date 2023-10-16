/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class DetectorGeometry: contains detector geometry information
//
/******************************************************************************/

#include "DetectorGeometry.h"

#include <iostream>
using std::cout;
using std::endl;

DetectorGeometry::DetectorGeometry() {
  calcGeom();
}


/******************************************************************************/
// GEOMETRY HELPER FUNCTIONS AND DERIVED GEOMETRY PARAMETERS
/******************************************************************************/

// calculate length along electrode at r=r_out for electrode
// starting at r=r_in and inclined in phi by alpha
// for r_in>r_out, when the solution exists, it's always > 0
double DetectorGeometry::getL(double alpha, double r_in, double r_out) {
  if (r_out > r_in)
    return sqrt(r_out*r_out - r_in*r_in*sin(alpha)*sin(alpha)) - r_in*cos(alpha);
  else
    return r_in*cos(alpha) - sqrt(r_out*r_out - r_in*r_in*sin(alpha)*sin(alpha));
}

// calculate phi offset between point at r=r_out with respect 
// to point at r=r_in if they are on the same electrode 
// inclined in phi by alpha
double DetectorGeometry::getPhi(double alpha, double r_in, double r_out) {
  double L = getL(alpha, r_in, r_out);
  if (r_out > r_in)
    return TMath::ASin(L/r_out * sin(alpha));
  else
    return -TMath::ASin(L/r_out * sin(alpha));
}

// calculate radius of point on electrode starting at 
// r=r_in, inclined in phi by alpha, and at distance L 
// from beginning of electrode
double DetectorGeometry::getR(double alpha, double r_in, double L) {
  return sqrt( (r_in+L*cos(alpha)) * (r_in+L*cos(alpha))
	       + (L*sin(alpha)) * (L*sin(alpha)) );
}


// calculate the derived geometry parameters of the detector
void DetectorGeometry::calcGeom() {

  cout << "******************************************************************************" << endl;
  cout << "Calculating the geometry parameters" << endl;
  cout << "******************************************************************************" << endl << endl;

  drNom = std::vector<double>({
      15.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 
      35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm, 35.*mm
    });
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
  cout <<  "0 " << r[0] << " " << units << endl;
  for (int iLayer=1; iLayer<=nLayers; iLayer++) {
    r[iLayer] = getR(alpha, r[0], L[iLayer]);
    dr[iLayer-1] = r[iLayer]-r[iLayer-1];
    cout <<  iLayer  <<  " " << r[iLayer] << " " << units << endl;
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
}


/******************************************************************************/
// HELPER FUNCTIONS related to the readout
/******************************************************************************/

// extract layer number from cellID
ULong_t DetectorGeometry::Layer(ULong_t cellID) {
  const ULong_t mask = (1<<8) -1;
  return (cellID >> 11) & mask;
}

// extract module number from cellID
ULong_t DetectorGeometry::Module(ULong_t cellID) {
  const ULong_t mask = (1<<11) -1;
  return (cellID >> 19) & mask;
}

// extract theta bin from cellID
ULong_t DetectorGeometry::ThetaBin(ULong_t cellID) {
  const ULong_t mask = (1<<10) -1;
  return (cellID >> 30) & mask;
}
