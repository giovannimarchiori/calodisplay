/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class CaloCluster: contains cluster summary information
//
/******************************************************************************/

#ifndef CALOCLUSTER_H
#define CALOCLUSTER_H

/******************************************************************************/
// dependencies
/******************************************************************************/
#include <array>
#include <vector>
#include "TVector3.h"

/******************************************************************************/
// class declaration
/******************************************************************************/
class CaloCluster
{
public:
  // the parts of the calorimeters, with the energy and barycenter of the
  // cluster in each of their layers
  enum Part { kECalBarrel, kECalEndcap, kHCalBarrel, kHCalEndcap, kMuon, kNParts };

  // constructor
  CaloCluster() = default;

  unsigned int getIndex() const { return m_index; };
  void setIndex(unsigned int index) { m_index = index; };

  float getEnergy() const { return m_energy; };
  void setEnergy(float energy) { m_energy = energy; };

  // energy in a part, and in the ECal, HCal and muon tagger (summed over the layers)
  float getEnergy(Part part) const;
  float getEnergyECal() const { return getEnergy(kECalBarrel) + getEnergy(kECalEndcap); };
  float getEnergyHCal() const { return getEnergy(kHCalBarrel) + getEnergy(kHCalEndcap); };
  float getEnergyMuon() const { return getEnergy(kMuon); };

  TVector3 getBarycenter() const { return m_barycenter; };
  void setBarycenterXYZ(float x, float y, float z) { m_barycenter = TVector3(x, y, z); };

  // energy and barycenter in each layer of a part
  unsigned int getNLayers(Part part) const { return m_layers[part].energy.size(); };
  // 0 if the layer is out of range
  float getEnergyInLayer(Part part, unsigned int layer) const;
  // the origin if the layer is out of range
  TVector3 getBarycenterInLayer(Part part, unsigned int layer) const;
  void setEnergyVsLayers(Part part, const std::vector<float> &energyVsLayer) { m_layers[part].energy = energyVsLayer; };
  void setBarycenterVsLayers(Part part, const std::vector<TVector3> &barycenterVsLayer) { m_layers[part].barycenter = barycenterVsLayer; };

  void print() const;

private:
  unsigned int m_index = 0;
  float m_energy = 0.;
  TVector3 m_barycenter;

  struct Layers
  {
    std::vector<float> energy;
    std::vector<TVector3> barycenter;
  };
  std::array<Layers, kNParts> m_layers;

  // direction ...
  // shower shapes...
};

#endif
