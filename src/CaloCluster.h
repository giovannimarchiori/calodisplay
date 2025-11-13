/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
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
#include <vector>
#include "TVector3.h"

/******************************************************************************/
// class declaration
/******************************************************************************/
class CaloCluster
{
public:
  // constructor
  CaloCluster() = default;

  unsigned int getIndex() const { return m_index; };
  void setIndex(unsigned int index) { m_index = index; };

  float getEnergy() const { return m_energy; };
  void setEnergy(float energy) { m_energy = energy; };

  float getEnergyECal() const;
  float getEnergyHCal() const;
  float getEnergyMuon() const;

  TVector3 getBarycenter() const { return m_barycenter; };
  void setBarycenterXYZ(float x, float y, float z) { m_barycenter = TVector3(x, y, z); };

  unsigned int getNLayersECalBarrel() const { return m_energyVsECalBarrelLayer.size(); };
  float getEnergyInECalBarrelLayer(unsigned int layer) const;
  void setEnergyVsECalBarrelLayers(const std::vector<float> &energyVsECalBarrelLayer);
  TVector3 getBarycenterInECalBarrelLayer(unsigned int layer) const;
  void setBarycenterVsECalBarrelLayers(const std::vector<TVector3> &barycenterVsECalBarrelLayer);

  unsigned int getNLayersECalEndCap() const { return m_energyVsECalEndCapLayer.size(); };
  float getEnergyInECalEndCapLayer(unsigned int layer) const;
  void setEnergyVsECalEndCapLayers(const std::vector<float> &energyVsECalEndCapLayer);
  TVector3 getBarycenterInECalEndCapLayer(unsigned int layer) const;
  void setBarycenterVsECalEndCapLayers(const std::vector<TVector3> &barycenterVsECalEndCapLayer);
  
  unsigned int getNLayersHCal() const { return m_energyVsHCalLayer.size(); };
  float getEnergyInHCalLayer(unsigned int layer) const;
  void setEnergyVsHCalLayers(const std::vector<float> &energyVsHCalLayer);
  TVector3 getBarycenterInHCalLayer(unsigned int layer) const;
  void setBarycenterVsHCalLayers(const std::vector<TVector3> &barycenterVsHCalLayer);

  unsigned int getNLayersMuon() const { return m_energyVsMuonLayer.size(); };
  float getEnergyInMuonLayer(unsigned int layer) const;
  void setEnergyVsMuonLayers(const std::vector<float> &energyVsMuonLayer);
  TVector3 getBarycenterInMuonLayer(unsigned int layer) const;
  void setBarycenterVsMuonLayers(const std::vector<TVector3> &barycenterVsMuonLayer);

  void print() const;

private:
  unsigned int m_index;
  float m_energy;
  TVector3 m_barycenter;
  std::vector<float> m_energyVsECalBarrelLayer;
  std::vector<TVector3> m_barycenterVsECalBarrelLayer;
  std::vector<float> m_energyVsECalEndCapLayer;
  std::vector<TVector3> m_barycenterVsECalEndCapLayer;
  std::vector<float> m_energyVsHCalLayer;
  std::vector<TVector3> m_barycenterVsHCalLayer;
  std::vector<float> m_energyVsMuonLayer;
  std::vector<TVector3> m_barycenterVsMuonLayer;
  // direction ...
  // shower shapes...
};

#endif
