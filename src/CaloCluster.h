/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
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

  TVector3 getBarycenter() const { return m_barycenter; };
  void setBarycenterXYZ(float x, float y, float z) { m_barycenter = TVector3(x, y, z); };

  unsigned int getNLayersECal() const { return m_energyVsECalLayer.size(); };
  float getEnergyInECalLayer(unsigned int layer) const;
  void setEnergyVsECalLayers(const std::vector<float> &energyVsECalLayer);
  TVector3 getBarycenterInECalLayer(unsigned int layer) const;
  void setBarycenterVsECalLayers(const std::vector<TVector3> &barycenterVsECalLayer);

  unsigned int getNLayersHCal() const { return m_energyVsHCalLayer.size(); };
  float getEnergyInHCalLayer(unsigned int layer) const;
  void setEnergyVsHCalLayers(const std::vector<float> &energyVsHCalLayer);
  TVector3 getBarycenterInHCalLayer(unsigned int layer) const;
  void setBarycenterVsHCalLayers(const std::vector<TVector3> &barycenterVsHCalLayer);

  void print() const;

private:
  unsigned int m_index;
  float m_energy;
  TVector3 m_barycenter;
  std::vector<float> m_energyVsECalLayer;
  std::vector<TVector3> m_barycenterVsECalLayer;
  std::vector<float> m_energyVsHCalLayer;
  std::vector<TVector3> m_barycenterVsHCalLayer;
  // direction ...
  // shower shapes...
};

#endif