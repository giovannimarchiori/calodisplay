/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class CaloCluster: contains cluster summary information
//
/******************************************************************************/

/******************************************************************************/
// dependencies
/******************************************************************************/
#include "CaloCluster.h"

#include <iostream>
using std::cout;
using std::endl;

/******************************************************************************/
// class definition
/******************************************************************************/

float CaloCluster::getEnergyECal() const
{
  float E(0.);
  for (float Elayer : m_energyVsECalBarrelLayer)
  {
    E += Elayer;
  }
  for (float Elayer : m_energyVsECalEndCapLayer)
  {
    E += Elayer;
  }
  return E;
}

float CaloCluster::getEnergyHCal() const
{
  float E(0.);
  for (float Elayer : m_energyVsHCalLayer)
  {
    E += Elayer;
  }
  return E;
}

float CaloCluster::getEnergyMuon() const
 {
  float E(0.);
  for (float Elayer : m_energyVsMuonLayer)
  {
    E += Elayer;
  }
  return E;
}

float CaloCluster::getEnergyInECalBarrelLayer(unsigned int layer) const
{
  if (layer > m_energyVsECalBarrelLayer.size())
  {
    cout << "CaloCluster::getEnergyInECalBarrelLayer : layer outside of range, returning 0.0" << endl;
    return 0;
  }
  else
  {
    return m_energyVsECalBarrelLayer[layer];
  }
}

float CaloCluster::getEnergyInECalEndCapLayer(unsigned int layer) const
{
  if (layer > m_energyVsECalEndCapLayer.size())
  {
    cout << "CaloCluster::getEnergyInECalEndCapLayer : layer outside of range, returning 0.0" << endl;
    return 0;
  }
  else
  {
    return m_energyVsECalEndCapLayer[layer];
  }
}

float CaloCluster::getEnergyInHCalLayer(unsigned int layer) const
{
  if (layer > m_energyVsHCalLayer.size())
  {
    cout << "CaloCluster::getEnergyInHCalLayer : layer outside of range, returning 0.0" << endl;
    return 0;
  }
  else
  {
    return m_energyVsHCalLayer[layer];
  }
}

float CaloCluster::getEnergyInMuonLayer(unsigned int layer) const
{
  if (layer > m_energyVsMuonLayer.size())
  {
    cout << "CaloCluster::getEnergyInMuonLayer : layer outside of range, returning 0.0" << endl;
    return 0;
  }
  else
  {
    return m_energyVsMuonLayer[layer];
  }
}

void CaloCluster::setEnergyVsECalBarrelLayers(const std::vector<float> &energyVsECalBarrelLayer)
{
  m_energyVsECalBarrelLayer = std::vector<float>(energyVsECalBarrelLayer); // use copy constructor
}

void CaloCluster::setEnergyVsECalEndCapLayers(const std::vector<float> &energyVsECalEndCapLayer)
{
  m_energyVsECalEndCapLayer = std::vector<float>(energyVsECalEndCapLayer); // use copy constructor
}

void CaloCluster::setEnergyVsHCalLayers(const std::vector<float> &energyVsHCalLayer)
{
  m_energyVsHCalLayer = std::vector<float>(energyVsHCalLayer); // use copy constructor
}

void CaloCluster::setEnergyVsMuonLayers(const std::vector<float> &energyVsMuonLayer)
{
  m_energyVsMuonLayer = std::vector<float>(energyVsMuonLayer); // use copy constructor
}

TVector3 CaloCluster::getBarycenterInECalBarrelLayer(unsigned int layer) const
{
  if (layer > m_barycenterVsECalBarrelLayer.size())
  {
    cout << "CaloCluster::getBarycenterInECalBarrelLayer : layer outside of range, returning origin" << endl;
    return TVector3(0.0, 0.0, 0.0);
  }
  else
  {
    return m_barycenterVsECalBarrelLayer[layer];
  }
}

TVector3 CaloCluster::getBarycenterInECalEndCapLayer(unsigned int layer) const
{
  if (layer > m_barycenterVsECalEndCapLayer.size())
  {
    cout << "CaloCluster::getBarycenterInECalEndCapLayer : layer outside of range, returning origin" << endl;
    return TVector3(0.0, 0.0, 0.0);
  }
  else
  {
    return m_barycenterVsECalEndCapLayer[layer];
  }
}

TVector3 CaloCluster::getBarycenterInHCalLayer(unsigned int layer) const
{
  if (layer > m_barycenterVsHCalLayer.size())
  {
    cout << "CaloCluster::getBarycenterInHCalLayer : layer outside of range, returning origin" << endl;
    return TVector3(0.0, 0.0, 0.0);
  }
  else
  {
    return m_barycenterVsHCalLayer[layer];
  }
}

TVector3 CaloCluster::getBarycenterInMuonLayer(unsigned int layer) const
{
  if (layer > m_barycenterVsMuonLayer.size())
  {
    cout << "CaloCluster::getEnergyInMuonLayer : layer outside of range, returning origin" << endl;
    return TVector3(0.0, 0.0, 0.0);
  }
  else
  {
    return m_barycenterVsMuonLayer[layer];
  }
}

void CaloCluster::setBarycenterVsECalBarrelLayers(const std::vector<TVector3> &barycenterVsECalBarrelLayer)
{
  m_barycenterVsECalBarrelLayer = std::vector<TVector3>(barycenterVsECalBarrelLayer); // use copy constructor
}

void CaloCluster::setBarycenterVsECalEndCapLayers(const std::vector<TVector3> &barycenterVsECalEndCapLayer)
{
  m_barycenterVsECalEndCapLayer = std::vector<TVector3>(barycenterVsECalEndCapLayer); // use copy constructor
}

void CaloCluster::setBarycenterVsHCalLayers(const std::vector<TVector3> &barycenterVsHCalLayer)
{
  m_barycenterVsHCalLayer = std::vector<TVector3>(barycenterVsHCalLayer); // use copy constructor
}

void CaloCluster::setBarycenterVsMuonLayers(const std::vector<TVector3> &barycenterVsMuonLayer)
{
  m_barycenterVsMuonLayer = std::vector<TVector3>(barycenterVsMuonLayer); // use copy constructor
}

void CaloCluster::print() const
{
  std::cout << "Cluster: " << m_index << std::endl;
  std::cout << "Energy: " << m_energy << " GeV" << endl;
  std::cout << "Barycenter (r, theta, phi): " << m_barycenter.Pt() << " " << m_barycenter.Theta() << " " << m_barycenter.Phi() << std::endl;
  std::cout << "Energy in ECAL: " << getEnergyECal() << " GeV" << endl;
  std::cout << "Energy in HCAL: " << getEnergyHCal() << " GeV" << endl;
  std::cout << "Energy in MUON: " << getEnergyMuon() << " GeV" << endl;
  std::cout << "Energy and barycenter in each layer of the ECAL barrel: " << std::endl;
  for (unsigned int i=0; i<m_energyVsECalBarrelLayer.size(); i++)
  {
    std::cout << i << "   " << m_energyVsECalBarrelLayer[i] << " GeV  "  <<   m_barycenterVsECalBarrelLayer[i].Pt() << " " << m_barycenterVsECalBarrelLayer[i].Theta() << " " << m_barycenterVsECalBarrelLayer[i].Phi() << std::endl;
  }
  std::cout << "Energy and barycenter in each layer of the ECAL endcap: " << std::endl;
  for (unsigned int i=0; i<m_energyVsECalEndCapLayer.size(); i++)
  {
    std::cout << i << "   " << m_energyVsECalEndCapLayer[i] << " GeV  "  <<   m_barycenterVsECalEndCapLayer[i].Pt() << " " << m_barycenterVsECalEndCapLayer[i].Theta() << " " << m_barycenterVsECalEndCapLayer[i].Phi() << std::endl;
  }
  std::cout << "Energy and barycenter in each layer of the HCAL: " << std::endl;
  for (unsigned int i=0; i<m_energyVsHCalLayer.size(); i++)
  {
    std::cout << i << "   " << m_energyVsHCalLayer[i] << " GeV  "  <<   m_barycenterVsHCalLayer[i].Pt() << " " << m_barycenterVsHCalLayer[i].Theta() << " " << m_barycenterVsHCalLayer[i].Phi() << std::endl;
  }
  std::cout << "Energy and barycenter in each layer of the MUON: " << std::endl;
  for (unsigned int i=0; i<m_energyVsMuonLayer.size(); i++)
  {
    std::cout << i << "   " << m_energyVsMuonLayer[i] << " GeV  "  <<   m_barycenterVsMuonLayer[i].Pt() << " " << m_barycenterVsMuonLayer[i].Theta() << " " << m_barycenterVsMuonLayer[i].Phi() << std::endl;
  }
  // ideally divide by cm (or mm)
  // and report output in cm
}
