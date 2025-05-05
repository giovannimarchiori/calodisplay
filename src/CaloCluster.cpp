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
  for (float Elayer : m_energyVsECalLayer)
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

float CaloCluster::getEnergyInECalLayer(unsigned int layer) const
{
  if (layer > m_energyVsECalLayer.size())
  {
    cout << "CaloCluster::getEnergyInECalLayer : layer outside of range, returning 0.0" << endl;
    return 0;
  }
  else
  {
    return m_energyVsECalLayer[layer];
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

void CaloCluster::setEnergyVsECalLayers(const std::vector<float> &energyVsECalLayer)
{
  m_energyVsECalLayer = std::vector<float>(energyVsECalLayer); // use copy constructor
}

void CaloCluster::setEnergyVsHCalLayers(const std::vector<float> &energyVsHCalLayer)
{
  m_energyVsHCalLayer = std::vector<float>(energyVsHCalLayer); // use copy constructor
}

void CaloCluster::setEnergyVsMuonLayers(const std::vector<float> &energyVsMuonLayer)
{
  m_energyVsMuonLayer = std::vector<float>(energyVsMuonLayer); // use copy constructor
}

TVector3 CaloCluster::getBarycenterInECalLayer(unsigned int layer) const
{
  if (layer > m_barycenterVsECalLayer.size())
  {
    cout << "CaloCluster::getEnergyInECalLayer : layer outside of range, returning origin" << endl;
    return TVector3(0.0, 0.0, 0.0);
  }
  else
  {
    return m_barycenterVsECalLayer[layer];
  }
}

TVector3 CaloCluster::getBarycenterInHCalLayer(unsigned int layer) const
{
  if (layer > m_barycenterVsHCalLayer.size())
  {
    cout << "CaloCluster::getEnergyInHCalLayer : layer outside of range, returning origin" << endl;
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

void CaloCluster::setBarycenterVsECalLayers(const std::vector<TVector3> &barycenterVsECalLayer)
{
  m_barycenterVsECalLayer = std::vector<TVector3>(barycenterVsECalLayer); // use copy constructor
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
  std::cout << "Energy and barycenter in each layer of the ECAL: " << std::endl;
  for (unsigned int i=0; i<m_energyVsECalLayer.size(); i++)
  {
    std::cout << i << "   " << m_energyVsECalLayer[i] << " GeV  "  <<   m_barycenterVsECalLayer[i].Pt() << " " << m_barycenterVsECalLayer[i].Theta() << " " << m_barycenterVsECalLayer[i].Phi() << std::endl;
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
