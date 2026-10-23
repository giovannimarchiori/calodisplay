/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
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

namespace
{
// names of the parts, for the printouts
const char *partName[CaloCluster::kNParts] = {"ECAL barrel", "ECAL endcap", "HCAL barrel", "HCAL endcap", "MUON"};
}

/******************************************************************************/
// class definition
/******************************************************************************/

float CaloCluster::getEnergy(Part part) const
{
  float E(0.);
  for (float Elayer : m_layers[part].energy)
  {
    E += Elayer;
  }
  return E;
}

float CaloCluster::getEnergyInLayer(Part part, unsigned int layer) const
{
  if (layer >= m_layers[part].energy.size())
  {
    cout << "CaloCluster::getEnergyInLayer : layer " << layer << " of the " << partName[part]
         << " outside of range, returning 0.0" << endl;
    return 0;
  }
  return m_layers[part].energy[layer];
}

TVector3 CaloCluster::getBarycenterInLayer(Part part, unsigned int layer) const
{
  if (layer >= m_layers[part].barycenter.size())
  {
    cout << "CaloCluster::getBarycenterInLayer : layer " << layer << " of the " << partName[part]
         << " outside of range, returning origin" << endl;
    return TVector3(0.0, 0.0, 0.0);
  }
  return m_layers[part].barycenter[layer];
}

void CaloCluster::print() const
{
  std::cout << "Cluster: " << m_index << std::endl;
  std::cout << "Energy: " << m_energy << " GeV" << endl;
  std::cout << "Barycenter (r, theta, phi): " << m_barycenter.Pt() << " " << m_barycenter.Theta() << " " << m_barycenter.Phi() << std::endl;
  std::cout << "Energy in ECAL: " << getEnergyECal() << " GeV" << endl;
  std::cout << "Energy in HCAL: " << getEnergyHCal() << " GeV" << endl;
  std::cout << "Energy in MUON: " << getEnergyMuon() << " GeV" << endl;
  for (int part = 0; part < kNParts; part++)
  {
    const Layers &layers = m_layers[part];
    std::cout << "Energy and barycenter in each layer of the " << partName[part] << ": " << std::endl;
    for (unsigned int i = 0; i < layers.energy.size(); i++)
    {
      std::cout << i << "   " << layers.energy[i] << " GeV  " << layers.barycenter[i].Pt() << " " << layers.barycenter[i].Theta() << " " << layers.barycenter[i].Phi() << std::endl;
    }
  }
  // ideally divide by cm (or mm)
  // and report output in cm
}
