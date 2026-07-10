// ElectronExtractor.cpp

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>
#include <limits>

constexpr unsigned int INVALID_INDEX =
    std::numeric_limits<unsigned int>::max();

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include <TTreeReader.h>
#include <TTreeReaderArray.h>


//============================================================
// Lightweight MC particle representation
//============================================================

struct Particle {

    int pdg = 0;
    double mass = 0.0;
    double time = 0.0;

    TVector3 vertex;
    TVector3 endpoint;
    TVector3 momentum;

    unsigned int daughters_begin = 0;
    unsigned int daughters_end   = 0;

    unsigned int parents_begin = 0;
    unsigned int parents_end   = 0;
};


struct ElectronSegment {

    TVector3 vertex;
    TVector3 endpoint;

    TLorentzVector p4;

    unsigned int sourceIndex = INVALID_INDEX;
};


//============================================================
// Helpers
//============================================================

void printParticle(
		   unsigned int idx,
		   const std::vector<Particle>& particles,
		   const TTreeReaderArray<Int_t>& daughters,
		   int depth=0)
{
  if(idx>=particles.size())
    return;
  if(depth>20)
    {
      std::cout<<"max depth reached\n";
      return;
    }
  const auto& p=particles[idx];
  std::cout
    << std::string(depth*2,' ')
    << idx
    << " PDG="
    << p.pdg
    << " daughters=";
  for(unsigned int i=p.daughters_begin;
      i<p.daughters_end;
      ++i)
    {
      if(i>=daughters.GetSize())
	continue;
      std::cout << daughters[i] << " ";
    }
  std::cout<<"\n";
  for(unsigned int i=p.daughters_begin;
      i<p.daughters_end;
      ++i)
    {
      if(i>=daughters.GetSize())
	continue;
      int child=daughters[i];
      if(child<0)
	continue;
      unsigned int uchild = static_cast<unsigned int>(child);
      if(uchild>=particles.size())
	continue;
      printParticle(
		    uchild,
		    particles,
		    daughters,
		    depth+1);
    }
}

static TLorentzVector makeP4(const Particle& p)
{
    double E =
        std::sqrt(
            p.momentum.Mag2()
            +
            p.mass*p.mass
        );

    return TLorentzVector(
        p.momentum.X(),
        p.momentum.Y(),
        p.momentum.Z(),
        E
    );
}


static TLorentzVector makePhotonP4(const Particle& p)
{
    double E = p.momentum.Mag();

    return TLorentzVector(
        p.momentum.X(),
        p.momentum.Y(),
        p.momentum.Z(),
        E
    );
}


//------------------------------------------------------------
// Return photons ordered by distance along electron direction
//------------------------------------------------------------

static std::vector<unsigned int>
orderedPhotonDaughters(
		       unsigned int electronIndex,
		       const std::vector<Particle>& particles,
		       const std::vector<int>& daughters)
{
  const auto& e = particles[electronIndex];

  std::vector<unsigned int> photons;
  for(unsigned int i=e.daughters_begin; i<e.daughters_end; ++i)
    {
      if(i>=daughters.size())
	continue;

      int idx = daughters[i];
      if(idx<0)
	continue;
      unsigned int uidx = static_cast<unsigned int>(idx);
      if(uidx>=particles.size())
	continue;

      if(particles[uidx].pdg==22)
	photons.push_back(uidx);
    }

  TVector3 direction = e.momentum.Unit();
  TVector3 origin = e.vertex;
  
  std::sort(photons.begin(), photons.end(),
	    [&](unsigned int a, unsigned int b)
	    {
	      double da =
                (particles[a].vertex-origin).Dot(direction);
	      
	      double db =
                (particles[b].vertex-origin).Dot(direction);
	      
	      return da < db;
	    });
  
  return photons;
}



//------------------------------------------------------------
// Find daughter electron
//------------------------------------------------------------

static unsigned int electronDaughter(
				     unsigned int electronIndex,
				     const std::vector<Particle>& particles,
				     const std::vector<int>& daughters)
{
  const auto& e = particles[electronIndex];

  for(unsigned int i=e.daughters_begin;
      i<e.daughters_end;
      ++i)
    {
      if(i>=daughters.size())
	continue;

        int idx=daughters[i];
	if(idx<0)
	  continue;
	unsigned int uidx =
	  static_cast<unsigned int>(idx);
        if(uidx>=particles.size())
	  continue;

        if(particles[uidx].pdg==11)
	  return uidx;
    }

  return INVALID_INDEX;
}



//============================================================
// Recursive electron traversal
//============================================================

static void followElectron(
			   unsigned int electronIndex,
			   const std::vector<Particle>& particles,
			   const std::vector<int>& daughters,
			   std::vector<ElectronSegment>& output,
			   std::set<unsigned int>& visited)
{
    if(electronIndex>=particles.size())
        return;

    // protect against malformed loops
    if(visited.count(electronIndex))
    {
        std::cerr
            << "WARNING: electron loop detected at "
            << electronIndex
            << std::endl;
        return;
    }

    visited.insert(electronIndex);

    const Particle& e =
        particles[electronIndex];

    TLorentzVector currentP4 =
        makeP4(e);


    //--------------------------------------------------------
    // Look for daughters
    //--------------------------------------------------------

    unsigned int nextElectron =
      electronDaughter(
		       electronIndex,
		       particles,
		       daughters);


    std::vector<unsigned int> photons =
      orderedPhotonDaughters(
			     electronIndex,
			     particles,
			     daughters);
    


    //--------------------------------------------------------
    // Case 1:
    // no photon, no continuation
    //--------------------------------------------------------

    if(photons.empty() &&
       nextElectron==INVALID_INDEX)
    {
        output.push_back(
            {
                e.vertex,
                e.endpoint,
                currentP4,
                electronIndex
            });

        return;
    }



    //--------------------------------------------------------
    // Case 2:
    // photon(s) present
    //--------------------------------------------------------

    TVector3 currentVertex =
        e.vertex;


    for(unsigned int phIndex : photons)
    {
        const Particle& gamma =
            particles[phIndex];

        output.push_back(
            {
                currentVertex,
                gamma.vertex,
                currentP4,
                electronIndex
            });


        //
        // Only subtract photons if there is no
        // daughter electron continuation
        //

        if(nextElectron==INVALID_INDEX)
        {
            currentP4 -=
                makePhotonP4(gamma);
        }


        currentVertex =
            gamma.vertex;
    }


    //--------------------------------------------------------
    // Continue electron chain
    //--------------------------------------------------------

    if(nextElectron != INVALID_INDEX)
    {

        followElectron(
            nextElectron,
            particles,
            daughters,
            output,
            visited);

    }
    else
    {

        //
        // no daughter electron:
        // create final segment
        //

        output.push_back(
            {
                currentVertex,
                e.endpoint,
                currentP4,
                electronIndex
            });
    }
}



//============================================================
// Identify primary electron chains
//============================================================

static bool hasElectronParent(
			      unsigned int index,
			      const std::vector<Particle>& particles,
			      const std::vector<int>& parents)
{
  const auto& p = particles[index];

  for(unsigned int i=p.parents_begin; i<p.parents_end; ++i)
    {
      if(i>=parents.size())
	continue;
      
      int parent = parents[i];
      if(parent < 0)
	continue;
      unsigned int uparent =
	static_cast<unsigned int>(parent);
      if(uparent<particles.size() &&
	 particles[uparent].pdg==11)
	return true;
    }

  return false;
}



//============================================================
// Main extraction function
//============================================================

std::vector<ElectronSegment> extractElectronSegments(
        const std::vector<Particle>& particles,
        const std::vector<int>& parents,
        const std::vector<int>& daughters)
{

    std::vector<ElectronSegment> result;
    result.reserve(10);
    
    std::set<unsigned int> visited;
    for(size_t i=0;i<particles.size();++i)
    {

        if(particles[i].pdg!=11)
            continue;


        //
        // start only from electrons which are not
        // descendants of another electron
        //

        if(hasElectronParent(
                i,
                particles,
                parents))
            continue;



        followElectron(
            i,
            particles,
            daughters,
            result,
            visited);
    }


    return result;
}


int test() {
  TFile* f = TFile::Open("data/allegro_o1_v03_evts_10_pdg_11_MomentumMinMax_10.0_10.0_GeV_ThetaMinMax_90.0_90.0_PhiMinMax_0_360_HCal_ON_digi_reco.root");
  if(!f || f->IsZombie())
    {
      std::cerr<<"Cannot open file\n";
      return 1;
    }

  TTreeReader fReader("events", f);
  
  TTreeReaderArray<Int_t> MCParticles_PDG = {fReader, "MCParticles.PDG"};
  TTreeReaderArray<Int_t> MCParticles_generatorStatus = {fReader, "MCParticles.generatorStatus"};
  TTreeReaderArray<Int_t> MCParticles_simulatorStatus = {fReader, "MCParticles.simulatorStatus"};
  TTreeReaderArray<Float_t> MCParticles_charge = {fReader, "MCParticles.charge"};
  TTreeReaderArray<Float_t> MCParticles_time = {fReader, "MCParticles.time"};
  TTreeReaderArray<Double_t> MCParticles_mass = {fReader, "MCParticles.mass"};
  TTreeReaderArray<Double_t> MCParticles_vertex_x = {fReader, "MCParticles.vertex.x"};
  TTreeReaderArray<Double_t> MCParticles_vertex_y = {fReader, "MCParticles.vertex.y"};
  TTreeReaderArray<Double_t> MCParticles_vertex_z = {fReader, "MCParticles.vertex.z"};
  TTreeReaderArray<Double_t> MCParticles_endpoint_x = {fReader, "MCParticles.endpoint.x"};
  TTreeReaderArray<Double_t> MCParticles_endpoint_y = {fReader, "MCParticles.endpoint.y"};
  TTreeReaderArray<Double_t> MCParticles_endpoint_z = {fReader, "MCParticles.endpoint.z"};
  TTreeReaderArray<Double_t> MCParticles_momentum_x = {fReader, "MCParticles.momentum.x"};
  TTreeReaderArray<Double_t> MCParticles_momentum_y = {fReader, "MCParticles.momentum.y"};
  TTreeReaderArray<Double_t> MCParticles_momentum_z = {fReader, "MCParticles.momentum.z"};
  TTreeReaderArray<UInt_t> MCParticles_parents_begin = {fReader, "MCParticles.parents_begin"};
  TTreeReaderArray<UInt_t> MCParticles_parents_end = {fReader, "MCParticles.parents_end"};
  TTreeReaderArray<Int_t> _MCParticles_parents_index = {fReader, "_MCParticles_parents.index"};
  TTreeReaderArray<UInt_t> MCParticles_daughters_begin = {fReader, "MCParticles.daughters_begin"};
  TTreeReaderArray<UInt_t> MCParticles_daughters_end = {fReader, "MCParticles.daughters_end"};
  TTreeReaderArray<Int_t> _MCParticles_daughters_index = {fReader, "_MCParticles_daughters.index"};
  
  Long64_t nEntries = fReader.GetEntries();
  for(Long64_t event=0; event<nEntries; ++event)
    {
      std::cout << std::endl;
      std::cout << "######################" << std::endl;
      std::cout << "EVENT " << event << std::endl;
      std::cout << "######################" << std::endl;
      fReader.SetEntry(event);
      
      std::vector<Particle> particles;
      particles.reserve(MCParticles_PDG.GetSize());
      for(size_t i=0;i<MCParticles_PDG.GetSize();i++)
	{
	  Particle p;
	  p.pdg = MCParticles_PDG[i];
	  p.mass = MCParticles_mass[i];
	  p.time = MCParticles_time[i];
	  p.vertex.SetXYZ(
			  MCParticles_vertex_x[i],
			  MCParticles_vertex_y[i],
			  MCParticles_vertex_z[i]);
	  p.endpoint.SetXYZ(
			    MCParticles_endpoint_x[i],
			    MCParticles_endpoint_y[i],
			   MCParticles_endpoint_z[i]);
	  p.momentum.SetXYZ(
			    MCParticles_momentum_x[i],
			    MCParticles_momentum_y[i],
			    MCParticles_momentum_z[i]);
	  p.daughters_begin =
	    MCParticles_daughters_begin[i];
	  p.daughters_end =
	    MCParticles_daughters_end[i];
	  p.parents_begin =
	    MCParticles_parents_begin[i];
	  p.parents_end =
	    MCParticles_parents_end[i];
	  particles.push_back(p);
	}

      // for debug, print particles
      printParticle(
		    0,
		    particles,
		    _MCParticles_daughters_index);

      std::vector<int> parents;
      std::vector<int> daughters;
      parents.reserve(_MCParticles_parents_index.GetSize());
      daughters.reserve(_MCParticles_daughters_index.GetSize());
      for (unsigned int i=0; i<_MCParticles_parents_index.GetSize(); ++i)
	parents.push_back(_MCParticles_parents_index[i]);
      for (unsigned int i=0; i<_MCParticles_daughters_index.GetSize(); ++i)
	daughters.push_back(_MCParticles_daughters_index[i]);

      auto electrons =
	extractElectronSegments(
				particles,
				parents,
				daughters);
      
      for(const auto& e : electrons)
	{
	  std::cout
	    << "Electron "
	    << e.sourceIndex << std::endl
	    << "from "
	    << e.vertex.X() << " "
	    << e.vertex.Y() << " "
	    << e.vertex.Z() << " (R="
	    << e.vertex.Mag() 
	    << ") to "
	    << e.endpoint.X() << " "
	    << e.endpoint.Y() << " "
	    << e.endpoint.Z() << " (R="
	    << e.endpoint.Mag()
	    << ") p = "
	    << e.p4.X() << " "
	    << e.p4.Y() << " "
	    << e.p4.Z() << " (|p|="
	    << e.p4.P() << ")"
	    << std::endl;
	}
    }
    return 0;
}
