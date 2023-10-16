/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventDisplay: creates and populates the event display
//
/******************************************************************************/

#include "EventDisplay.h"
#include "DetectorGeometry.h"

#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPRegexp.h>

#include <TGButton.h>
#include <TGTab.h>

#include <TGeoTube.h>
#include <TGeoMatrix.h>

#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TEveStraightLineSet.h>
#include <TEveTrackPropagator.h>
#include <TEveProjectionAxes.h>

#include <filesystem>

// return the sign of a float
int sgn(float val) {
  return (val > 0.) - (val < 0.);
}


/******************************************************************************/
// HELPER FUNCTIONS related to the graphic system
/******************************************************************************/

// derived TGLAnnotation class that overrides MouseEnter method
// to avoid editing of annotation
TGLConstAnnotation::TGLConstAnnotation(TGLViewerBase *parent, const char *text, Float_t posx, Float_t posy) :
  TGLAnnotation(parent, text, posx, posy)
{
  ;
}

Bool_t TGLConstAnnotation::MouseEnter(TGLOvlSelectRecord& /*rec*/)
{
  fActive = kFALSE;
  return kTRUE;
}


void EventDisplay::loadEvent(int event) {
  
  if (event != -1) eventId = event;
  
  printf("Loading event %d ...\n", eventId);
  textEntry->SetTextColor(0xff0000);
  textEntry->SetText(Form("Loading event %d ...",eventId));

  eventReader->loadEvent(eventId);
  
  TString partType;
  double pmax=0.0;
  int ipmax=-1;
  for (unsigned int i = 0; i < eventReader->genParticles_generatorStatus->GetSize(); i ++) {
    if ( (* eventReader->genParticles_generatorStatus)[i] != 1 ) continue;
    float _px = (* eventReader->genParticles_momentum_x)[i];
    float _py = (* eventReader->genParticles_momentum_y)[i];
    float _pz = (* eventReader->genParticles_momentum_z)[i];
    float _p = sqrt( _px*_px + _py*_py + _pz*_pz );
    //cout << p << endl;
    if (_p==0.) continue;
    if (_p>pmax) {
      pmax = _p;
      ipmax = i;
    }
  }
  int pdgID = (* eventReader->genParticles_PDG)[ipmax];
  if (pdgID == 111) partType = "pi0";
  else if (pdgID == 22) partType = "y";
  else if (pdgID == 11) partType = "e-";
  else if (pdgID == -11) partType = "e+";
  else partType = "unknown";

  const double cm = geomReader->cm;
  const double mm = geomReader->mm;
  const double rMin = geomReader->rMin;
  const double rMax = geomReader->rMax;
  const double alpha = geomReader->alpha;
  const double thetaGrid = geomReader->thetaGrid;
  const double gridPhi = geomReader->gridPhi;
  const double phiMin = geomReader->phiMin;
  
  //
  // particles
  //
  if (drawParticles) {
    cout << "Creating particles" << endl;
    if (particles == nullptr) {
      particles = new TEveTrackList("particles");
      TEveTrackPropagator* trkProp = particles->GetPropagator();
      trkProp->SetMagField( 0.01 );
      particles->SetMainColor(kYellow);
      particles->SetLineWidth(2);
      gEve->AddElement(particles);
    }
    else
      particles->DestroyElements();
    
    // handle differently the e/gamma vs pi0 particle guns
    // (for pi0, need to look for the two photons in secondary particle list
    // unfortunately cannot just use the latter to show particles since
    // info about endpoint is broken (identical to vertex)
    float m = (* eventReader->genParticles_mass)[ipmax];
    float px = (* eventReader->genParticles_momentum_x)[ipmax];
    float py = (* eventReader->genParticles_momentum_y)[ipmax];
    float pz = (* eventReader->genParticles_momentum_z)[ipmax];
    float p = sqrt( px*px + py*py + pz*pz );
    
    double t = (* eventReader->genParticles_time)[ipmax];
    double x1 = (* eventReader->genParticles_vertex_x)[ipmax] * mm;
    double y1 = (* eventReader->genParticles_vertex_y)[ipmax] * mm;
    double z1 = (* eventReader->genParticles_vertex_z)[ipmax] * mm;
    double r1 = sqrt(x1*x1+y1*y1);
    double sintheta = sqrt(px*px + py*py)/p;
    double x2 = x1 + px/p * rMax / sintheta;
    double y2 = y1 + py/p * rMax / sintheta;
    double z2 = z1 + pz/p * rMax / sintheta;
    double r2 = sqrt(x2*x2+y2*y2);
    
    TEveMCTrack mct;
    mct.SetPdgCode( pdgID );
    mct.SetMomentum( px, py, pz, sqrt(p*p + m*m) );
    mct.SetProductionVertex( x1, y1, z1, t );
    TEveTrack* track = new TEveTrack(&mct, particles->GetPropagator());
    track->SetAttLineAttMarker(particles);
    track->SetElementTitle(Form("p = %.3f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
				p, acos(pz/p), atan2(py,px),
				x1 / cm, y1 / cm, z1 / cm));
    particles->AddElement(track);
    
    // if the particle is a pi0, also draw the two photons, and set the endpoint
    // of the pi0 track
    if (pdgID == 111) {
      bool decayVtxSet = false;
      for (unsigned int i = 0; i < eventReader->SimParticleSecondaries_PDG->GetSize(); i ++) {
	pdgID = (* eventReader->SimParticleSecondaries_PDG)[i];
	// keep only photons
	if (pdgID!=22) continue;
	//if ( (* eventReader->SimParticleSecondaries_generatorStatus)[i] != 1 ) continue;
	px = (* eventReader->SimParticleSecondaries_momentum_x)[i];
	py = (* eventReader->SimParticleSecondaries_momentum_y)[i];
	pz = (* eventReader->SimParticleSecondaries_momentum_z)[i];
	m = (* eventReader->SimParticleSecondaries_mass)[i];
	p = sqrt( px*px + py*py + pz*pz );
	float e = sqrt(p*p + m*m);
	// cout << "p = "<< p << endl;
	if (p<ParticleEnergyThreshold) continue;
	
	// cout << "PDG = "<< pdgID << endl;
	t = (* eventReader->SimParticleSecondaries_time)[i];
	x1 = (* eventReader->SimParticleSecondaries_vertex_x)[i] * mm;
	y1 = (* eventReader->SimParticleSecondaries_vertex_y)[i] * mm;
	z1 = (* eventReader->SimParticleSecondaries_vertex_z)[i] * mm;
	r1 = sqrt(x1*x1+y1*y1);
	// the two photons from a pi0 in the origin must come from small R
	if (r1>1.) continue;
	sintheta = sqrt(px*px + py*py)/p;
	x2 = x1 + px/p * rMax / sintheta;
	y2 = y1 + py/p * rMax / sintheta;
	z2 = z1 + pz/p * rMax / sintheta;
	//double x2 = (* eventReader->SimParticleSecondaries_endpoint_x)[i] * mm;
	//double y2 = (* eventReader->SimParticleSecondaries_endpoint_y)[i] * mm;
	//double z2 = (* eventReader->SimParticleSecondaries_endpoint_z)[i] * mm;
	r2 = sqrt(x2*x2+y2*y2);
	// cout << "x1 y1 z1 x2 y2 z2 = "
	//      << x1 << " " << y1 << " " << z1 << " "
	//      << x2 << " " << y2 << " " << z2 << endl;
	// set pi0 decay point
	if (!decayVtxSet) {
	  TEveVectorF v; v[0] = x1; v[1]=y1; v[2]=z1;
	  TEvePathMark mark(TEvePathMark::kDecay, v);
	  track->AddPathMark(mark);
	  decayVtxSet = true;
	}
	TEveMCTrack mct;
	mct.SetPdgCode( pdgID );
	mct.SetMomentum( px, py, pz, sqrt(p*p + m*m) );
	mct.SetProductionVertex( x1, y1, z1, t );
	TEveTrack* track = new TEveTrack(&mct, particles->GetPropagator());
	track->SetAttLineAttMarker(particles);
	track->SetElementTitle(Form("p = .3%f GeV\ntheta = %f\nphi = %f\nx = %f cm\ny = %f cm\nz= %f cm",
				    p, acos(pz/p), atan2(py,px),
				    x1 / cm, y1 / cm, z1 / cm));
	
	particles->AddElement(track);
      }
    }
    particles->MakeTracks();
  }

  //
  // hits
  //
  if (drawHits) {
    cout << "Creating hits" << endl;
    if (hits == nullptr) {
      hits = new TEvePointSet();
      hits->SetName(Form("hits (E>%.1f GeV)",CellEnergyThreshold));
      hits->SetMarkerStyle(4);
      hits->SetMarkerSize(1);
      hits->SetMarkerColor(kRed);
      gEve->AddElement(hits);
    }
    else
      hits->Reset();
    for (unsigned int i = 0; i < eventReader->ECalBarrelPositionedHits_position_x->GetSize(); i ++) {
      float E = (*eventReader->ECalBarrelPositionedHits_energy)[i];
      if (E<HitEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->ECalBarrelPositionedHits_cellID)[i];
      // ULong_t layer = DetectorGeometry:::Layer(cellID);
      hits->SetNextPoint( (*eventReader->ECalBarrelPositionedHits_position_x)[i] * mm ,
			  (*eventReader->ECalBarrelPositionedHits_position_y)[i] * mm ,
			  (*eventReader->ECalBarrelPositionedHits_position_z)[i] * mm );
    }
  }
  
  //
  // cells
  //
  if (drawCells) {
    cout << "Creating cells" << endl; 
    if (cells == nullptr) {
      cells = new TEvePointSet();
      cells->SetName(Form("cells (E>%.1f GeV)",CellEnergyThreshold));
      cells->SetMarkerStyle(4);
      cells->SetMarkerSize(3);
      cells->SetMarkerColor(kBlue);
      gEve->AddElement(cells);
    }
    else
      cells->Reset();
    for (unsigned int i = 0; i < eventReader->ECalBarrelPositionedCells_position_x->GetSize(); i ++) {
      float E = (*eventReader->ECalBarrelPositionedCells_energy)[i];
      if (E<CellEnergyThreshold) continue;
      // ULong_t cellID = (*eventReader->ECalBarrelPositionedCells_cellID)[i];
      // ULong_t layer = DetectorGeometry::Layer(cellID);
      cells->SetNextPoint( (*eventReader->ECalBarrelPositionedCells_position_x)[i] * mm ,
			   (*eventReader->ECalBarrelPositionedCells_position_y)[i] * mm ,
			     (*eventReader->ECalBarrelPositionedCells_position_z)[i] * mm );
    }
  }
  
  //
  // cells merged 
  //
  if (drawMergedCells) {
    cout << "Creating merged cells" << endl;
    if (cells_merged == nullptr) {
      cells_merged = new TEvePointSet();
      cells_merged->SetName("cells_merged");
      cells_merged->SetMarkerStyle(4);
      cells_merged->SetMarkerSize(3);
      cells_merged->SetMarkerColor(kBlue);
      gEve->AddElement(cells_merged);
    }
    else
      cells_merged->Reset();
    for (unsigned int i = 0; i < eventReader->ECalBarrelPositionedCells2_position_x->GetSize(); i ++) {
      float E = (*eventReader->ECalBarrelPositionedCells2_energy)[i];
      // if (E<minCellE) continue;
      // ULong_t cellID = (*eventReader->ECalBarrelPositionedCells_cellID)[i];
      // ULong_t layer = DetectorGeometry::Layer(cellID);
      cells_merged->SetNextPoint( (*eventReader->ECalBarrelPositionedCells2_position_x)[i] * mm ,
				  (*eventReader->ECalBarrelPositionedCells2_position_y)[i] * mm ,
				  (*eventReader->ECalBarrelPositionedCells2_position_z)[i] * mm );
    }
  }


  //
  // clusters
  //
  if (drawTopoClusters) {
    cout << "Creating clusters" << endl;

    // centers of the clusters
    if (clusters == nullptr) {
      clusters = new TEvePointSet();
      clusters->SetName(Form("clusters (E>%.1f GeV)",TopoClusterEnergyThreshold));
      clusters->SetMarkerStyle(4);
      clusters->SetMarkerSize(6);
      clusters->SetMarkerColor(kGreen);
      gEve->AddElement(clusters);
    }
    else
      clusters->Reset();
      
    for (unsigned int i = 0; i < eventReader->CorrectedCaloTopoClusters_position_x->GetSize(); i ++) {
      float E = (*eventReader->CorrectedCaloTopoClusters_energy)[i];
      if (E < TopoClusterEnergyThreshold) continue;
      // cluster positions are in cm and hits/cells in mm ...
      clusters->SetNextPoint( (*eventReader->CorrectedCaloTopoClusters_position_x)[i] * cm ,
			      (*eventReader->CorrectedCaloTopoClusters_position_y)[i] * cm ,
			      (*eventReader->CorrectedCaloTopoClusters_position_z)[i] * cm ); 
    }

      
    // the DestroyElements method does not work for TExeQuadSet as for the TEvePointSet..
    // so I have to destroy and recreate the quad and box sets
    for (auto qs : qs_rhoz) {
      if (qs)
	qs->Destroy();
    }
    qs_rhoz.clear();
      
    for (auto qs : qs_rhophi) {
      if (qs)
	qs->Destroy();
    }
    qs_rhophi.clear();

    for (auto _bs : bs) {
      if (_bs)
	_bs->Destroy();
    }
    bs.clear();

    TEveRGBAPalette *pal = new TEveRGBAPalette(0, 1000);
      
    // clusters in 3D
    if (topoclusters_3D==nullptr) {
      topoclusters_3D = new TEveElementList("Clusters in 3D (no E cut)");
      gEve->AddElement(topoclusters_3D);
    }
    else
      topoclusters_3D->DestroyElements();
      
    // clusters in 2D
    if (topoclusters_rhoz==nullptr) {
      topoclusters_rhoz = new TEveElementList(Form("Clusters in rho-z (E>%.1f GeV)",
						   TopoClusterEnergyThreshold));
      // add to scene or manager?
      // rhoZProjManager->AddElement(topoclusters_rhoz);
      rhoZEventSceneManual->AddElement(topoclusters_rhoz);
      gEve->AddToListTree(topoclusters_rhoz,false);
    }
    else
      topoclusters_rhoz->DestroyElements();

    if (topoclusters_rhophi==nullptr) {
      topoclusters_rhophi = new TEveElementList(Form("Clusters in rho-phi (E>%.1f GeV)",
						     TopoClusterEnergyThreshold));
      // add to scene or manager? -- scene that is not auto-projected!
      // rhoPhiProjManager->AddElement(topoclusters_rhophi);
      rhoPhiEventSceneManual->AddElement(topoclusters_rhophi);
      gEve->AddToListTree(topoclusters_rhophi,false);
    }
    else
      topoclusters_rhophi->DestroyElements();

    for (unsigned int i = 0; i < eventReader->CorrectedCaloTopoClusters_energy->GetSize(); i ++) {
      float energy = (*eventReader->CorrectedCaloTopoClusters_energy)[i];
      float xcl = (*eventReader->CorrectedCaloTopoClusters_position_x)[i];
      float ycl = (*eventReader->CorrectedCaloTopoClusters_position_y)[i];
      float zcl = (*eventReader->CorrectedCaloTopoClusters_position_z)[i];
      float rcl = sqrt(xcl*xcl + ycl*ycl);
      float phicl = atan2(ycl, xcl);
      float thetacl = atan2(rcl, zcl);

      if (energy < TopoClusterEnergyThreshold) {
	qs_rhoz.push_back(nullptr);
	qs_rhophi.push_back(nullptr);
      }
      else {
	TEveQuadSet* aqs = new TEveQuadSet(TEveQuadSet::kQT_FreeQuad, false, 32,
					   Form("cluster %d", (int) i));
	aqs->SetMainTransparency(80);
	// by calling SetOwnIds(kTRUE) the digit-set becomes
	// the owner of the assigned objects and deletes
	// them on destruction.
	aqs->SetOwnIds(kTRUE);
	aqs->SetPalette(pal);
	aqs->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",
			   energy,
			   rcl,
			   thetacl,
			   phicl));
	qs_rhoz.push_back(aqs);
	topoclusters_rhoz->AddElement(aqs);
	  
	TEveQuadSet* aqs2 = new TEveQuadSet(TEveQuadSet::kQT_FreeQuad, false, 32,
					    Form("cluster %d", (int) i));
	aqs2->SetMainTransparency(80);
	aqs2->SetOwnIds(kTRUE);
	aqs2->SetPalette(pal);
	aqs2->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",			    
			    energy,
			    rcl,
			    thetacl,
			    phicl));
	qs_rhophi.push_back(aqs2);
	topoclusters_rhophi->AddElement(aqs2);
      }
	
      TEveBoxSet* _bs = new TEveBoxSet(Form("cluster %d", (int) i),"");
      _bs->Reset(TEveBoxSet::kBT_FreeBox, false, 32);
      _bs->SetMainTransparency(80);
      _bs->SetOwnIds(kTRUE);
      _bs->SetPalette(pal);
      _bs->SetTitle(Form("E = %f GeV\nR = %f cm\ntheta = %f\nphi = %f",
			 energy,
			 rcl,
			 thetacl,
			 phicl));
      bs.push_back(_bs);
      topoclusters_3D->AddElement(_bs);
    }
      
    // loop over cells and attach them to clusters
    for (unsigned int i = 0; i < eventReader->PositionedCaloTopoClusterCells_energy->GetSize(); i ++) {
      int icl = -1;
      for (unsigned int j = 0; j < eventReader->CorrectedCaloTopoClusters_energy->GetSize(); j ++) {
	if (i >= (*eventReader->CorrectedCaloTopoClusters_hits_begin)[j] &&
	    i < (*eventReader->CorrectedCaloTopoClusters_hits_end)[j]) {
	  icl = j;
	  break;
	}
      }
      if (icl==-1) continue; // should never happen..
      ULong_t cellID = (*eventReader->PositionedCaloTopoClusterCells_cellID)[i];
      int layer = (int) DetectorGeometry::Layer(cellID);
      float x_center = (*eventReader->PositionedCaloTopoClusterCells_position_x)[i] * mm;
      float y_center = (*eventReader->PositionedCaloTopoClusterCells_position_y)[i] * mm;
      float z_center = (*eventReader->PositionedCaloTopoClusterCells_position_z)[i] * mm;
      float r_center = sqrt(x_center*x_center + y_center*y_center);
      float r_in = geomReader->r[layer];
      float r_out = geomReader->r[layer+1];
      float theta_center = atan2(r_center, z_center);
      //float phi_center = atan2(y_center, x_center);
      float energy = (*eventReader->PositionedCaloTopoClusterCells_energy)[i];
	
      // cluster cells in rho-z projection
      float verts[12];
      verts[0] = r_in / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
      verts[1] = r_in * sgn(y_center);
      verts[3] = r_in / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
      verts[4] = r_in * sgn(y_center);
      verts[6] = r_out / tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
      verts[7] = r_out * sgn(y_center);
      verts[9] = r_out / tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
      verts[10] = r_out * sgn(y_center);
      verts[2] = verts[5] = verts[8] = verts[11] = 0.;
      if (qs_rhoz[icl]!=nullptr) {
	qs_rhoz[icl]->AddQuad(verts);
	qs_rhoz[icl]->QuadValue( (int) (1000 * energy) );
	qs_rhoz[icl]->QuadId( new TNamed(Form("Cell %lu", cellID), "Dong!") );
      }
      // cluster cells in rho-phi projection
      int module = (int) DetectorGeometry::Module(cellID);
      double Lin = geomReader->getL(alpha, rMin, geomReader->r[layer]);
      double Lout = geomReader->getL(alpha, rMin, geomReader->r[layer+1]);
      double deltaL = rMin*sin(gridPhi/2.0)*sin(alpha);
      for (int j = 0; j < geomReader->mergedModules[layer]; j++) {
	int iModule = module + j;
	double phi0 = phiMin + iModule*gridPhi;
	verts[0] = rMin*cos(phi0) + (Lin+deltaL)*cos(phi0+alpha);
	verts[1] = rMin*sin(phi0) + (Lin+deltaL)*sin(phi0+alpha);
	verts[3] = rMin*cos(phi0) + (Lout+deltaL)*cos(phi0+alpha);
	verts[4] = rMin*sin(phi0) + (Lout+deltaL)*sin(phi0+alpha);
	verts[6] = rMin*cos(phi0+gridPhi) + (Lout-deltaL)*cos(phi0+gridPhi+alpha);
	verts[7] = rMin*sin(phi0+gridPhi) + (Lout-deltaL)*sin(phi0+gridPhi+alpha);
	verts[9] = rMin*cos(phi0+gridPhi) + (Lin-deltaL)*cos(phi0+gridPhi+alpha);
	verts[10] = rMin*sin(phi0+gridPhi) + (Lin-deltaL)*sin(phi0+gridPhi+alpha);
	verts[2] = verts[5] = verts[8] = verts[11] = 0.;
	if (qs_rhophi[icl]!=nullptr) {
	  qs_rhophi[icl]->AddQuad(verts);
	  qs_rhophi[icl]->QuadValue( (int) (1000 * energy) );
	  qs_rhophi[icl]->QuadId(new TNamed(Form("Cell %lu", cellID), "Dong!"));
	}
      }

      // cluster cells in 3D
      float verts3D[24];
      for (int j=0; j < geomReader->mergedModules[layer]; j++) {
	int iModule = module + j;
	double phi0 = phiMin + iModule * gridPhi;
	  
	verts3D[0] = rMin*cos(phi0) + (Lin+deltaL)*cos(phi0+alpha);
	verts3D[1] = rMin*sin(phi0) + (Lin+deltaL)*sin(phi0+alpha);
	verts3D[2] = r_in/tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
	  
	verts3D[3] = rMin*cos(phi0+gridPhi) + (Lin-deltaL)*cos(phi0+gridPhi+alpha);
	verts3D[4] = rMin*sin(phi0+gridPhi) + (Lin-deltaL)*sin(phi0+gridPhi+alpha);
	verts3D[5] = r_in/tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
	  
	verts3D[6] = rMin*cos(phi0+gridPhi) + (Lin-deltaL)*cos(phi0+gridPhi+alpha);
	verts3D[7] = rMin*sin(phi0+gridPhi) + (Lin-deltaL)*sin(phi0+gridPhi+alpha);
	verts3D[8] = r_in/tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
	  
	verts3D[9] = rMin*cos(phi0) + (Lin+deltaL)*cos(phi0+alpha);
	verts3D[10] = rMin*sin(phi0) + (Lin+deltaL)*sin(phi0+alpha);
	verts3D[11] = r_in/tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
	  
	verts3D[12] = rMin*cos(phi0) + (Lout+deltaL)*cos(phi0+alpha);
	verts3D[13] = rMin*sin(phi0) + (Lout+deltaL)*sin(phi0+alpha);
	verts3D[14] = r_out/tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
	  
	verts3D[15] = rMin*cos(phi0+gridPhi) + (Lout-deltaL)*cos(phi0+gridPhi+alpha);
	verts3D[16] = rMin*sin(phi0+gridPhi) + (Lout-deltaL)*sin(phi0+gridPhi+alpha);
	verts3D[17] = r_out/tan(theta_center - thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
	  
	verts3D[18] = rMin*cos(phi0+gridPhi) + (Lout-deltaL)*cos(phi0+gridPhi+alpha);
	verts3D[19] = rMin*sin(phi0+gridPhi) + (Lout-deltaL)*sin(phi0+gridPhi+alpha);
	verts3D[20] = r_out/tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
	  
	verts3D[21] = rMin*cos(phi0) + (Lout+deltaL)*cos(phi0+alpha);
	verts3D[22] = rMin*sin(phi0) + (Lout+deltaL)*sin(phi0+alpha);
	verts3D[23] = r_out/tan(theta_center + thetaGrid * geomReader->mergedCells_Theta[layer]/2.);
	  
	bs[icl]->AddBox(verts3D);
	bs[icl]->DigitValue( (int) (1000 * energy) );
	//bs->BoxId(new TNamed(Form("Cell %lu", cellID), "Dong!"));
      }
    }
  }

  if (eventLabel == nullptr) {
    eventLabel = new TGLConstAnnotation(gEve->GetDefaultGLViewer(),
					Form("%s, %.1f GeV\nEvent %d",
					     partType.Data(), pmax, eventId), 0.1, 0.9);
    eventLabel->SetTextSize(0.05);// % of window diagonal
    eventLabel->SetAllowClose(false);
  }
  else  {
    eventLabel->SetText(Form("%s, %.1f GeV\nEvent %d", partType.Data(), pmax, eventId));
  }
    
  TEveElement* top = (TEveElement*) gEve->GetCurrentEvent();
  // if nothing is drawn, top is null
  if (top) {
    rhoPhiEventScene->DestroyElements();
    rhoPhiProjManager->ImportElements(top, rhoPhiEventScene);

    // slow??
    //TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(pal, 0.55, 0.1, 0.4, 0.05);
    //rhoPhiGLView->AddOverlayElement(po);

    rhoZEventScene->DestroyElements();
    rhoZProjManager->ImportElements(top, rhoZEventScene);
      
    //TEveRGBAPaletteOverlay *po2 = new TEveRGBAPaletteOverlay(pal, 0.55, 0.1, 0.4, 0.05);
    //rhoZGLView->AddOverlayElement(po2);
  }
  
  cout << "Done" << endl << endl;

  textEntry->SetTextColor((Pixel_t)0x000000);
  textEntry->SetText(Form("Event %d loaded", eventId));

}

EventDisplay::EventDisplay() {
}

void EventDisplay::startDisplay(int initialEvent) {

  // calculate the geometry parameters
  geomReader = new DetectorGeometry;

  cout << "******************************************************************************" << endl;
  cout << "Displaying the geometry" << endl;
  cout << "******************************************************************************" << endl << endl;

  // create the eve manageer
  TEveManager::Create();

  // see palettes here: https://root.cern.ch/doc/master/classTColor.html
  // gStyle->SetPalette(kAvocado);
  gStyle->SetPalette(kSienna);
  
  // first tab
  gEve->GetDefaultGLViewer()->SetGuideState(TGLUtil::kAxesOrigin, false, false, 0);
  gEve->GetDefaultGLViewer()->DrawGuides();
  gEve->GetDefaultViewer()->SetElementName("3D view");
  gEve->GetDefaultGLViewer()->CurrentCamera().RotateRad(-.7, 0.5);

  // Create the geometry and the readout
  TEveElementList* geom = new TEveElementList("Geometry");
  TEveElementList* PCBs = new TEveElementList("PCBs");
  TEveElementList* actives = new TEveElementList("Active elements");
  TEveElementList* passives = new TEveElementList("Passive elements");
  TEveElementList* readout = new TEveElementList("Readout");
  if (useG4geom) {
    //auto fGeom = TFile::Open(geomFile.c_str(), "CACHEREAD");
    cout << "Reading Geant4 geometry from file " << geomFile << endl;
    auto fGeom = TFile::Open(geomFile.c_str(), "READ");
    if (!fGeom) return;
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) fGeom->Get("ECalBarrel_vol");
    barrel = TEveGeoShape::ImportShapeExtract(gse, 0);
    barrel->SetMainTransparency(70);
    fGeom->Close();
    delete fGeom;
    barrel->SetPickableRecursively(kTRUE);
    geom->AddElement(barrel);
    TPRegexp re;
    // set transparency of the subvolumes of the bath
    re = TPRegexp("LAr_bath*");
    TEveElement* bath = barrel->FindChild(re);
    TEveElement::List_t matches;
    re = TPRegexp("ECAL_Cryo*");
    barrel->FindChildren(matches, re);
    for (auto a : matches) a->SetMainTransparency(70);
    re = TPRegexp("services*");
    barrel->FindChildren(matches, re);
    for (auto a : matches) a->SetMainTransparency(70);
    // make lists of elements inside bath to turn on/off simultaneously
    if (bath) {
      TEveElementList* newbath = new TEveElementList("LAr_bath");
      barrel->AddElement(newbath);
      TEveElement::List_t matches;
      re = TPRegexp("PCB*");
      bath->FindChildren(matches, re);
      for (auto a : matches) PCBs->AddElement(a);
      newbath->AddElement(PCBs);
      TEveElement::List_t matches2;
      re = TPRegexp("active*");
      bath->FindChildren(matches2, re);
      for (auto a : matches2) actives->AddElement(a);
      newbath->AddElement(actives); 
      TEveElement::List_t matches3;     
      re = TPRegexp("passive*");
      bath->FindChildren(matches3, re);
      for (auto a : matches3) passives->AddElement(a);
      newbath->AddElement(passives);
      barrel->RemoveElement(bath);
      // hide elements inside bath by default because they are slow in 3D
      newbath->SetRnrSelfChildren(true, false);
    }
    gEve->AddGlobalElement(geom);
    gEve->AddToListTree(geom, true);
  }
  else {
    cout << "Creating simplified geometry based on calculated dimensions " << endl;
    // the barrel envelope
    barrel = new TEveGeoShape("Barrel");
    barrel->SetShape(new TGeoTube(geomReader->rMin, geomReader->rMax, geomReader->zMax));
    barrel->SetMainColor(kCyan);
    barrel->SetMainTransparency(90);
    barrel->SetNSegments(128);
    // geom->AddElement(barrel);
    gEve->AddGlobalElement(barrel);
    // the barrel layers
    TEveElementList* layers = new TEveElementList("layers");
    barrel->AddElement(layers);
    TEveGeoShape* b;
    for (int iLayer = 0; iLayer < geomReader->nLayers; iLayer++) {
      b = new TEveGeoShape(Form("Barrel layer %d", iLayer));
      b->SetShape(new TGeoTube(geomReader->r[iLayer], geomReader->r[iLayer+1], geomReader->zMax));
      b->SetMainColor(kCyan);
      b->SetMainTransparency(90);
      b->SetNSegments(128);
      layers->AddElement(b);
    }    
    // the electrodes
    TEveElementList* modules = new TEveElementList("modules");
    barrel->AddElement(modules);
    for (int iModule = 0; iModule < geomReader->nModules; iModule++) {
      double phi0 = iModule * geomReader->gridPhi - geomReader->gridPhi/12.; // small extra shift is due to finite width of element (?)
      double phi = phi0 + geomReader->dPhiAvg;
      b = new TEveGeoShape(Form("Module %d", iModule));
      b->SetShape(new TGeoBBox(geomReader->Ltot/2, 0.01, (geomReader->zMax - geomReader->zMin)/2.0));
      b->SetMainColor(kGray);
      b->SetMainTransparency(95);
      TGeoRotation* rot = new TGeoRotation();
      rot->SetAngles((geomReader->alpha + iModule * geomReader->gridPhi)*180./TMath::Pi(), 0., 0.);
      TGeoCombiTrans* c1 = new TGeoCombiTrans(geomReader->rAvg*cos(phi), geomReader->rAvg*sin(phi), 0.0, rot);
      b->SetTransMatrix(*c1);
      modules->AddElement(b);
    }
    gEve->AddToListTree(barrel, true);
  }

  gEve->AddToListTree(readout,true);

  // create second tab (R-phi view)
  rhoPhiView = gEve->SpawnNewViewer("Projection Rho-Phi");
  // two scenes, for geometry and event
  rhoPhiScene = gEve->SpawnNewScene("Rho-Phi geometry",
				    "Scene holding projected geometry data for the RhoPhi view.");
  rhoPhiView->AddScene(rhoPhiScene);
  if (evtFile!="") {
    rhoPhiEventScene = gEve->SpawnNewScene("RhoPhi Event Data",
					   "Scene holding projected event-data for the RhoPhi view.");
    rhoPhiView->AddScene(rhoPhiEventScene);
  }

  rhoPhiEventSceneManual = gEve->SpawnNewScene("RhoPhi Event Data 2",
					       "Scene holding hand-crafted event-data for the RhoPhi view.");
  rhoPhiView->AddScene(rhoPhiEventSceneManual);
  rhoPhiGLView = rhoPhiView->GetGLViewer();
  // set camera orientation
  rhoPhiGLView->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  // create 3D->2D projection manager for rho-phi
  rhoPhiProjManager = new TEveProjectionManager();
  rhoPhiProjManager->SetProjection(TEveProjection::kPT_RPhi);
  auto axes = new TEveProjectionAxes(rhoPhiProjManager);
  axes->SetElementName("Rho-Phi projection axes");
  rhoPhiScene->AddElement(axes);
  if (useG4geom)
    rhoPhiProjManager->ImportElements(geom, rhoPhiScene);
  else
    rhoPhiProjManager->ImportElements(barrel, rhoPhiScene);

  // the merged module grid
  TEveStraightLineSet* gridmod = new TEveStraightLineSet("phi readout merged");
  gridmod->SetLineColor(kViolet+2);
  gridmod->SetLineWidth(8);
  for (int iLayer=0; iLayer < geomReader->nLayers; iLayer++) {
    double Lin = geomReader->getL(geomReader->alpha, geomReader->rMin, geomReader->r[iLayer]);
    double Lout = geomReader->getL(geomReader->alpha, geomReader->rMin, geomReader->r[iLayer+1]);
    for (int iModule=0; iModule < geomReader->nModules; iModule++) {
      if (iModule % geomReader->mergedModules[iLayer] != 0) continue;
      double phi0 = geomReader->phiMin + iModule * geomReader->gridPhi;
      double x1 = geomReader->rMin*cos(phi0) + Lin * cos(phi0 + geomReader->alpha);
      double y1 = geomReader->rMin*sin(phi0) + Lin * sin(phi0 + geomReader->alpha);
      double x2 = geomReader->rMin*cos(phi0) + Lout * cos(phi0 + geomReader->alpha);
      double y2 = geomReader->rMin*sin(phi0) + Lout * sin(phi0 + geomReader->alpha);
      gridmod->AddLine(x1, y1, 0., x2, y2, 0.);
    }
  }
  // add to scene
  rhoPhiScene->AddElement(gridmod);
  readout->AddElement(gridmod);
  
  // third tab (R-z view)
  rhoZView = gEve->SpawnNewViewer("Projection Rho-Z");
  rhoZScene = gEve->SpawnNewScene("Rho-Z geometry",
				  "Scene holding projected geometry data for the RhoZ view.");
  rhoZView->AddScene(rhoZScene);
  if (evtFile!="") {
    rhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
					 "Scene holding projected event-data for the RhoZ view.");
    rhoZView->AddScene(rhoZEventScene);
  }
  rhoZEventSceneManual = gEve->SpawnNewScene("RhoZ Event Data 2",
					     "Scene holding hand-crafted event-data for the RhoZ view.");
  rhoZView->AddScene(rhoZEventSceneManual);
  rhoZGLView = rhoZView->GetGLViewer();
  rhoZGLView->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);

  rhoZProjManager = new TEveProjectionManager();
  rhoZProjManager->SetProjection(TEveProjection::kPT_RhoZ);
  auto axes2 = new TEveProjectionAxes(rhoZProjManager);
  axes2->SetElementName("Rho-Z projection axes");
  rhoZScene->AddElement(axes2);
  if (useG4geom)
    rhoZProjManager->ImportElements(geom, rhoZScene);
  else
    rhoZProjManager->ImportElements(barrel, rhoZScene);

  // the theta readout grid
  TEveStraightLineSet* grid = new TEveStraightLineSet("theta readout");
  grid->SetLineColor(kViolet);
  for (int iTheta=0; iTheta <= geomReader->nThetaBins; iTheta++) {
    double theta = geomReader->thetaMin + iTheta * geomReader->thetaGrid;
    double r1 = geomReader->rMin;
    double r2 = geomReader->rMax;
    double z1 = r1*cos(theta)/sin(theta);
    double z2 = r2*cos(theta)/sin(theta);
    if (z1 < geomReader->zMax && z1 > geomReader->zMin) {
      if (z2 > geomReader->zMax) {
	z2 = geomReader->zMax;
	r2 = z2*sin(theta)/cos(theta);
      }
      if (z2 < geomReader->zMin) {
	z2 = geomReader->zMin;
	r2 = z2*sin(theta)/cos(theta);
      }
      grid->AddLine(z1,  r1, 0., z2,  r2, 0.);
      grid->AddLine(z1, -r1, 0., z2, -r2, 0.);
    }
  }
  rhoZScene->AddElement(grid);
  readout->AddElement(grid);
  
  // the merged grid
  TEveStraightLineSet* grid2 = new TEveStraightLineSet("theta readout merged");
  grid2->SetLineColor(kViolet+2);
  grid2->SetLineWidth(8);
  for (int iLayer = 0; iLayer < geomReader->nLayers; iLayer++) {
    double r1 = geomReader->r[iLayer];
    double r2 = geomReader->r[iLayer+1];
    for (int iTheta = 0; iTheta <= geomReader->nThetaBins; iTheta++) {
      if (iTheta % geomReader->mergedCells_Theta[iLayer] != 0) continue;
      double theta = geomReader->thetaMin + iTheta * geomReader->thetaGrid;
      double z1 = r1*cos(theta)/sin(theta);
      double z2 = r2*cos(theta)/sin(theta);
      double r2tmp = r2;
      if (z1<geomReader->zMax && z1>geomReader->zMin) {
	if (z2 > geomReader->zMax) {
	  z2 = geomReader->zMax;
	  r2tmp = z2*sin(theta)/cos(theta);
	}
	if (z2 < geomReader->zMin) {
	  z2 = geomReader->zMin;
	  r2tmp = z2*sin(theta)/cos(theta);
	}
	grid2->AddLine(z1, r1, 0.,
		       z2, r2tmp, 0.);
	grid2->AddLine(z1, -r1, 0.,
		       z2, -r2tmp, 0.);
      }
    }
  }
  rhoZScene->AddElement(grid2);
  readout->AddElement(grid2);

  gEve->Redraw3D(true);

  if (evtFile!="") {

  //
  // create the gui for event navigation
  //
  makeGui();

  
  //
  // display the data
  //
  cout << endl;
  cout << "******************************************************************************" << endl;
  cout << "Setting up the event reader" << endl;
  cout << "******************************************************************************" << endl << endl;

  cout << "Reading event data from file " << evtFile << endl << endl;

  // setup the reader
  TFile* f = TFile::Open(evtFile.c_str(), "READ");
  eventReader = new EventReader(f);
  nEvents = eventReader->nEvents;
  
  // load and display the requested event
  cout << endl;
  cout << "******************************************************************************" << endl;
  cout << "Reading the events" << endl;
  cout << "******************************************************************************" << endl << endl;
  loadEvent(initialEvent);
  }
  
  
  // Set the 3D view as the active tab and rotate the camera
  gEve->GetBrowser()->GetTabRight()->SetTab(0);

  // Draw
  gEve->Redraw3D(true);
}


/******************************************************************************/
// GUI
/******************************************************************************/

void EventDisplay::fwd()
{
  if (eventId < nEvents - 1) {
    ++eventId;
    loadEvent(eventId);
  } else {
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText("Already at last event");
    printf("\nAlready at last event.\n");
  }
}

void EventDisplay::bck()
{
  if (eventId > 0) {
    --eventId;
    loadEvent(eventId);
  } else {
    textEntry->SetTextColor(0xff0000);
    textEntry->SetText("Already at first event");
    printf("\nAlready at first event.\n");
  }
}

void EventDisplay::makeGui()
{
  // Create minimal GUI for event navigation.

  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);

  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("GUI");
  frmMain->SetCleanup(kDeepCleanup);

  TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
  {
    TGPictureButton* b = 0;

    TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
    if (!std::filesystem::exists(icondir.Data())) {
      icondir = Form("%s/share/root/icons/", gSystem->Getenv("ROOTSYS"));
    }
    
    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 10, 2, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "bck()");
    
    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
    hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 10, 10, 10));
    b->Connect("Clicked()", "EventDisplay", this, "fwd()");
    
    textEntry = new TGTextEntry(hf);
    textEntry->SetEnabled(kFALSE);
    hf->AddFrame(textEntry, new TGLayoutHints(kLHintsLeft | kLHintsCenterY  |
					      kLHintsExpandX, 2, 10, 10, 10));
    
  }
  frmMain->AddFrame(hf);
  
  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
  
  browser->StopEmbedding();
  browser->SetTabTitle("Event Control", 0);
}
