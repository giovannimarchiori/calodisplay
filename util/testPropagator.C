#include "TEveTrackPropagator.h"
#include "TEveVector.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>

void PropagateExactPoint() {
    // ==========================================
    // 1. INITIAL CONDITIONS & PHYSICAL CONSTANTS
    // ==========================================
    TEveVectorD start_pos(0.0, 0.0, 0.0);       // Origin (cm)
    TEveVectorD init_mom(1.0, 0.5, 2.0);        // Initial momentum (GeV/c)
    Int_t charge = 1;                           // Positron/Proton (+1)
    Double_t Bz = 2.0;                          // Constant field of 2.0 Tesla along Z
    
    // Constant conversion factor for GeV/c, Tesla, and cm
    Double_t k = 0.00299792458;
    Double_t pT = TMath::Sqrt(init_mom.fX * init_mom.fX + init_mom.fY * init_mom.fY);
    Double_t R = pT / (k * Bz * charge);        // Helix radius in cm

    // ==========================================
    // 2. ANALYTICALLY GENERATE A POINT ON THE TRAJECTORY
    // ==========================================
    // Let's find where the particle is after rotating by exactly 45 degrees (pi/4 radians)
    Double_t delta_phi = TMath::Pi() / 4.0; 
    
    // Calculate the Z displacement for this rotation: delta_z = delta_phi * p_z / (k * B * q)
    Double_t delta_z = (delta_phi * init_mom.fZ) / (k * Bz * charge);
    
    // Center of the helix circle in the XY plane
    Double_t xc = start_pos.fX - (init_mom.fY / (k * Bz * charge));
    Double_t yc = start_pos.fY + (init_mom.fX / (k * Bz * charge));
    
    // Initial phase angle in XY plane
    Double_t phi_0 = TMath::ATan2(start_pos.fY - yc, start_pos.fX - xc);
    
    // Exact target position on the helix after delta_phi rotation
    TEveVectorD target_pos;
    target_pos.fX = xc + R * TMath::Cos(phi_0 + delta_phi);
    target_pos.fY = yc + R * TMath::Sin(phi_0 + delta_phi);
    target_pos.fZ = start_pos.fZ + delta_z;

    // Exact analytical momentum at this target point
    Double_t analytical_px = init_mom.fX * TMath::Cos(delta_phi) - init_mom.fY * TMath::Sin(delta_phi);
    Double_t analytical_py = init_mom.fX * TMath::Sin(delta_phi) + init_mom.fY * TMath::Cos(delta_phi);
    Double_t analytical_pz = init_mom.fZ;

    // ==========================================
    // 3. EXTRAPOLATION VIA ROOT (TEveTrackPropagator)
    // ==========================================
    TEveTrackPropagator* propagator = new TEveTrackPropagator();
    propagator->SetMagField(0.0, 0.0, Bz);
    // propagator->SetStepper(TEveTrackPropagator::kHelix);
    
    propagator->InitTrack(start_pos, charge);
    
    TEveVectorD propagated_mom = init_mom; // Will be modified in-place
    Bool_t success = propagator->GoToVertex(target_pos, propagated_mom);

    // ==========================================
    // 4. DISPLAY COMPARISON
    // ==========================================
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "========================================================\n";
    std::cout << " TARGET GENERATED ON THE HELIX PATH\n";
    std::cout << " Target Coordinates (cm): (" << target_pos.fX << ", " << target_pos.fY << ", " << target_pos.fZ << ")\n";
    std::cout << "========================================================\n";
    
    if (!success) {
        std::cout << "ROOT Propagator failed to reach the designated vertex!\n";
    }

    std::cout << "Component | ROOT Extrapolation | Analytical Formula | Difference \n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << " Px       | " << std::setw(18) << propagated_mom.fX 
              << " | " << std::setw(18) << analytical_px 
              << " | " << std::setw(10) << (propagated_mom.fX - analytical_px) << "\n";
              
    std::cout << " Py       | " << std::setw(18) << propagated_mom.fY 
              << " | " << std::setw(18) << analytical_py 
              << " | " << std::setw(10) << (propagated_mom.fY - analytical_py) << "\n";
              
    std::cout << " Pz       | " << std::setw(18) << propagated_mom.fZ 
              << " | " << std::setw(18) << analytical_pz 
              << " | " << std::setw(10) << (propagated_mom.fZ - analytical_pz) << "\n";
    std::cout << "========================================================\n";
    
    delete propagator;
}
