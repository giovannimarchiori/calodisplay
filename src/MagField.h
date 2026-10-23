/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class MagField: magnetic field of an ideal solenoid (uniform field inside,
// return field outside), for the propagation of the tracks
//
/******************************************************************************/

#ifndef __MAGFIELD_H__
#define __MAGFIELD_H__

#include "TEveTrackPropagator.h"

class MagField: public TEveMagField
{
private:
   bool m_magnetIsOn;
   bool m_reverse;
   double Bin;
   double Bout;
   double zmax;
   double Rmax;
   
public:
   MagField(double Bi, double Bo, double z, double R);
   ~MagField() override;
   void setMagnetState(bool state);
   bool isMagnetOn() const;
   void setReverseState(bool state);
   bool isReverse() const;
   TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t z) const override;
};

#endif
