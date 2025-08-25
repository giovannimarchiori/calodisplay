#include "MagField.h"
#include <iostream>

MagField::MagField(double Bi, double Bo, double z, double R):
  m_magnetIsOn(true),
  Bin(Bi),
  Bout(Bo),
  zmax(z),
  Rmax(R)
{};

MagField::~MagField()
{};

void MagField::setMagnetState( bool state )
{
  if (state != m_magnetIsOn)
    {
      if ( state )
	std::cout << "Magnet state is changed to ON" << std::endl;
      else
	std::cout << "Magnet state is changed to OFF" << std::endl;
    }
  m_magnetIsOn = state;
}

bool MagField::isMagnetOn() const
{
  return m_magnetIsOn;
}

void MagField::setReverseState(bool state)
{
  m_reverse = state; }

bool MagField::isReverse() const
{
  return m_reverse;

}

TEveVectorD MagField::GetFieldD(Double_t x, Double_t y, Double_t z) const
{
  if (!m_magnetIsOn)
    return TEveVectorD(0,0,0);

  if ( TMath::Abs(z) > zmax) {
    return TEveVectorD(0,0,0);
  }

  double R = sqrt(x*x+y*y);
  double field = (R<Rmax) ? Bin : Bout;
  if (m_reverse) field = -field;
  return TEveVectorD(0,0,field);
};

