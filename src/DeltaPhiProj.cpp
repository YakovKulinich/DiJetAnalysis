#include <TAxis.h>
#include "DeltaPhiProj.h"

DeltaPhiProj::DeltaPhiProj( int a0, int a1, int a2, int a3 ){
  m_vAxisI[0] = a0;
  m_vAxisI[1] = a1;
  m_vAxisI[2] = a2;
  m_vAxisI[3] = a3;
}

DeltaPhiProj::~DeltaPhiProj(){
  for( auto& ta : m_vTAxis ){ delete ta; }
}
