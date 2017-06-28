#include <iostream>

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

std::vector<int> DeltaPhiProj::GetMappedBins( const std::vector< int >&  vAxisBins ){
  std::vector< int > vMappedaAxisBins( vAxisBins.size(), 0 );

  for( uint axis = 0; axis < vAxisBins.size(); axis++ ){
    int mappedAxis = GetAxisI( axis );
    vMappedaAxisBins[ mappedAxis ] = vAxisBins[ axis ];
  }

  return vMappedaAxisBins;
}

bool DeltaPhiProj::CorrectPhaseSpace( const std::vector< int >&  vAxisBins ){
  std::vector< int > vMappedaAxisBins = GetMappedBins( vAxisBins );
  /*
  std::cout << "Ystar1 bin " << vMappedaAxisBins[0] << std::endl;
  std::cout << "Ystar2 bin " << vMappedaAxisBins[1] << std::endl;
  std::cout << "Pt1    bin " << vMappedaAxisBins[2] << std::endl;
  std::cout << "Pt2    bin " << vMappedaAxisBins[3] << std::endl;
  */

  // rules for correct phase space
  // if ystar1 is not forward (bin 1) we are
  // not interested
  if( vMappedaAxisBins[0] != 1 ){ return false; }
  // if pt2 > pt1 we are not interested. By definition
  // leading jet has higher pt.
  if( vMappedaAxisBins[3] > vMappedaAxisBins[2] ){ return false; }
  
  return true;
}
