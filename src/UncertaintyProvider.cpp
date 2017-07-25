#include "UncertaintyTool.h"
#include <TLorentzVector.h>

UncertaintyTool::UncertaintyTool()
  : UncertaintyTool( 0 ){}

UncertaintyTool::UncertaintyTool( int uc )
  : m_uc( uc ), m_uncertaintyTool( NULL ){}  

UncertaintyTool::~UncertaintyTool(){}

void UncertaintyTool::Initialize(){
  int pos_uc = std::abs( m_uc );

  if( pos_uc == 1 ) {
    m_uncertaintyTool = new JERUncertaintyTool();
  } else if( pos_uc > 1 && pos_uc < 17 ){
    m_uncertaintyTool = new JESUncertaintyTool();
  }
}

void UncertaintyTool::ApplyUncertainty( TLorentzVector& jet, double factor ){
  m_uncertaintyTool( jet, factor );
}

void JERUncertaintyTool::ApplyUncertainty( TLorentzVector& jet, double factor ){

}

void JESUncertaintyTool::ApplyUncertainty( TLorentzVector& jet, double factor ){

}
