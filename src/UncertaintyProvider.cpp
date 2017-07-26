#include "UncertaintyProvider.h"

#include <TLorentzVector.h>

//--------------------------------
//        Uncertainty Tool
//--------------------------------
UncertaintyTool::UncertaintyTool( int uc ) : m_uc( uc ){
  m_sign = uc > 0 ? 1 : -1;
}

UncertaintyTool::~UncertaintyTool(){}

//--------------------------------
//      JER Uncertainty Tool 
//--------------------------------

JERUncertaintyTool::JERUncertaintyTool( int uc ) : UncertaintyTool( uc ){}

JERUncertaintyTool::~JERUncertaintyTool(){}

void JERUncertaintyTool::ApplyUncertainty( TLorentzVector& jet, double factor )
{}

//--------------------------------
//      JES Uncertainty Tool
//--------------------------------

JESUncertaintyTool::JESUncertaintyTool( int uc )
  : UncertaintyTool( uc ){}

JESUncertaintyTool::~JESUncertaintyTool(){}

void JESUncertaintyTool::ApplyUncertainty( TLorentzVector& jet, double factor ){
  jet.SetPtEtaPhiM( jet.Pt()*( 1 + factor * m_sign ),
		    jet.Eta(), jet.Phi(), jet.M() );
}

//--------------------------------
//      Uncertainty Provider  
//--------------------------------

UncertaintyProvider::UncertaintyProvider() : UncertaintyProvider( 0 ){}

UncertaintyProvider::UncertaintyProvider( int uc ) : m_uncertaintyTool( NULL ){
  
  int pos_uc = std::abs( uc );

  if( pos_uc == 1 ) {
    m_uncertaintyTool = new JERUncertaintyTool( uc );
  } else if( pos_uc > 1 && pos_uc < 17 ){
    m_uncertaintyTool = new JESUncertaintyTool( uc );
  }
}  

UncertaintyProvider::~UncertaintyProvider(){
  delete m_uncertaintyTool; m_uncertaintyTool = NULL;
}

void UncertaintyProvider::ApplyUncertainty( TLorentzVector& jet, double factor ){
  m_uncertaintyTool->ApplyUncertainty( jet, factor );
}
