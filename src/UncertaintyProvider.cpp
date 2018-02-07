#include "UncertaintyProvider.h"
#include "MyRoot.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <iostream>

//--------------------------------
//        Uncertainty Tool
//--------------------------------
UncertaintyTool::UncertaintyTool( int uc, bool is_pPb )
  : m_uc( uc ), m_is_pPb( is_pPb ) {

  m_sign   = uc > 0 ? 1 : -1;
  m_system = m_is_pPb ? "pPb" : "pp";
  rand     = new TRandom3();
}

UncertaintyTool::~UncertaintyTool(){}

void UncertaintyTool::RegisterUFactors( std::vector< std::vector< float > >* p_vSysUncert )
{ m_p_vSysUncert = p_vSysUncert; }


double UncertaintyTool::GetUncertaintyWeight( const TLorentzVector& jet1,
					      const TLorentzVector& jet2 )
{ return 1; }

int UncertaintyTool::GetEtaUJERBin(float eta){

  int yBin=0;
  
  if (fabs(eta) < 0.3) yBin =0;
  if (fabs(eta) > 0.3 && fabs(eta) < 0.8) yBin =1;
  if (fabs(eta) > 0.8 && fabs(eta) < 1.2) yBin =2;
  if (fabs(eta) > 1.2 && fabs(eta) < 2.1) yBin =3;
  if (fabs(eta) > 2.1 && fabs(eta) < 2.8) yBin =4;
  if (fabs(eta) > 2.8 && fabs(eta) < 3.6) yBin =5;
  if (fabs(eta) > 3.6) yBin =6;

  return yBin;
}

double UncertaintyTool::GetYstar( const TLorentzVector& jet )
{ return m_is_pPb ? jet.Rapidity() + constants::BETAZ : jet.Rapidity(); }


//--------------------------------
//   Unfolding Uncertainty Tool
//--------------------------------

UnfoldingUncertaintyTool::UnfoldingUncertaintyTool( int uc, bool is_pPb )
  : UncertaintyTool( uc, is_pPb ){

  // get the pp MC dPhi
  // TFile* fMC = TFile::Open( "data/dPhi_reco_pp_mc_pythia8.root" );
}

UnfoldingUncertaintyTool::~UnfoldingUncertaintyTool(){}

void UnfoldingUncertaintyTool::ApplyUncertainties( std::vector< TLorentzVector >& recoJets,
						   std::vector< TLorentzVector >& truthJets ){}

double UnfoldingUncertaintyTool::GetUncertaintyWeight( const TLorentzVector& jet1,
						       const TLorentzVector& jet2 )
{ return 1; }

//--------------------------------
//      Angular Uncertainty Tool 
//--------------------------------

AngularUncertaintyTool::AngularUncertaintyTool( int uc, bool is_pPb )
  : UncertaintyTool( uc, is_pPb ){

  // get the pythia8 eta sigm+mean file 
  TFile* f_eta_p8 = new TFile( "data/recoTruthDeta_pp_mc_pythia8.root", "read" );
  hAngularUncertEta = static_cast< TH2D* >( f_eta_p8->Get( "h_recoTruthDeta_sigma" ) );
  hAngularUncertEta->SetName( "hAngularUncertEta" );
  hAngularUncertEta->SetDirectory(0);
  f_eta_p8->Close();
  delete f_eta_p8;
  
  // clone this before modifying it.
  hAngularResEta = static_cast< TH2D* >( hAngularUncertEta->Clone( "hAngularResEta" ) );

  // get the reference herwig eta sigm+mean file 
  TFile* f_eta_ref = new TFile( "data/recoTruthDeta_pp_mc_herwig.root", "read" );
  TH2D* hAngularUncertEtaTmp = static_cast< TH2D* >( f_eta_ref->Get( "h_recoTruthDeta_sigma" ) );
  hAngularUncertEtaTmp->SetDirectory(0);
  f_eta_ref->Close();
  delete f_eta_ref;

  // now subtract the pythia8 and the reference (herwig)
  hAngularUncertEta->Add( hAngularUncertEtaTmp, -1 );

  // get the pythia8 phi sigm+mean file 
  TFile* f_phi_p8 = new TFile( "data/recoTruthDphi_pp_mc_pythia8.root", "read" );
  hAngularUncertPhi = static_cast< TH2D* >( f_phi_p8->Get( "h_recoTruthDphi_sigma" ) );
  hAngularUncertPhi->SetName( "hAngularUncertPhi" );
  hAngularUncertPhi->SetDirectory(0);
  f_phi_p8->Close();
  delete f_phi_p8;

  // clone this before modifying it.
  hAngularResPhi = static_cast< TH2D* >( hAngularUncertPhi->Clone( "hAngularResPhi" ) );

    // get the reference herwig phi sigm+mean file 
  TFile* f_phi_ref = new TFile( "data/recoTruthDphi_pp_mc_herwig.root", "read" );
  TH2D* hAngularUncertPhiTmp = static_cast< TH2D* >( f_phi_ref->Get( "h_recoTruthDphi_sigma" ) );
  hAngularUncertPhiTmp->SetDirectory(0);
  f_phi_ref->Close();
  delete f_phi_ref;
  
  // now subtract the pythia8 and the reference (herwig)
  hAngularUncertPhi->Add( hAngularUncertPhiTmp, -1 );

}

AngularUncertaintyTool::~AngularUncertaintyTool(){

  delete hAngularUncertEta;
  delete hAngularUncertPhi;
  delete hAngularResEta;
  delete hAngularResPhi;
}

void AngularUncertaintyTool::ApplyUncertainties( std::vector< TLorentzVector >& recoJets,
						 std::vector< TLorentzVector >& truthJets ){

  for( uint iJet = 0; iJet < recoJets.size(); iJet++ ){
    TLorentzVector& recoJet  = recoJets [ iJet ];
    TLorentzVector& truthJet = truthJets[ iJet ];

    float recoJetPt    = recoJet.Pt()/1000.;
    float recoJetEta   = recoJet.Eta();
    float recoJetYstar = GetYstar( recoJet );
    float recoJetPhi   = recoJet.Phi();
    float recoJetM     = recoJet.M();
  
    float truthJetEta   = truthJet.Eta();
    float truthJetPhi   = truthJet.Phi();
 
    float uncertaintyEta  = hAngularUncertEta->
      GetBinContent( hAngularUncertEta->FindBin( recoJetYstar, recoJetPt ) );
    float uncertaintyPhi  = hAngularUncertPhi->
      GetBinContent( hAngularUncertPhi->FindBin( recoJetYstar, recoJetPt ) );

    float etaRes = hAngularResEta->
      GetBinContent( hAngularResEta->FindBin( recoJetYstar, recoJetPt ) );
    float phiRes = hAngularResPhi->
      GetBinContent( hAngularResPhi->FindBin( recoJetYstar, recoJetPt ) );
 
    float smearingFactorSystEta =
      sqrt( pow( etaRes + uncertaintyEta, 2 ) - pow( etaRes, 2 ) );
    float smearingFactorSystPhi =
      sqrt( pow( phiRes + uncertaintyPhi, 2 ) - pow( phiRes, 2 ) );

    float correctionEta = rand->Gaus(0., smearingFactorSystEta);
    float correctionPhi = rand->Gaus(0., smearingFactorSystPhi);
          
    float recoJetEtaNew = recoJetEta + truthJetEta * correctionEta;
    float recoJetPhiNew = recoJetPhi + truthJetPhi * correctionPhi;

    /*
    std::cout << "----" << std::endl;
    std::cout << uncertaintyEta << " ... " << etaRes << " ... "
	      << correctionEta << " : " << recoJetEta << " -> " << recoJetEtaNew << std::endl;
    std::cout << uncertaintyPhi << " ... " << phiRes << " ... "
	      << correctionPhi << " : " << recoJetPhi << " -> " << recoJetPhiNew << std::endl;
    */
    
    recoJet.SetPtEtaPhiM
      ( recoJetPt * 1000, recoJetEtaNew, recoJetPhiNew, recoJetM );
  }
}

//--------------------------------
//      JER Uncertainty Tool 
//--------------------------------

JERUncertaintyTool::JERUncertaintyTool( int uc, bool is_pPb )
  : UncertaintyTool( uc, is_pPb ){

  // get the analysis jer 
  TFile* f_JER = new TFile
    ( Form( "data/recoTruthRpt_%s_mc_pythia8.root", m_system.c_str() ), "read" );

  hJER = static_cast<TH2D*>( f_JER->Get( "h_recoTruthRpt_sigma" ) );
  hJER->SetDirectory(0);
  f_JER->Close();

  // get the hi provided jer uncertainty factors
  TFile* f_HI_JER = new TFile(  "data/hi_jer.root" ,"read" );

  for( int ybin = 0; ybin < 7; ybin++ ){
    m_vJERhistos.push_back( (TH1D*)f_HI_JER->Get
			    (Form( "delta_sigma_hi_eta%i" ,ybin )));
    m_vJERhistos.back()->SetDirectory(0);
  }

  f_HI_JER->Close();
}

JERUncertaintyTool::~JERUncertaintyTool(){
  
  delete hJER;
  for( auto& h : m_vJERhistos ){ delete h; }
}

void JERUncertaintyTool::ApplyUncertainties( std::vector< TLorentzVector >& recoJets,
					     std::vector< TLorentzVector >& truthJets ){

  for( uint iJet = 0; iJet < recoJets.size(); iJet++ ){
    TLorentzVector& recoJet  = recoJets [ iJet ];
    TLorentzVector& truthJet = truthJets[ iJet ];

  
    float recoJetPt    = recoJet.Pt()/1000.;
    float recoJetEta   = recoJet.Eta();
    float recoJetYstar = GetYstar( recoJet );
    float recoJetPhi   = recoJet.Phi();
    float recoJetM     = recoJet.M();
  
    float truthJetPt   = truthJet.Pt()/1000.;

    int   etaBin       = GetEtaUJERBin( recoJetEta );
    float uncertainty  = m_vJERhistos[ etaBin ]->Interpolate( recoJetPt );

    // float JER       = hJER->Interpolate( recoJetYstar, recoJetPt );
    int   jerBin = hJER->FindBin( recoJetYstar, recoJetPt );
    
    float JER          = hJER->GetBinContent( jerBin ); 
    float smearingFactorSyst =
      std::sqrt( std::pow( JER + uncertainty, 2 ) - pow( JER, 2 ) );

    float correction   = rand->Gaus(0., smearingFactorSyst);
          
    float recoJetPtNew = recoJetPt + truthJetPt * correction;
    if ( recoJetPtNew < 10){
      recoJetPtNew = 1;
      recoJetM = 0;
    }
  
    recoJet.SetPtEtaPhiM
      ( recoJetPtNew * 1000., recoJetEta, recoJetPhi, recoJetM );
  }
}

//--------------------------------
//      HIJES Uncertainty Tool 
//--------------------------------

HIJESUncertaintyTool::HIJESUncertaintyTool( int uc, bool is_pPb )
  : UncertaintyTool( uc, is_pPb ){

  // load histograms for JES
  //HI JES <-> crosscalibration
  TFile* f_HI_JES = new TFile(  "data/cc_sys_090816.root" ,"read" );
  for( int ybin = 0; ybin < 7; ybin++ ){
    m_vJEShistos.push_back( (TH1D*)f_HI_JES->Get(Form("fsys_rel_%i",ybin)));
    m_vJEShistos.back()->SetDirectory(0);
  }
  
  f_HI_JES->Close();
}

HIJESUncertaintyTool::~HIJESUncertaintyTool(){
  for( auto& h : m_vJEShistos ){ delete h; }
}

void HIJESUncertaintyTool::ApplyUncertainties( std::vector< TLorentzVector >& recoJets,
					       std::vector< TLorentzVector >& truthJets ){

  for( uint iJet = 0; iJet < recoJets.size(); iJet++ ){
    TLorentzVector& recoJet  = recoJets[ iJet ];
  
    double eta = recoJet.Eta();
    double pT  = recoJet.Pt()/1000.;
  
    int etaBin = GetEtaUJERBin( eta );
    // need to subtract 1 after, the values are around 1.
    double histoFactor = m_vJEShistos[ etaBin ]->Interpolate( pT );
    histoFactor -= 1;
    
    recoJet.SetPtEtaPhiM( recoJet.Pt() * ( 1 + histoFactor * m_sign ),
			  recoJet.Eta(), recoJet.Phi(), recoJet.M() );

  }
}

//--------------------------------
//      JES Uncertainty Tool
//--------------------------------

JESUncertaintyTool::JESUncertaintyTool( int uc, bool is_pPb )
  : UncertaintyTool( uc, is_pPb ){}

JESUncertaintyTool::~JESUncertaintyTool(){}

void JESUncertaintyTool::ApplyUncertainties( std::vector< TLorentzVector >& recoJets,
					     std::vector< TLorentzVector >& truthJets ){

  for( uint iJet = 0; iJet < recoJets.size(); iJet++ ){
    TLorentzVector& recoJet = recoJets[ iJet ];
    double factor = (*m_p_vSysUncert)[ iJet ][ std::abs( m_uc ) - 1 ];
      
    recoJet.SetPtEtaPhiM( recoJet.Pt() * ( 1 + factor * m_sign ),
			  recoJet.Eta(), recoJet.Phi(), recoJet.M() );
  }
}

//--------------------------------
//      Uncertainty Provider  
//--------------------------------

UncertaintyProvider::UncertaintyProvider() : UncertaintyProvider( 0, false ){}

UncertaintyProvider::UncertaintyProvider( int uc, bool is_pPb ) : m_uncertaintyTool( NULL ){
  
  int pos_uc = std::abs( uc );

  if( pos_uc > 0  && pos_uc <= 18 ){
    m_uncertaintyTool = new JESUncertaintyTool      ( uc, is_pPb );
  } else if( pos_uc == 19 ) {
    m_uncertaintyTool = new HIJESUncertaintyTool    ( uc, is_pPb );
  } else if( pos_uc == 20 ) {
    m_uncertaintyTool = new JERUncertaintyTool      ( uc, is_pPb );
  } else if( pos_uc == 21 ) {
    m_uncertaintyTool = new AngularUncertaintyTool  ( uc, is_pPb );
  } else if( pos_uc == 22 ) {
    m_uncertaintyTool = new UnfoldingUncertaintyTool( uc, is_pPb );
  } 
}  

UncertaintyProvider::~UncertaintyProvider(){
  delete m_uncertaintyTool; m_uncertaintyTool = NULL;
}

void UncertaintyProvider::RegisterUFactors( std::vector< std::vector< float > >* p_vSysUncert )
{ m_uncertaintyTool->RegisterUFactors( p_vSysUncert ); }

void UncertaintyProvider::ApplyUncertainties( std::vector< TLorentzVector >& recoJets,
					      std::vector< TLorentzVector >& truthJets ){

  m_uncertaintyTool->ApplyUncertainties( recoJets, truthJets );
}

double UncertaintyProvider::GetUncertaintyWeight( const TLorentzVector& jet1,
						  const TLorentzVector& jet2 )
{

  return m_uncertaintyTool->GetUncertaintyWeight( jet1, jet2 );
}
