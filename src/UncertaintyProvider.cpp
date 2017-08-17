#include "UncertaintyProvider.h"
#include "MyRoot.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <iostream>

//--------------------------------
//        Uncertainty Tool
//--------------------------------
UncertaintyTool::UncertaintyTool( int uc, bool is_pPb )
  : m_uc( uc ), m_is_pPb( is_pPb ) {
  m_sign = uc > 0 ? 1 : -1;
}

UncertaintyTool::~UncertaintyTool(){}

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
//      JER Uncertainty Tool 
//--------------------------------

JERUncertaintyTool::JERUncertaintyTool( int uc, bool is_pPb )
  : UncertaintyTool( uc, is_pPb ){
  std::string system = m_is_pPb ? "pPb" : "pp";

  // get the analysis jer 
  TFile* f_JER = new TFile
    ( Form( "data/recoTruthRpt_%s_mc_pythia8.root", system.c_str() ), "read" );

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

  rand = new TRandom3();
}

JERUncertaintyTool::~JERUncertaintyTool(){}

void JERUncertaintyTool::ApplyUncertainty( TLorentzVector& recoJet,
					   TLorentzVector& truthJet,
					   double factor ){

  float recoJetPt    = recoJet.Pt()/1000.;
  float recoJetEta   = recoJet.Eta();
  float recoJetYstar = GetYstar( recoJet );
  float recoJetPhi   = recoJet.Phi();
  float recoJetM     = recoJet.M();
  
  float truthJetPt   = truthJet.Pt()/1000.;

  int   etaBin       = GetEtaUJERBin( recoJetEta );
  float uncertainty  = m_vJERhistos[ etaBin ]->Interpolate( recoJetPt );

  float JER          = hJER->Interpolate( recoJetYstar, recoJetPt ); 
  float smearingFactorSyst = sqrt(pow(JER+uncertainty,2)-pow(JER,2));

  float correction   = rand->Gaus(0., smearingFactorSyst);
          
  float recoJetPtNew = recoJetPt + truthJetPt * correction;
  if ( recoJetPtNew < 10){
    recoJetPtNew = 1;
    recoJetM = 0;
  }
  
  recoJet.SetPtEtaPhiM( recoJetPtNew * 1000., recoJetEta,
			recoJetPhi  , recoJetM );
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

HIJESUncertaintyTool::~HIJESUncertaintyTool(){}

void HIJESUncertaintyTool::ApplyUncertainty( TLorentzVector& recoJet,
					     TLorentzVector& truthJet,
					     double factor){
  
  double eta = recoJet.Eta();
  double pT  = recoJet.Pt()/1000.;
  
  int etaBin = GetEtaUJERBin( eta );
  // need to subtract 1 after, the values are around 1.
  double histoFactor = m_vJEShistos[ etaBin ]->Interpolate( pT );
  histoFactor -= 1;
    
  recoJet.SetPtEtaPhiM( recoJet.Pt() * ( 1 + histoFactor * m_sign ),
			recoJet.Eta(), recoJet.Phi(), recoJet.M() );
}

//--------------------------------
//      JES Uncertainty Tool
//--------------------------------

JESUncertaintyTool::JESUncertaintyTool( int uc, bool is_pPb )
  : UncertaintyTool( uc, is_pPb ){}

JESUncertaintyTool::~JESUncertaintyTool(){}

void JESUncertaintyTool::ApplyUncertainty( TLorentzVector& recoJet,
					   TLorentzVector& truthJet,
					   double factor){
  recoJet.SetPtEtaPhiM( recoJet.Pt() * ( 1 + factor * m_sign ),
			recoJet.Eta(), recoJet.Phi(), recoJet.M() );
}

//--------------------------------
//      Uncertainty Provider  
//--------------------------------

UncertaintyProvider::UncertaintyProvider() : UncertaintyProvider( 0, false ){}

UncertaintyProvider::UncertaintyProvider( int uc, bool is_pPb ) : m_uncertaintyTool( NULL ){
  
  int pos_uc = std::abs( uc );

  if( pos_uc > 0  && pos_uc <= 18 ){
    m_uncertaintyTool = new JESUncertaintyTool( uc, is_pPb );
  } else if( pos_uc == 19 ) {
    m_uncertaintyTool = new HIJESUncertaintyTool( uc, is_pPb );
  } else if( pos_uc == 20 ) {
    m_uncertaintyTool = new JERUncertaintyTool( uc, is_pPb );
  } 
}  

UncertaintyProvider::~UncertaintyProvider(){
  delete m_uncertaintyTool; m_uncertaintyTool = NULL;
}

void UncertaintyProvider::ApplyUncertainty( TLorentzVector& recoJet,
					    TLorentzVector& truthJet,
					    double factor ){

  m_uncertaintyTool->ApplyUncertainty( recoJet, truthJet, factor );
}
