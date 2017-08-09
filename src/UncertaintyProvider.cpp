#include "UncertaintyProvider.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <iostream>

//--------------------------------
//        Uncertainty Tool
//--------------------------------
UncertaintyTool::UncertaintyTool( int uc ) : m_uc( uc ){
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


//--------------------------------
//      HIJER Uncertainty Tool 
//--------------------------------

HIJERUncertaintyTool::HIJERUncertaintyTool( int uc ) : UncertaintyTool( uc ){
  // load histograms for JES
  //HI JES <-> crosscalibration
  TFile* f_HI_JES = new TFile(  "data/cc_sys_090816.root" ,"read" );
  for( int ybin = 0; ybin < 7; ybin++ ){
    m_vJEShistos.push_back( (TH1D*)f_HI_JES->Get(Form("fsys_rel_%i",ybin)));
    m_vJEShistos.back()->SetDirectory(0);
  }
  f_HI_JES->Close();
}

HIJERUncertaintyTool::~HIJERUncertaintyTool(){}

void HIJERUncertaintyTool::ApplyUncertainty( TLorentzVector& jet, double factor )
{
  double eta = jet.Eta();
  double pT  = jet.Pt()/1000.;
  
  int etaBin = GetEtaUJERBin( eta );
  // need to subtract 1 after, the values are around 1.
  double histoFactor = m_vJEShistos[ etaBin ]->Interpolate( pT );
  histoFactor -= 1;
    
  jet.SetPtEtaPhiM( jet.Pt() * ( 1 + histoFactor * m_sign ),
		    jet.Eta(), jet.Phi(), jet.M() );
}


//--------------------------------
//      JER Uncertainty Tool 
//--------------------------------

JERUncertaintyTool::JERUncertaintyTool( int uc ) : UncertaintyTool( uc ){
  TFile* f_JER = new TFile( "data/recoTruthRpt_pp_mc_pythia8.root", "read" );
  hJER = static_cast<TH2D*>( f_JER->Get("h_recoTruthRpt_sigma") );
  hJER->SetDirectory(0);
  f_JER->Close();

  TFile* f_HI_JER = new TFile(  "data/hi_jer.root" ,"read" );
  for( int ybin = 0; ybin < 7; ybin++ ){
    m_vJERhistos.push_back( (TH1D*)f_HI_JER->Get(Form("delta_sigma_hi_eta%i",ybin)));
    m_vJERhistos.back()->SetDirectory(0);
  }
  f_HI_JER->Close();

  rand = new TRandom3();
}

JERUncertaintyTool::~JERUncertaintyTool(){}

void JERUncertaintyTool::ApplyUncertainty( TLorentzVector& jet, double factor )
{}

void JERUncertaintyTool::ApplyUncertainty( TLorentzVector& recoJet,
					   TLorentzVector& truthJet,
					   double factor ){

  float jetPtRecon   = recoJet.Pt()/1000.;
  float jetEtaRecon  = recoJet.Eta();
  float jetPhiRecon  = recoJet.Phi();
  float jetMRecon    = recoJet.M();
  
  float jetPtTruth   = truthJet.Pt()/1000.;

  int   etaBin       = GetEtaUJERBin( jetEtaRecon );
  float uncertainty  = m_vJERhistos[ etaBin ]->Interpolate( jetPtRecon );

  float JER          = hJER->Interpolate( jetEtaRecon, jetPtRecon ); 
  float smearingFactorSyst = sqrt(pow(JER+uncertainty,2)-pow(JER,2));

  float correction   = rand->Gaus(0., smearingFactorSyst);
          
  float jetPtRecoNew = jetPtRecon + jetPtTruth * correction;
  if ( jetPtRecoNew < 10){
    jetPtRecoNew = 1;
    jetMRecon = 0;
  }

  recoJet.SetPtEtaPhiM( jetPtRecoNew, jetEtaRecon,
			jetPhiRecon , jetMRecon );
}

//--------------------------------
//      JES Uncertainty Tool
//--------------------------------

JESUncertaintyTool::JESUncertaintyTool( int uc )
  : UncertaintyTool( uc ){}

JESUncertaintyTool::~JESUncertaintyTool(){}

void JESUncertaintyTool::ApplyUncertainty( TLorentzVector& jet, double factor ){
  jet.SetPtEtaPhiM( jet.Pt() * ( 1 + factor * m_sign ),
		    jet.Eta(), jet.Phi(), jet.M() );
}

//--------------------------------
//      Uncertainty Provider  
//--------------------------------

UncertaintyProvider::UncertaintyProvider() : UncertaintyProvider( 0 ){}

UncertaintyProvider::UncertaintyProvider( int uc ) : m_uncertaintyTool( NULL ){
  
  int pos_uc = std::abs( uc );

  if( pos_uc > 0  && pos_uc <= 18 ){
    m_uncertaintyTool = new JESUncertaintyTool( uc );
  } else if( pos_uc == 19 ) {
    m_uncertaintyTool = new JERUncertaintyTool( uc );
  } else if( pos_uc == 20 ) {
    m_uncertaintyTool = new HIJERUncertaintyTool( uc );
  } 
}  

UncertaintyProvider::~UncertaintyProvider(){
  delete m_uncertaintyTool; m_uncertaintyTool = NULL;
}

void UncertaintyProvider::ApplyUncertainty( TLorentzVector& jet, double factor ){
  m_uncertaintyTool->ApplyUncertainty( jet, factor );
}
