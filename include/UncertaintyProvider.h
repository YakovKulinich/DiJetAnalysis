#ifndef UNCERTAINTYPROVIDER_H
#define UNCERTAINTYPROVIDER_H

#include <vector>

class TH1;
class TH2;
class TRandom3;
class TLorentzVector;

//--------------------------------
//        Uncertainty Tool
//--------------------------------
class UncertaintyTool{

 public:
  UncertaintyTool( int, bool );
  virtual ~UncertaintyTool();
  
  virtual void ApplyUncertainty( TLorentzVector&, TLorentzVector&, double ) {}

  int GetEtaUJERBin( float eta );

  double GetYstar( const TLorentzVector& jet );
  
 protected:
  int  m_uc;
  bool m_is_pPb;
  
  int m_sign;
};

//--------------------------------
//      Angular Uncertainty Tool 
//--------------------------------
class AngularUncertaintyTool : public UncertaintyTool{
 public:
  AngularUncertaintyTool( int, bool );
  ~AngularUncertaintyTool();

  void ApplyUncertainty( TLorentzVector&, TLorentzVector&, double = 0 );

 private:
  TRandom3* rand;
  
  TH2* hAngularUncertEta;
  TH2* hAngularUncertPhi;

  TH2* hAngularResEta;
  TH2* hAngularResPhi;
};

//--------------------------------
//      JER Uncertainty Tool 
//--------------------------------
class JERUncertaintyTool : public UncertaintyTool{
 public:
  JERUncertaintyTool( int, bool );
  ~JERUncertaintyTool();

  void ApplyUncertainty( TLorentzVector&, TLorentzVector&, double = 0 );

 private:
  TRandom3* rand;
  
  TH2* hJER;
  std::vector<TH1*> m_vJERhistos;
};

//--------------------------------
//      HIJES Uncertainty Tool 
//--------------------------------
class HIJESUncertaintyTool : public UncertaintyTool{
 public:
  HIJESUncertaintyTool( int, bool );
  ~HIJESUncertaintyTool();

  void ApplyUncertainty( TLorentzVector&, TLorentzVector&, double = 0 );

 private:
  std::vector<TH1*> m_vJEShistos;
};

//--------------------------------
//      JES Uncertainty Tool
//--------------------------------
class JESUncertaintyTool : public UncertaintyTool{
 public:
  JESUncertaintyTool( int, bool );
  ~JESUncertaintyTool();
  
  void ApplyUncertainty( TLorentzVector&, TLorentzVector&, double = 0 );
};

//--------------------------------
//      Uncertainty Provider  
//--------------------------------
class UncertaintyProvider{

 public:
  UncertaintyProvider();
  UncertaintyProvider( int, bool );
  virtual ~UncertaintyProvider();

  void ApplyUncertainty( TLorentzVector&, TLorentzVector&, double = 0 );
  
 private:
  UncertaintyTool*  m_uncertaintyTool;
};


#endif
