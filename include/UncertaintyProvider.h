#ifndef UNCERTAINTYPROVIDER_H
#define UNCERTAINTYPROVIDER_H

#include <string>
#include <vector>

class TH1;
class TH2;
class TH3;
class THnSparse;
class TRandom3;
class TLorentzVector;
class HIJESUncertaintyProvider;

//--------------------------------
//        Uncertainty Tool
//--------------------------------
class UncertaintyTool{

 public:
  UncertaintyTool( int, bool );
  virtual ~UncertaintyTool();

  void RegisterUFactors( std::vector< std::vector< float > >* );
  
  virtual void ApplyUncertainties( std::vector< TLorentzVector >&,
				   std::vector< TLorentzVector >& ) = 0;

  int GetEtaUJERBin( float eta );

  double GetYstar( const TLorentzVector& jet );
  
 protected:
  int  m_uc;
  bool m_is_pPb;

  std::vector< std::vector< float > >* m_p_vSysUncert;
  
  int m_sign;
  std::string m_system;

  TRandom3* rand;
};

//--------------------------------
//      Angular Uncertainty Tool 
//--------------------------------
class AngularUncertaintyTool : public UncertaintyTool{
 public:
  AngularUncertaintyTool( int, bool );
  ~AngularUncertaintyTool();

  void ApplyUncertainties( std::vector< TLorentzVector >&,
			   std::vector< TLorentzVector >& );
  
 private:
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

  void ApplyUncertainties( std::vector< TLorentzVector >&,
			   std::vector< TLorentzVector >& );

 private:
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

  void ApplyUncertainties( std::vector< TLorentzVector >&,
			   std::vector< TLorentzVector >& );

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
  
  void ApplyUncertainties( std::vector< TLorentzVector >&,
			   std::vector< TLorentzVector >& );

 private:
   HIJESUncertaintyProvider*  m_hiJetUncertaintyTool;
};

//--------------------------------
//      Uncertainty Provider  
//--------------------------------
class UncertaintyProvider{

 public:
  UncertaintyProvider();
  UncertaintyProvider( int, bool );
  virtual ~UncertaintyProvider();

  void RegisterUFactors( std::vector< std::vector< float > >* );
  
  void ApplyUncertainties( std::vector< TLorentzVector >&,
			   std::vector< TLorentzVector >& );

 private:
  UncertaintyTool*  m_uncertaintyTool;
};


#endif
