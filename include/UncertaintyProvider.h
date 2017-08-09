#ifndef UNCERTAINTYPROVIDER_H
#define UNCERTAINTYPROVIDER_H

#include <vector>

class TH1;
class TH2;
class TRandom3;
class TLorentzVector;

class UncertaintyTool{

 public:
  UncertaintyTool( int );
  virtual ~UncertaintyTool();

  virtual void ApplyUncertainty( TLorentzVector&, double ) = 0;
  virtual void ApplyUncertainty( TLorentzVector&,
				 TLorentzVector&, double ) {}

  int GetEtaUJERBin( float eta );

 protected:
  int m_uc;
  int m_sign;
};

class HIJERUncertaintyTool : public UncertaintyTool{
 public:
  HIJERUncertaintyTool( int );
  ~HIJERUncertaintyTool();

  void ApplyUncertainty(  TLorentzVector&, double = 0 );

 private:
  std::vector<TH1*> m_vJEShistos;
};

class JERUncertaintyTool : public UncertaintyTool{
 public:
  JERUncertaintyTool( int );
  ~JERUncertaintyTool();

  void ApplyUncertainty(  TLorentzVector&, double = 0 );
  virtual void ApplyUncertainty( TLorentzVector&,
				 TLorentzVector&, double = 0 );

 private:
  TRandom3* rand;
  
  TH2* hJER;
  std::vector<TH1*> m_vJERhistos;
};

class JESUncertaintyTool : public UncertaintyTool{
 public:
  JESUncertaintyTool( int );
  ~JESUncertaintyTool();

  void ApplyUncertainty(  TLorentzVector&, double = 0 );
};

class UncertaintyProvider{

 public:
  UncertaintyProvider();
  UncertaintyProvider( int );
  virtual ~UncertaintyProvider();

  void ApplyUncertainty( TLorentzVector&, double = 0 );
  
 private:
  UncertaintyTool*  m_uncertaintyTool;
};


#endif
