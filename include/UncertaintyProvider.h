#ifndef UNCERTAINTYPROVIDER_H
#define UNCERTAINTYPROVIDER_H

class TLorentzVector;

class UncertaintyTool{

 public:
  UncertaintyTool( int );
  virtual ~UncertaintyTool();

  virtual void ApplyUncertainty( TLorentzVector& jet, double ) = 0;

  int m_uc;
  int m_sign;
};

class JERUncertaintyTool : public UncertaintyTool{
 public:
  JERUncertaintyTool( int );
  ~JERUncertaintyTool();

  void ApplyUncertainty(  TLorentzVector& jet, double = 0 );
};


class JESUncertaintyTool : public UncertaintyTool{
 public:
  JESUncertaintyTool( int );
  ~JESUncertaintyTool();

  void ApplyUncertainty(  TLorentzVector& jet, double = 0 );
};

class UncertaintyProvider{

 public:
  UncertaintyProvider();
  UncertaintyProvider( int );
  virtual ~UncertaintyProvider();

 private:
  void ApplyUncertainty( TLorentzVector& jet, double = 0 );
  
 private:
  UncertaintyTool*  m_uncertaintyTool;
};


#endif
