#ifndef UNCERTAINTYTOOL_H
#define UNCERTAINTYTOOL_H

class TLorentzVector;

class JERUncertaintyTool{
 public:
  JERUncertaintyTool();
  ~JERUncertaintyTool();

  ApplyUncertainty();
};

class UncertaintyTool{

 public:
  UncertaintyTool();
  UncertaintyTool( int );
  virtual ~UncertaintyTool();

  void Initialize();

 private:
  void ApplyUncertainty( TLorentzVector& jet, double = 0 ) = 0;
  
 private:
  int m_uc;
  UncertaintyFunction m_fUncertaintyMethod;
};

#endif
