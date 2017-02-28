#ifndef JETPAIR_H
#define JETPAIR_H

#include <TLorentzVector.h>

class JetPair{
 public:
  JetPair();
  JetPair( TLorentzVector*, TLorentzVector*, double );
  ~JetPair();

  TLorentzVector* RecoJet() { return m_recoJet; }
  TLorentzVector* TruthJet(){ return m_truthJet; }
  
 private:
  
  TLorentzVector* m_recoJet;
  TLorentzVector* m_truthJet;
  double          m_deltaR;
};

#endif
