#include "JetPair.h"

JetPair::JetPair()
  : m_recoJet( NULL ), m_truthJet( NULL ), m_deltaR( 0 ) {}


JetPair::JetPair( TLorentzVector* rJ, TLorentzVector* tJ, double dR)
  : m_recoJet( rJ ), m_truthJet( tJ ), m_deltaR( dR ) {}

JetPair::~JetPair(){}
