#ifndef DELTAPHIPROJ_H
#define DELTAPHIPROJ_H

#include <vector>
#include <string>

class TAxis;

enum PA { YS1, YS2, PT1, PT2 };

class DeltaPhiProj{

 public:
  DeltaPhiProj( int, int, int, int );
  ~DeltaPhiProj();

  int          GetAxisI    ( int i ){ return m_vAxisI[i]; } 
  std::string& GetAxisName ( int i ){ return m_vAxisName [ GetAxisI(i) ]; }
  std::string& GetAxisLabel( int i ){ return m_vAxisLabel[ GetAxisI(i) ]; }  
  TAxis*       GetTAxis    ( int i ){ return m_vTAxis    [ GetAxisI(i) ]; }  

  std::string& GetDefaultAxisName ( int i ){ return m_vAxisName [ i ]; }
  std::string& GetDefaultAxisLabel( int i ){ return m_vAxisLabel[ i ]; }  
  
  void AddTAxis( TAxis* ta ){ m_vTAxis.push_back( ta ); }

  std::vector<int> GetMappedBins( const std::vector<int>& );
  
  bool CorrectPhaseSpace( const std::vector<int>& );

 private:
  std::vector< int > m_vAxisI
    { 0, 1, 2, 3 }; // default values
  std::vector< std::string > m_vAxisName
    { "Ystar1" , "Ystar2", "Pt1", "Pt2", };
  std::vector< std::string > m_vAxisLabel
    { "#it{y}_{1}*", "#it{y}_{2}*", "#it{p}_{T}^{1}", "#it{p}_{T}^{2}" };
  std::vector< TAxis* > m_vTAxis;
};

#endif 
