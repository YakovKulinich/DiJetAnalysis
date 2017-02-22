#ifndef __THMULF_H__
#define __THMULF_H__

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>

#include <TH1.h>
#include <TClonesArray.h>

#ifndef __THMULF_VARBINNING__
#define __THMULF_VARBINNING__
#endif

class TAxis;
class TTreeFormula;

/** A multi-dimension histogram.
    Definition of the axis are done through AddAxis() methods.
    Note the slightly unusual Fill method, where the first parameter
    is the weight (usually 1.0).
    If you use weights != 1.0, use Sumw2() prior to any filling, otherwise
    the errors won't be computed correctly.
    The projection to something manageable is mainly done through 
    Draw() method. Currently supported projection dimensions are 1,2 and 3.
    Note also that, thanks to the derivation from TH1, you can seamlessly
    use Fun4AllServer:registerHisto() method with THmulf objects (and thus use
    Fun4AllServer::dumpHistos() later on).
*/

class THmulf : public TH1
{
 public:
  THmulf();
  THmulf(const char* name, const char* title);
  virtual ~THmulf();

  /// Add an axis.
  bool AddAxis(const char* name, const char* title, 
	       int nbins, float low, float up);

  /// Add an axis.  This version allows variable bin widths
  bool AddAxis(const char* name, const char* title, 
	       int nbins, double* lowedges);

  /// Make a new axis, given another axis.
  bool AddAxis(TAxis* axis);

  /// Add all the axis contained in th into this.
  bool AddAxis(THmulf* th);

  /// Fill with a weigth.
  int Fill(float w,
	   float x0, float x1 = 0, float x2 = 0, float x3 = 0, float x4 = 0,
	   float x5 = 0, float x6 = 0, float x7 = 0, float x8 = 0, float x9 = 0,
	   float x10 = 0, float x11 = 0, float x12 = 0, float x13 = 0);

  /// Core filling method.
  int Fill(float* p, float weight = 1);

  /// Add c1*th to this.
  //bool Add(THmulf* th, double c1 = 1);

  // Add with throw a bad_cast if h1 can't be cast to a THmulf
  virtual bool Add(const TH1* h1, double c1 = 1.0); 

  /// Find the bin corresponding the tuple p[]
  int FindBin(float* p);

  /// Get the bin content.
  virtual Stat_t GetBinContent(Int_t bin) const;
  
  /** Get the bin error.
      The error is sqrt(sumofweights) if Sumw2 has been used,
      sqrt(bincontent) otherwise.
  */
  virtual Stat_t GetBinError(Int_t bin) const;

  // Project this THmulf into a subset THmulf based on the 
  // variables specified in varexp.
  THmulf* Project(const char* varexp, const char* selection = "",
		  const char* option = "");

  // Project varexp into a 1, 2, or 3 dimensional Histogram, with 
  // selection applied.  If the hist named name exists in the current
  // directory or pad, it is reset and refilled.
  TH1* Projection(const char* name, const char* varexp, const char* selection = "",
		  const char* option = "" /* unused */);

  /// Draw.
  //TH1* Draw(const char* varexp, const char* selection = "", char* option = "");
  //void Draw(const char* varexp, const char* selection = "", char* option = "");
  virtual TH1* Draw(const char* varexp, const char* selection, const char* option="");

  /// Sets the range for a given axis.
  bool SetRange(const char* axisname, int min, int max);

  /// Get a given axis.
  TAxis* GetAxis(const char* axisname);

  /// Our list of internal axis.
  TClonesArray* GetAxis()
  {
    return _axis_array;
  };

  /// Reset.
  virtual void Reset(Option_t* opt = "");

  /// Print.
  virtual void Print(const char* opt = "") const;

  /// Add 1.0 to bin i.
  void AddBinContent(Int_t i);

  /** Add w to bin i. If Sumw2 has been called, also update
      the bin error.
  */
  void AddBinContent(Int_t i, Stat_t w);

  /// Must be called if Fill with weights is to be used.
  void Sumw2();

  /// Our dimension (i.e. our number of axis).
  int GetDimension() const { return _naxis; }

  /// How many Fill were performed.
  Stat_t GetEntries() const { return _nentries; }

  // Return how many bins there are 
  int GetNContent() const { return _ncont; }

  // Return the raw contents of the THmulf
  void GetBinContents(float* a, int len);

  // Return the pointer to the cached histogram
  const TH1* GetHist() const
  {
    if (_hist_initialized) return _hist;
    else return 0;
  }

  void SetOptDebug(bool debug = true) { _debug=debug; }

  bool GetOptDebug() const { return _debug; }

  int SetVerbosity(int lvl = 0) {
    int old = _verbosity;
    _verbosity = lvl;
    return old;
  }

  // Make a topological comparison between two THmulfs
  // Returns 0 if they are the "same"
  static int CompareTopology(THmulf& a, THmulf& b, int verbose=0);

  // These functions are supplied by the base TH1, but are overridden in THmulf
  // with slightly different signatures.  Because the functions below are not
  // actually well-defined for THmulf, this isn't exactly an error.  We do, however,
  // need to supply the proper polymorphic versions of them to allow inheritance
  // to behave in an expected manner.  Since these functions are not technically useful
  // for a THmulf, they raise exceptions if ever used.
  virtual Int_t GetBin(int, int, int) const { bad_call("THmulf::GetBin(int,int,int)"); return 0; }
  virtual Int_t FindBin(double, double, double) { bad_call("THmulf::FindBin(double,double,double)"); return 0;}
  virtual Stat_t GetBinContent(int, int, int) const { bad_call("THmulf::GetBinContent(int,int,int)"); return 0;}
  virtual Stat_t GetBinContent(int, int) const { bad_call("THmulf::GetBinContent(int,int)"); return 0;}
  virtual Stat_t GetBinError(int, int, int) const { bad_call("THmulf::GetBinError(int,int,int)"); return 0;}
  virtual Stat_t GetBinError(int, int) const { bad_call("THmulf::GetBinError(int,int)"); return 0;}
  virtual Int_t Fill(const char*, double) { bad_call("THmulf::Fill(const char*,double)"); return 0; }
  virtual Int_t Fill(double, double) { bad_call("THmulf::Fill(double,double)"); return 0; }
  virtual Int_t Fill(double) { bad_call("THmulf::Fill(double)"); return 0; }
  virtual bool Add(const TH1*, const TH1*, double, double) { bad_call("THmulf::Add(const TH1*,const TH1*,double,double)"); return 0;}
  virtual bool Add(TF1*, double, const Option_t*) { bad_call("THmulf::Add(TF1*,double,const Option_t*)"); return 0;}

  // Draw() is a little different.  THmulf::Draw() projects as well as draws, and thus
  // has an entirely different signature as well as behavior.  The approach here is simply
  // to print an error message and make a gentle suggestion to the user.
  virtual void Draw(const Option_t* s = "") {
    std::cout << "Calling THmulf::Draw(\"" << s
  	      << "\") prohibited, try THmulf::Draw(\"" << s 
  	      << "\",\"\") instead."
  	      << std::endl;
  }

  /** Get our internal bin correspond to the array of TH1-like bins ibin[].
      (reverse of DecodeBin).
  */
  int GetBin(int* ibin);

  /// From our internal bin to TH1-like bins ibin[] (reverse of GetBin)
  void DecodeBin(int bin, int* ibin) const;

private:
  void bad_call(const char* s) const {
    char msg[1024];
    snprintf(msg,1024,"ERROR: Calling %s prohibited",s);
    std::cout << msg << std::endl;
    throw std::runtime_error(msg);
  }

  bool AnalyzeVar(const char* varexp);
  bool CreateHist(const char* varexp, const char* title = "", const char* optionalName=0);
  
  int FillHist1D(int iaxis);
  int FillHist2D(int iaxis);
  int FillHist3D(int iaxis);

  int WFillHist1D(int iaxis);
  int WFillHist2D(int iaxis);
  int WFillHist3D(int iaxis);
  
  bool CreateTHmulf(const char* varexp, const char* title = "");
  int FillTHmulfFastLoop(int iaxis);

  bool _hist_initialized;            //!
  TH1* _hist;                        //! Created histogram by Draw()
  bool _thmul_initialized;           //!
  THmulf* _thmul;                    //! Created THmulf by Project()
  bool _debug;                       //! Debug option
  int _verbosity;                    //! verbosity level

  // These are used for data storage.
  int _naxis;                        //  Number of DIMENSION(axis)
  TClonesArray* _axis_array;         //->Array of axis
  int _nentries;                     //  Number of entries
  int _ncont;                        //  Number of contents
  float* _cont;                      //[_ncont]
  bool _initialized;                 //  initialization of _cont[]
  bool Initialize();

  static const int MAXDIMENSION=25;

  // These are used only in Draw() and Project() method.
  int _fill_dim;                            //!
  int _fill_dim_iaxis[MAXDIMENSION];        //! _fill_dim size < _naxis(DIMENSION)
  TAxis* _fill_dim_axis[MAXDIMENSION];      //! _fill_dim size < _naxis(DIMENSION)
  float _fill_dim_cont[MAXDIMENSION];       //! _fill_dim size < _naxis(DIMENSION)
  TTreeFormula* _fill_formula;              //!
  int _fill_axis_ibin[MAXDIMENSION];        //! _naxis size
  float _fill_axis_cont[MAXDIMENSION];      //! _naxis size
  
  //  ClassDef(THmulf, 3) // Multi-DIMENSIONAL Histogram Class
    
};

#endif
