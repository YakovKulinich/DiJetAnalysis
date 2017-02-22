#include "THmulf.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>

#include <TTreeFormula.h>
#include <TTree.h>
#include <TH2.h>
#include <TH3.h>
#include <TStopwatch.h>
#include <TROOT.h>

// ClassImp(THmulf);

using namespace std;

//=============================================================================
THmulf::THmulf() :
  TH1(),
  _hist_initialized(false),
  _hist(0), 
  _thmul_initialized(false),
  _thmul(0), 
  _debug(false),
  _verbosity(0),
  _naxis(0),
  _axis_array(new TClonesArray("TAxis",MAXDIMENSION)),
  _nentries(0),
  _ncont(0), 
  _cont(0),
  _initialized(false),
  _fill_dim(0),
  _fill_formula(0)
{}

//============================================================================
THmulf::THmulf(const char* name, const char* title) :
  TH1(name, title, 0, 0, 0),
  _hist_initialized(false),
  _hist(0), 
  _thmul_initialized(false),
  _thmul(0), 
  _debug(false),
  _verbosity(0),
  _naxis(0),
  _axis_array(new TClonesArray("TAxis",MAXDIMENSION)),
  _nentries(0),
  _ncont(0),
  _cont(0), 
  _initialized(false),
  _fill_dim(0),
  _fill_formula(0)
{}

//============================================================================
THmulf::~THmulf()
{
  delete[] _cont;
  _axis_array->Delete();
  delete _axis_array;
}

//============================================================================
bool 
THmulf::Initialize()
{
  if ( _naxis == 0 )
    {
      return false;
    }

  delete[] _cont;

  _ncont = 1;
  int iaxis = _naxis;
  while ( iaxis-- )
    {
      TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
      _ncont *= a->GetNbins();
    }

  _cont = new float[_ncont];

  memset(_cont,0,_ncont*sizeof(float));
  
  _initialized = true;
  return true;
}

//============================================================================
bool 
THmulf::AddAxis(const char* name, const char* title, 
		int nbins, float low, float up)
{
  if ( up < low )
    {
      return false;
    }
  if ( _naxis + 1 >= MAXDIMENSION )
    {
      return false;
    }
  int iaxis = _axis_array->GetEntries();

  while ( iaxis-- )
    {
      TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
      if ( strcmp(name,a->GetName() ) == 0 )
        {
	  std::cerr << " " << GetName() 
		    << "::AddAxis() receive same axis name." 
		    << std::endl;
          return false;
        }
    }

  TClonesArray& array = *_axis_array;
  TAxis* axis = new(array[_naxis]) TAxis(nbins, low, up);
  axis->SetName(name);
  axis->SetTitle(title);

  _naxis++;

  return Initialize();
}

bool
THmulf::AddAxis(const char* name, const char* title, 
		int nbins, double* lowedges) 
{
  // Check that we won't overflow the allocated memory
  //
  if ( _naxis + 1 >= MAXDIMENSION ) return false;

  // Iterate over the axes and check that an axis with this name
  // hasn't already been created.
  //
  int iaxis = _axis_array->GetEntries();
  while ( iaxis-- )
    {
      //TAxis& a = ((TAxis*)_axis_array)[iaxis];
      TAxis *a = (TAxis*)_axis_array->operator[](iaxis);
      if ( strcmp(name,a->GetName() ) == 0 )
        {
	  std::cerr << " " << GetName() 
		    << "::AddAxis() received axis name already exists." 
		    << std::endl;
          return false;
        }
    }

  TClonesArray& array = *_axis_array;
  TAxis* axis = new(array[_naxis]) TAxis(nbins, lowedges);
  axis->SetName(name);
  axis->SetTitle(title);

  _naxis++;

  return Initialize();
}

//============================================================================
bool
THmulf::AddAxis(TAxis* axis )
{
  const Int_t nbins = axis->GetNbins();
  std::vector<double> edges(nbins);
  axis->GetLowEdge(&edges[0]);
  edges.push_back(axis->GetBinUpEdge(nbins)); // need upper edge of last bin for lower edge of overflow

  return AddAxis(axis->GetName(), axis->GetTitle(),
                 nbins, &edges[0]);
}

//============================================================================
bool 
THmulf::AddAxis(THmulf* th)
{
  TClonesArray* array = th->_axis_array;

  for ( int iaxis = 0 ; iaxis < array->GetEntries(); iaxis++ )
    {
      TAxis* a = (TAxis*)array->operator[](iaxis);
      if (! AddAxis(a) )
	{
	  return false;
	}
    }
  return true;
}

//============================================================================
void
THmulf::Sumw2()
{
  if (!_initialized)
    {
      Initialize();
    }

  if (fSumw2.fN) 
    {
      Warning("Sumw2","Sum of squares of weights structure already created");
      return;
    }
  
  fSumw2.Set(_ncont);

  for (Int_t bin=0; bin<_ncont; bin++) 
    {
      fSumw2.fArray[bin] = GetBinContent(bin);
    }
}

//============================================================================
int 
THmulf::Fill(float *p, float weight)
{
  if (!_initialized )
    {
      Initialize();
    }

  int bin = FindBin(p);

  if ( bin < 0 || bin >= _ncont )
    {
      // Unlike ROOT TH1 or TH2, we do not deal with
      // underflow or overflow.
      // Simply returns false if p[] values are not within
      // our bin limits.

      return false;
    }

  // Fill the content of the bin...
  assert(bin>=0&&bin<_ncont);
  _cont[bin] += weight;

  // ...and the sum of the weigths, if available
  if ( fSumw2.fN )
    {
      assert(bin>=0 && bin<fSumw2.fN);
      fSumw2.fArray[bin] += weight*weight;
    }

  ++_nentries;

  return true;
}

//============================================================================
int 
THmulf::Fill(float w, float x0, float x1, float x2, float x3, float x4,
	     float x5, float x6, float x7, float x8, float x9,
	     float x10, float x11, float x12, float x13)
{
  float f[14];
  f[0] = x0;
  f[1] = x1;
  f[2] = x2;
  f[3] = x3;
  f[4] = x4;
  f[5] = x5;
  f[6] = x6;
  f[7] = x7;
  f[8] = x8;
  f[9] = x9;
  f[10] = x10;
  f[11] = x11;
  f[12] = x12;
  f[13] = x13;
  return this->Fill(f, w);
}

//============================================================================
//void 
bool
THmulf::Add(const TH1* h, double c1)
{
  const THmulf* th = dynamic_cast<const THmulf*>(h);

  if ( ! th ) throw std::bad_cast();

  if ( _ncont != th->_ncont )
    {
      return false;
      //return;
    }
  int n = _ncont;
  while ( n-- )
    {
      _cont[n] += c1*th->_cont[n];
      if ( fSumw2.fN )
	{
	  fSumw2[n] += c1*c1*th->fSumw2[n];
	}
    }

  _nentries += th->_nentries;
  return true;
  //return;
}

//============================================================================
void 
THmulf::AddBinContent(Int_t bin)
{
  if ( bin >= 0 && bin < _ncont )
    {
      ++_cont[bin];
    }
}

//============================================================================
void 
THmulf::AddBinContent(Int_t bin, Stat_t w)
{
  if ( bin >= 0 && bin < _ncont )
    {
      _cont[bin] += w;
      if ( fSumw2.fN )
	{
	  fSumw2[bin] += w*w;
	}
    }
}

//============================================================================
Stat_t 
THmulf::GetBinContent(Int_t bin) const
{
  if ( bin >= 0 && bin < _ncont )
    {
      return _cont[bin];
    }
  else
    {
      return 0;
    }
}

//_____________________________________________________________________________
Stat_t
THmulf::GetBinError(Int_t bin) const
{
  if ( bin >= 0 && bin < _ncont )
    {
      if ( fSumw2.fN )
	{
	  assert(bin<fSumw2.fN);
	  return sqrt(fSumw2[bin]);
	}
      else
	{
	  return sqrt(GetBinContent(bin));
	}
    }
  else
    {
      std::cerr << __FILE__ << ":" << __LINE__ 
		<< " bin=" << bin << " _ncont=" <<_ncont
		<< std::endl;
      return 0;
    }
}

//============================================================================
int 
THmulf::GetBin(int *ibin)
{
  // ibin[] array contains the (bin numbers of TAxis) minus one !
  int iaxis = _naxis;
  int nbin;
  int bin = 0;
  while ( iaxis-- )
    {
      TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
      nbin = a->GetNbins();
      bin = bin * nbin + ibin[iaxis];
    }
  return bin;
}

//============================================================================
int 
THmulf::FindBin(float *p)
{
  int iaxis = _naxis;
  int bin = 0;

  while ( iaxis-- )
    {
      TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
      int nbin = a->GetNbins();
      int ibin = a->FindFixBin(p[iaxis]);
      if ( ibin > 0 && ibin <= nbin )
	{
	  // For ROOT, bin=0 is underflow and bin=nbin+1 is overflow.
	  bin = bin*nbin + (ibin-1);
	}
      else
	{
	  return -1;
	}
    }
  return bin;
}

//============================================================================
bool 
THmulf::AnalyzeVar(const char* varexp)
{
  //---------------------------------------------------------------------------
  // ---------- Finding histgram's axis...
  int iaxis;
  size_t pre_nvar = 0;
  int exp_dim = 0;
  _fill_dim = 0;
  char var[128];

  for ( size_t nvar = 0; nvar <= strlen(varexp) ; ++nvar)
    {
      if ( varexp[nvar] == ':' )
	{
	  ++exp_dim;
	}

      if ( _debug ) 
	{
	  cout << " " << GetName() 
	       << "::AnalyzeVar() .. searching varexp : " 
	       << varexp[nvar] << endl;
	}

      if ( varexp[nvar] == ':' || 
	   varexp[nvar] == '>' || 
	   nvar == strlen(varexp) )
        {
          if ( _debug )
	    {
	      cout << " " << GetName() << "::AnalyzeVar() nvar-pre_nvar = " 
		   << nvar - pre_nvar << endl;
	    }
          strncpy(var, varexp + pre_nvar, nvar - pre_nvar);
          var[nvar - pre_nvar] = '\0';
          if ( _debug )
	    {
	      cout << " " << GetName() << "::AnalyzeVar() var = " 
		   << var << endl;
	    }
          iaxis = _naxis;
          while ( iaxis-- )
            {
              TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
              if ( strcmp(var,a->GetName()) == 0 )
                {
                  _fill_dim_iaxis[_fill_dim] = iaxis;
                  _fill_dim_axis[_fill_dim] = a;
                  if ( _debug )
		    {
		      cout << " " << GetName() 
			   << "::AnalyzeVar()  Found matched axis : " 
			   << a->GetName() << endl;
		    }
                  _fill_dim++;
                }
            }
          pre_nvar = nvar + 1;
        }
    }
  if ( _debug )
    {
      cout << " " << GetName() << "::AnalyzeVar() dimension = " 
	   << _fill_dim << endl;
    }

  if (_fill_dim == 0 )
    {
      return 0;
    }

  if (_fill_dim != exp_dim + 1 )
    {
      cout << " " << GetName() 
	   << "::AnalyzeVar() Error!! varexp are not correct " << endl;
      return 0;
    }
  return true;
}

//============================================================================
bool 
THmulf::CreateHist(const char* varexp, const char* title,
		   const char* optionalName)
{
  if ( _fill_dim > 3 || _fill_dim <= 0 )
    {
      cout << " " << GetName() 
	   << "::CreateHist() can't create histogram with " 
	   << _fill_dim << " DIMENSION " << endl;
      return 0;
    }

  //--------------------------------------------------------------------------
  // ---------- Find indicated histgram name followed by ">>" in varexp
  // If optional name specified, use that instead
  char hname[1024], htitle[1024];
  bool hexist = false;
  bool hnamed = false;
  if ( optionalName && strlen(optionalName) )
    {
      sprintf(hname,optionalName);
      _hist = (TH1*)gDirectory->Get(hname);
      hnamed = true;
      hexist = ( _hist !=0 );
    }
  else
    {
      for (int k = strlen(varexp) - 1 ; k > 0 ; k--)
	{
	  if (varexp[k] == '>' && varexp[k - 1] == '>')
	    {
	      if ( strlen(varexp) - k - 1 )
		{
		  strncpy( hname, varexp + k + 1, strlen(varexp) - k - 1 );
		  hname[strlen(varexp) - k - 1] = '\0';
		  _hist = (TH1*)gDirectory->Get(hname);
		  hnamed = true;
		  hexist = ( _hist !=0 );
		}
	    }
	}
    }

  //--------------------------------------------------------------------------
  // ---------- Creating histogram
  if (! hexist )
    {
      if (! hnamed )
	{
	  sprintf(hname, "%s_hist", GetName());
	}
      sprintf(htitle, "%s %s", GetTitle(), title);

      std::vector<int> nbins;
      std::vector<std::vector<double> > edges;
      for (int i=0; i<_fill_dim; ++i)
        {
	  const int nbins_ = _fill_dim_axis[i]->GetNbins();
          std::vector<double> lowedges(nbins_+1);
	  _fill_dim_axis[i]->GetLowEdge(&lowedges[0]);
	  lowedges[nbins_] =_fill_dim_axis[i]->GetBinLowEdge(nbins_+1);
          nbins.push_back(nbins_);
          edges.push_back(lowedges);
        }

      if (_fill_dim==1) {
        _hist = new TH1F(hname, htitle, nbins[0], &(edges[0])[0]);
      } else if (_fill_dim==2) {
        _hist = new TH2F(hname, htitle, nbins[1], &(edges[1][0]), nbins[0], &(edges[0][0]));
      } else if (_fill_dim==3) {
        _hist = new TH3F(hname, htitle, nbins[2], &(edges[2][0]), nbins[1], &(edges[1][0]), nbins[0], &(edges[0][0]));
      }
      _hist_initialized = true;
    }

  assert(_hist!=0);

  _hist->Reset();

  if ( fSumw2.fN )
    {
      _hist->Sumw2();
    }

  return true;
}

//============================================================================
bool 
THmulf::CreateTHmulf(const char* varexp, const char* title )
{
  //  if( _thmul_initialized ) delete _thmul;
  //--------------------------------------------------------------------------
  // ---------- Creating THmulf
  char hname[1024], htitle[1024];
  int dim;
  TString st;
  for ( dim = 0 ; dim < _fill_dim; dim++ )
    {
      st += '_';
      st += _fill_dim_axis[dim]->GetName();
    }
  sprintf(hname, "%s%s", GetName(), st.Data());
  sprintf(htitle, "%s %s", GetTitle(), title);
  _thmul = new THmulf(hname, htitle);

  for ( dim = 0 ; dim < _fill_dim; dim++ )
    {
      _thmul->AddAxis(_fill_dim_axis[dim]->GetName(), 
		      _fill_dim_axis[dim]->GetTitle(),
                      _fill_dim_axis[dim]->GetNbins(), 
		      _fill_dim_axis[dim]->GetXmin(), 
		      _fill_dim_axis[dim]->GetXmax());
    }
  _thmul_initialized = true;
  return true;
}

//============================================================================
THmulf* 
THmulf::Project(const char* varexp, const char* selection, const char* opt)
{
  if (! AnalyzeVar(varexp) )
    {
      return 0;
    }

  if (! CreateTHmulf(varexp, selection) )
    {
      return 0;
    }

  //--------------------------------------------------------------------------
  // ---------- Call TTreeFormula
  _fill_formula=0;
  TTree* tree=0;

  if ( strlen(selection) )
    {
      char hname[1024];
      int iaxis = _naxis;
      tree = new TTree();
      while ( iaxis-- )
	{
	  TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
	  sprintf(hname,"b_%s", a->GetName());
	  tree->Branch(hname, &_fill_axis_cont[iaxis], a->GetName());
	}
      _fill_formula = new TTreeFormula("tf", selection, tree);
    }

  //--------------------------------------------------------------------------
  // ---------- Fill THmulf

  TStopwatch timer;
  timer.Start();
  FillTHmulfFastLoop(_naxis - 1);
  cout.precision(3);
  if ( _verbosity > 1 )
    {
      cout << " " << GetName() << "::Project() FillTHmulfFastLoop(" 
	   << _naxis - 1 << ") took RealTime="
	   << timer.RealTime() << " CpuTime=" << timer.CpuTime() << "sec " << endl;
    }
  delete tree;
  return _thmul;
}

//============================================================================
TH1* 
THmulf::Draw(const char* varexp, const char* selection, const char* opt)
{
  //-------------------------------------------------------------------------
  // ---------- Call Initialization
  if (! AnalyzeVar(varexp) )
    {
      return 0;
      //return;
    }
  if (! CreateHist(varexp, selection) )
    {
      return 0;
      //eturn;
    }

  //---------------------------------------------------------------------------
  // ---------- Call TTreeFormula
  char hname[1024];
  int iaxis = _naxis;
  TTree* tree = 0;
  _fill_formula = 0;

  if ( strlen(selection) )
    {
      tree = new TTree();
      while ( iaxis-- )
	{
	  TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
	  sprintf(hname,"b_%s", a->GetName());
	  tree->Branch(hname, &_fill_axis_cont[iaxis], a->GetName());
	}
      _fill_formula = new TTreeFormula("tf", selection, tree);
    }

  //--------------------------------------------------------------------------
  // ---------- Call FillHistLoop....
  int status=0;

  TStopwatch timer;
  timer.Start();
  cout.precision(3);

  // Try to put the more ifs here instead of within the loop Fill*Hist*
  if ( _fill_dim == 1 )
    {     
      if (fSumw2.fN)
	{
	  status = WFillHist1D(_naxis-1);
	}
      else
	{
	  status = FillHist1D(_naxis-1);
	}
    }
  else if ( _fill_dim == 2 )
    {
      if (fSumw2.fN)
	{
	  status = WFillHist2D(_naxis-1);
	}
      else
	{
	  status = FillHist2D(_naxis-1);
	}
    }
  else if ( _fill_dim == 3 )
    {
      if ( fSumw2.fN )
	{
	  status = WFillHist3D(_naxis-1);
	}
      else
	{
	  status = FillHist3D(_naxis-1);
	}
    }
  else
    {
      assert(0==1);
    }

  if ( _verbosity > 1 )
    {
      cout << "took RealTime=" << timer.RealTime() << " CpuTime=" 
	   << timer.CpuTime() << "sec " << endl;
    }

  _hist->Draw(opt);

  delete _fill_formula;
  delete tree;

  if ( status > 0 )
    {
      return _hist;
    }
  else
    {
      return NULL;
    }
}

//============================================================================
int 
THmulf::FillTHmulfFastLoop(int iaxis)
{
  int bin;
  int status = 0;
  if ( iaxis >= 0 )
    {
      //------ Check this axis is for filling ?
      int idim = 0;
      bool b_dim = false;
      int d = _fill_dim;
      while ( d-- )
        {
          if ( _fill_dim_iaxis[d] == iaxis )
            {
              idim = d;
              b_dim = true;
            }
        }

      //------ looping
      _fill_axis_ibin[iaxis] = 
	((TAxis*) _axis_array->operator[](iaxis))->GetLast();

      bin = ((TAxis*) _axis_array->operator[](iaxis))->GetFirst();

      while ( _fill_axis_ibin[iaxis]-- >= bin )
        {
          if ( b_dim )
	    {
	      _fill_dim_cont[idim] = ((TAxis*) _axis_array->operator[](iaxis))->GetBinCenter(_fill_axis_ibin[iaxis] + 1);
	    }
          _fill_axis_cont[iaxis] = ((TAxis*) _axis_array->operator[](iaxis))->GetBinCenter(_fill_axis_ibin[iaxis] + 1);
          //----------------------------------------
          status += FillTHmulfFastLoop(iaxis - 1);
          //----------------------------------------
        }
    }
  else
    {
      //------ Event selection
      bool fill = true;
      if ( _fill_formula )
        {
          if (! _fill_formula->EvalInstance() )
	    {
	      fill = false;
	    }
        }
      //------ Filling
      if ( fill )
        {
          bin = GetBin(_fill_axis_ibin);
          if ( _cont[bin] != 0 )
	    {
	      _thmul->Fill(_fill_dim_cont, _cont[bin]);
	    }
          status = 1;
        }
    }
  return status;
}

//_____________________________________________________________________________
int 
THmulf::WFillHist1D(int iaxis)
{
  int bin;
  int status = 0;
  
  if ( iaxis >= 0 )
    {
      //------ looping
      TAxis* axis = static_cast<TAxis*>(_axis_array->operator[](iaxis));
      assert(axis!=0);
      _fill_axis_ibin[iaxis] = axis->GetLast();
      bin = axis->GetFirst();
      while ( _fill_axis_ibin[iaxis]-- >= bin )
        {
	  // Warning: once here, _fill_axis_ibin[iaxis] *has been*
	  // decremented by 1 already !

	  _fill_axis_cont[iaxis] = 
	    axis->GetBinCenter(_fill_axis_ibin[iaxis]+1);

          status += WFillHist1D(iaxis - 1);
        }
    }
  else
    {
      bool fill = true;

      if ( _fill_formula )
	{
	  if ( ! _fill_formula->EvalInstance() )
	    {
	      fill = false;
	    }
	}

      if (!fill) 
	{
	  return 0;
	}

      //------ Event selection
      int bin = GetBin(_fill_axis_ibin);
      assert(bin >= 0);
      int b = _fill_axis_ibin[_fill_dim_iaxis[0]]+1;
      int ibin = _hist->GetBin(b);	

      if ( _cont[bin] != 0 )
        {
	  double e1 = _hist->GetBinError(ibin);
	  _hist->Fill(_fill_axis_cont[_fill_dim_iaxis[0]],_cont[bin]);
	  double e = GetBinError(bin);
	  double newerror = sqrt(e*e+e1*e1);
	  _hist->SetBinError(ibin,newerror);
	}
      status=1;
    }
  return status;
}

//_____________________________________________________________________________
int 
THmulf::FillHist1D(int iaxis)
{
  int bin;
  int status = 0;
  
  if ( iaxis >= 0 )
    {
      //------ looping
      TAxis* axis = static_cast<TAxis*>(_axis_array->operator[](iaxis));
      assert(axis!=0);
      _fill_axis_ibin[iaxis] = axis->GetLast();
      bin = axis->GetFirst();
      while ( _fill_axis_ibin[iaxis]-- >= bin )
        {
	  // Warning: once here, _fill_axis_ibin[iaxis] *has been*
	  // decremented by 1 already !

	  _fill_axis_cont[iaxis] = 
	    axis->GetBinCenter(_fill_axis_ibin[iaxis]+1);

          status += FillHist1D(iaxis - 1);
        }
    }
  else
    {
      bool fill = true;

      if ( _fill_formula )
	{
	  if ( ! _fill_formula->EvalInstance() )
	    {
	      fill = false;
	    }
	}

      if (!fill) 
	{
	  return 0;
	}

      //------ Event selection
      int bin = GetBin(_fill_axis_ibin);
      assert(bin >= 0);

      if ( _cont[bin] != 0 )
        {
	  _hist->Fill(_fill_axis_cont[_fill_dim_iaxis[0]],_cont[bin]);
	}
      status=1;
    }
  return status;
}

//_____________________________________________________________________________
int 
THmulf::WFillHist2D(int iaxis)
{
  // You may note that this is one is very similar to FillHist1D.
  // Why 2 methods then ?
  // Well, for speed, we want to avoid as many if as possible in
  // this low level loop. That's it.

  int bin;
  int status = 0;
  
  if ( iaxis >= 0 )
    {
      //------ looping
      TAxis* axis = static_cast<TAxis*>(_axis_array->operator[](iaxis));
      assert(axis!=0);
      _fill_axis_ibin[iaxis] = axis->GetLast();
      bin = axis->GetFirst();
      while ( _fill_axis_ibin[iaxis]-- >= bin )
        {
	  // Warning: once here, _fill_axis_ibin[iaxis] *has been*
	  // decremented by 1 already !

	  _fill_axis_cont[iaxis] = 
	    axis->GetBinCenter(_fill_axis_ibin[iaxis]+1);

          status += FillHist2D(iaxis - 1);
        }
    }
  else
    {
      bool fill = true;

      if ( _fill_formula )
	{
	  if ( ! _fill_formula->EvalInstance() )
	    {
	      fill = false;
	    }
	}

      if (!fill) 
	{
	  return 0;
	}

      //------ Event selection
      int bin = GetBin(_fill_axis_ibin);
      assert(bin >= 0);
      int b1 = _fill_axis_ibin[_fill_dim_iaxis[1]]+1;
      int b2 = _fill_axis_ibin[_fill_dim_iaxis[0]]+1;

      if ( _cont[bin] != 0 )
        {
	  TH2* h2 = static_cast<TH2*>(_hist);
	  double e1 = h2->GetCellError(b1,b2);
	  h2->Fill(_fill_axis_cont[_fill_dim_iaxis[1]],
		   _fill_axis_cont[_fill_dim_iaxis[0]],
		   _cont[bin]);

	  double e = GetBinError(bin);
	  double newerror = sqrt(e*e+e1*e1);
	  h2->SetCellError(b1,b2,newerror);
	}
      status=1;
    }
  return status;
}

//_____________________________________________________________________________
int 
THmulf::FillHist2D(int iaxis)
{
  // You may note that this is one is very similar to FillHist1D.
  // Why 2 methods then ?
  // Well, for speed, we want to avoid as many if as possible in
  // this low level loop. That's it.

  int bin;
  int status = 0;
  
  if ( iaxis >= 0 )
    {
      //------ looping
      TAxis* axis = static_cast<TAxis*>(_axis_array->operator[](iaxis));
      assert(axis!=0);
      _fill_axis_ibin[iaxis] = axis->GetLast();
      bin = axis->GetFirst();
      while ( _fill_axis_ibin[iaxis]-- >= bin )
        {
	  // Warning: once here, _fill_axis_ibin[iaxis] *has been*
	  // decremented by 1 already !

	  _fill_axis_cont[iaxis] = 
	    axis->GetBinCenter(_fill_axis_ibin[iaxis]+1);

          status += FillHist2D(iaxis - 1);
        }
    }
  else
    {
      bool fill = true;

      if ( _fill_formula )
	{
	  if ( ! _fill_formula->EvalInstance() )
	    {
	      fill = false;
	    }
	}

      if (!fill) 
	{
	  return 0;
	}

      //------ Event selection
      int bin = GetBin(_fill_axis_ibin);
      assert(bin >= 0);
      if ( _cont[bin] != 0 )
        {
	  TH2* h2 = static_cast<TH2*>(_hist);
	  h2->Fill(_fill_axis_cont[_fill_dim_iaxis[1]],
		   _fill_axis_cont[_fill_dim_iaxis[0]],
		   _cont[bin]);
	}
      status=1;
    }
  return status;
}

//_____________________________________________________________________________
int 
THmulf::WFillHist3D(int iaxis)
{
  // You may note that this is one is very similar to FillHist1D.
  // Why 2 methods then ?
  // Well, for speed, we want to avoid as many if as possible in
  // this low level loop. That's it.

  int bin;
  int status = 0;
  
  if ( iaxis >= 0 )
    {
      //------ looping
      TAxis* axis = static_cast<TAxis*>(_axis_array->operator[](iaxis));
      assert(axis!=0);
      _fill_axis_ibin[iaxis] = axis->GetLast();
      bin = axis->GetFirst();
      while ( _fill_axis_ibin[iaxis]-- >= bin )
        {
	  // Warning: once here, _fill_axis_ibin[iaxis] *has been*
	  // decremented by 1 already !

	  _fill_axis_cont[iaxis] = 
	    axis->GetBinCenter(_fill_axis_ibin[iaxis]+1);

          status += WFillHist3D(iaxis - 1);
        }
    }
  else
    {
      bool fill = true;

      if ( _fill_formula )
	{
	  if ( ! _fill_formula->EvalInstance() )
	    {
	      fill = false;
	    }
	}

      if (!fill) 
	{
	  return 0;
	}

      //------ Event selection
      int bin = GetBin(_fill_axis_ibin);
      assert(bin >= 0);
      int b1 = _fill_axis_ibin[_fill_dim_iaxis[2]]+1;
      int b2 = _fill_axis_ibin[_fill_dim_iaxis[1]]+1;
      int b3 = _fill_axis_ibin[_fill_dim_iaxis[0]]+1;
      if ( _cont[bin] != 0 )
        {
	  TH3* h3 = static_cast<TH3*>(_hist);
	  int hbin = h3->GetBin(b1,b2,b3);
	  double e1 = h3->GetBinError(hbin);
	  h3->Fill(_fill_axis_cont[_fill_dim_iaxis[2]],
		   _fill_axis_cont[_fill_dim_iaxis[1]],
		   _fill_axis_cont[_fill_dim_iaxis[0]],
		   _cont[bin]);

	  double e = GetBinError(bin);
	  double newerror = sqrt(e*e+e1*e1);
	  h3->SetBinError(hbin,newerror);
	}
      status=1;
    }
  return status;
}

//_____________________________________________________________________________
int 
THmulf::FillHist3D(int iaxis)
{
  // You may note that this is one is very similar to FillHist1D.
  // Why 2 methods then ?
  // Well, for speed, we want to avoid as many if as possible in
  // this low level loop. That's it.

  int bin;
  int status = 0;
  
  if ( iaxis >= 0 )
    {
      //------ looping
      TAxis* axis = static_cast<TAxis*>(_axis_array->operator[](iaxis));
      assert(axis!=0);
      _fill_axis_ibin[iaxis] = axis->GetLast();
      bin = axis->GetFirst();
      while ( _fill_axis_ibin[iaxis]-- >= bin )
        {
	  // Warning: once here, _fill_axis_ibin[iaxis] *has been*
	  // decremented by 1 already !

	  _fill_axis_cont[iaxis] = 
	    axis->GetBinCenter(_fill_axis_ibin[iaxis]+1);

          status += FillHist3D(iaxis - 1);
        }
    }
  else
    {
      bool fill = true;

      if ( _fill_formula )
	{
	  if ( ! _fill_formula->EvalInstance() )
	    {
	      fill = false;
	    }
	}

      if (!fill) 
	{
	  return 0;
	}

      //------ Event selection
      int bin = GetBin(_fill_axis_ibin);
      assert(bin >= 0);
    
      if ( _cont[bin] != 0 )
        {
	  TH3* h3 = static_cast<TH3*>(_hist);	 
	  h3->Fill(_fill_axis_cont[_fill_dim_iaxis[2]],
		   _fill_axis_cont[_fill_dim_iaxis[1]],
		   _fill_axis_cont[_fill_dim_iaxis[0]],
		   _cont[bin]);
	}
      status=1;
    }
  return status;
}



//============================================================================
TAxis* 
THmulf::GetAxis(const char* axisname)
{
  int iaxis = _naxis;
  int naxis = -1;
  TAxis* a;
  while ( iaxis-- )
    {
      a = (TAxis*)_axis_array->operator[](iaxis);
      if ( strcmp(a->GetName(), axisname) == 0 )
	{
	  naxis = iaxis;
	}
    }
  if ( naxis != -1 )
    {
      return (TAxis*)_axis_array->operator[](naxis);
    }
  else
    {
      return NULL;
    }
}

//============================================================================
bool 
THmulf::SetRange(const char* axisname, int min, int max)
{
  TAxis* a = GetAxis(axisname);
  if ( a == NULL)
    {
      return false;
    }
  a->SetRange(min, max);
  return true;
}

//============================================================================
void 
THmulf::Reset(Option_t* opt)
{
  if ( _initialized )
    {
      memset(_cont,0,_ncont*sizeof(float));
      if ( fSumw2.fN )
	{
	  fSumw2.Reset();
	}
    }
}

//_____________________________________________________________________________
void
THmulf::DecodeBin(int bin, int* ibin) const
{
  // ibin must have a size >= GetDimension()
  // ibin will contain the bin numbers suitable for TAxis
  // (that is our bin numbers + 1) 

  int nbin = _ncont;

  int iaxis = _naxis;

  while ( iaxis-- )
    {
      TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
      nbin /= a->GetNbins();
      int x = bin/nbin;
      ibin[iaxis] = x+1;
      bin -= x*nbin;
 
    }  
}

//============================================================================
void 
THmulf::Print(const char* opt) const
{
  cout << " --------------------------------------------------- " << endl;
  cout << " THmulf class. name:" << GetName() 
       << " title:" << GetTitle() << endl;
  cout << " " << _axis_array->GetEntries() << " DIMENSIONS " << endl;
  if ( fSumw2.fN )
    {
      cout << " Has Sumw2" << endl;
    }
  int n;
  for (n = 0 ; n < _axis_array->GetEntries(); n++)
    {
      TAxis* a = (TAxis*) _axis_array->At(n);
      cout << " -- " << n << " -- " << a->GetName() 
	   << "\t\t" << a->GetTitle() << "\t"
	   << a->GetNbins() << "[" << a->GetXmin() 
	   << "," << a->GetXmax() << "] " << endl;
    }
  cout << " Number of entries : " << _nentries << endl;
  cout << " The total binning size : " << _ncont << endl;
  cout << " --------------------------------------------------- " << endl;
  
  std::string sopt = opt;
  if ( sopt=="CONTENT" || sopt=="content" )
    {
      int* ibin = new int[_naxis];

      for ( int i = 0; i < _ncont; ++i )
	{
	  if ( _cont[i] )
	    {
	      DecodeBin(i,ibin);
	      // I think this line was left over from debugging.  esp since
	      // it ignores the value of _naxis!!  Commenting out -- DLW.
	      //if ( ibin[0] == 28 && ibin[1] == 8 && ibin[2] == 1 )
		{
		  printf("_cont[%7d]=%e err=%e ",i,_cont[i],GetBinError(i));
		  for ( int j = 0; j < _naxis; ++j )
		    {
		      TAxis* a = static_cast<TAxis*>(_axis_array->At(j));
		      printf("%s[%3d]=%e ",a->GetName(),
			     ibin[j],
			     a->GetBinCenter(ibin[j]));
		    }
		  printf("\n");
		}
	    }
	}
      delete[] ibin;
    }
}

void
THmulf::GetBinContents(float* a, int len)
{
  int size = std::min(_ncont,len);
  std::copy(_cont,_cont+size,a);
}

static bool
compareAxes(TAxis* a, TAxis* b)
{
  if ( a->GetNbins() != b->GetNbins() ) return false;
  if ( a->GetXmax() != b->GetXmax() ) return false;
  if ( a->GetXmin() != b->GetXmin() ) return false;
  // Maybe check the edges of the bins instead?
  int nbins = a->GetNbins();
  Axis_t* a_edges = new Axis_t[nbins];
  Axis_t* b_edges = new Axis_t[nbins];
  boost::scoped_array<Axis_t> _a(a_edges);
  boost::scoped_array<Axis_t> _b(b_edges);
  a->GetLowEdge(a_edges);
  b->GetLowEdge(b_edges);
  if ( ! std::equal(a_edges,a_edges+nbins,b_edges) ) return false;
  return true;
}

static bool
compareContents(THmulf& a, THmulf& b, int& pos, float vals[2])
{
  int nbins1 = a.GetNContent();
  int nbins2 = b.GetNContent();
  if ( nbins1 != nbins2 ) return false;
  float* a_cont = new float[nbins1];
  float* b_cont = new float[nbins2];
  boost::scoped_array<float> _a(a_cont);
  boost::scoped_array<float> _b(b_cont);
  a.GetBinContents(a_cont,nbins1);
  b.GetBinContents(b_cont,nbins2);
  std::pair<float*,float*> match = std::mismatch(a_cont,a_cont+nbins1,b_cont);
  if ( match.first != a_cont+nbins1 ) {
    pos = int(match.first-a_cont);
    vals[0] = *match.first;
    vals[1] = *match.second;
    return false;
  }
  return true;
}

// Perform a topological comparison of two THmulfs
// This means the following checks are done:
// - identical dimensions
// - for each axis, same Nbins, Xmax, Xmin, and low edges
// - each element of _cont array has identical counts/values
//
// If 2 THmulfs pass all the above, they are considered the 
// same.  HOWEVER, it is certainly possible to construct a THmulf
// with the same geometry but different physics interpretation, albeit 
// this would be very rare.
// That's why I call it a topological comparison... it's worth what it's
// worth.  I think it will be most useful in checking output hists
// before/after subtle changes in things like cuts.
//
int
THmulf::CompareTopology(THmulf& a, THmulf& b, int verbose)
{
  // Compare number of dimensions
  if ( a.GetDimension() != b.GetDimension() ) {
    if ( verbose ) std::cout << "Dimension mismatch, " << a.GetDimension() << " != " << b.GetDimension()
			     << std::endl;
    return 1;
  }

  // Compare axes, one dimension at a time. Check the bin edges' values
  int ndim = a.GetDimension();
  TClonesArray* a_axes = a.GetAxis();
  TClonesArray* b_axes = b.GetAxis();
  for (int i=0; i!=ndim; ++i)
    {
      TAxis* a_axis = (TAxis*) a_axes->At(i);
      TAxis* b_axis = (TAxis*) b_axes->At(i);
      if ( ! compareAxes(a_axis,b_axis) ) {
	if ( verbose ) std::cout << "Axis mismatch, " << a_axis->GetName() << " != " << b_axis->GetName()
				 << std::endl;
	return 1;
      }
    }
  int pos = -1;
  float vals[2] = { 0.0 };
  if ( ! compareContents(a,b,pos,vals) ) {
    if ( verbose ) std::cout << "Contents mismatch, " << a.GetName() << " != " << b.GetName()
			     << ", Failed bin = " << pos
			     << ", values are = " << vals[0] << ", " << vals[1]
			     << std::endl;
    return 1;
  }

  return 0;
}


TH1*
THmulf::Projection(const char* name, const char* varexp, const char* selection, const char*)
{
  // Initialization
  if ( ! AnalyzeVar(varexp) ) return 0;

  // Create (or Reset) the output hist
  if ( ! CreateHist(varexp, selection, name) ) return 0;

  // Call TTreeFormula
  char hname[1024];
  int iaxis = _naxis;

  //delete _fill_formula;
  //_fill_formula = 0;

  TTree* tree = 0;

  if ( strlen(selection) )
    {
      tree = new TTree();
      while ( iaxis-- )
	{
	  TAxis* a = (TAxis*)_axis_array->operator[](iaxis);
	  sprintf(hname,"b_%s", a->GetName());
	  tree->Branch(hname, &_fill_axis_cont[iaxis], a->GetName());
	}
      _fill_formula = new TTreeFormula("tf", selection, tree);
    }

  boost::scoped_ptr<TTreeFormula> tmpform(_fill_formula); // facilitates clean up at end
  boost::scoped_ptr<TTree> tmptree(tree);  // facilitates clean up at end

  // ---------- Call FillHistLoop....
  int status = 0;

  TStopwatch timer;
  timer.Start();

  // Try to put the more ifs here instead of within the loop Fill*Hist*
  if ( _fill_dim == 1 )
    {     
      if (fSumw2.fN) status = WFillHist1D(_naxis-1);
      else status = FillHist1D(_naxis-1);
    }
  else if ( _fill_dim == 2 )
    {
      if (fSumw2.fN) status = WFillHist2D(_naxis-1);
      else status = FillHist2D(_naxis-1);
    }
  else if ( _fill_dim == 3 )
    {
      if ( fSumw2.fN ) status = WFillHist3D(_naxis-1);
      else status = FillHist3D(_naxis-1);
    }
  else
    {
      std::cout << "THmulf::Projection: Invalid number of dimensions (" << _fill_dim << ")" 
		<< std::endl;
      return 0;
    }

  std::cout.precision(3);
  if ( _verbosity > 1 )
    {
      cout << "THmulf::Projection: Projection took RealTime=" << timer.RealTime() << " CpuTime=" 
	   << timer.CpuTime() << "sec " << endl;
    }

  if ( status <= 0 ) return 0;
  
  return _hist;
}
