// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <iostream>
#include <fstream>
#include <cstddef>
#include <math.h>
#include <complex>
#include <boost/math/special_functions/zeta.hpp>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
// -------- Classes ----------------------------------- //
#include "process.h"	// Process//
// -------- Functions --------------------------------- //
std::complex<double> N(double);
std::complex<double> Psi(std::complex<double>);
double Xcos(double);
double Xsin(double);
std::complex<double> MellinPDFDiff(double*,int,double);
std::complex<double> GlobalDiff(double*,double,double);
std::complex<double> TraceHSDiff(double*,int,double,complex<double>*,int);
std::complex<double> TraceBornDiff(double*,int,double);
std::complex<double> TraceHSExpDiff(double *, int , double , complex<double> *,int);
std::complex<double> ColinearDiff(double *, int , double );
std::complex<double> ColinearExpDiff(double *, int , double);

struct my_f_params { double sc; int mel; double M2; };  

extern "C"{
  //void update_mur_(double&);
  void ml5_0_sloopmatrix_thres_(double[4][4],complex<double>(*),double&,double[2][2],int&);
}

// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi, xmax, xmin;
extern int chan;

//------------------------------------------------------------//
//--------------------Total Integrand-------------------------//
//------------------------------------------------------------//
double TotVegas(double *x, size_t dim, void *params)
{
  double res=0., bet, t, hgg,  u, Xbet;
  std::complex<double> temp=0, i(0.0,1.0), Nb;
  
  struct my_f_params * fparams = (struct my_f_params *)params;  
  double M2=fparams->M2, sc=fparams->sc;
  int mel=fparams->mel;
  
  Nb=N(x[0])*exp(-Psi(1.));
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;
  //t=-M2/2.*(1-Xbet);

  if(mel==0)
    {
      std::complex<double> tmp=0.;      
      tmp=TraceBornDiff(x,chan,M2)*2.;
      res=std::imag(tmp*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));
    }

  else if(mel==1) //Diff
    {
      int Len=chan+2;
      int tu=0;
      complex<double> *xx;
      xx=new complex<double>[4*2*2];
      complex<double> *xx2;
      xx2=new complex<double>[4*2*2];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[2][2]={(0.,0.)};

      ppart[0][0]=pow(M2,0.5)/2.;
      ppart[1][0]=pow(M2,0.5)/2.;
      ppart[0][3]=pow(M2,0.5)/2.;
      ppart[1][3]=-pow(M2,0.5)/2.;

      ppart[2][0]=pow(M2,0.5)/2.;
      ppart[3][0]=pow(M2,0.5)/2.;

      ppart[2][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][3]=pow(M2,0.5)/2.*Xbet;
      ppart[3][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);

      res=std::imag(TraceHSDiff(x,chan,M2,xx,tu)*ColinearDiff(x,chan,M2)*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));

      res-=std::imag((TraceHSExpDiff(x,chan,M2,xx,tu)+(1.+ColinearExpDiff(x,chan,M2))*TraceBornDiff(x,chan,M2))*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));
      
      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);

      tu=1;
      res+=std::imag(TraceHSDiff(x,chan,M2,xx2,tu)*ColinearDiff(x,chan,M2)*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));

      res-=std::imag((TraceHSExpDiff(x,chan,M2,xx2,tu)+(1.+ColinearExpDiff(x,chan,M2))*TraceBornDiff(x,chan,M2))*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));
    }

    else if(mel==2) //NLL 
    {
      int tu=0;
      complex<double> *xx;
      xx=new complex<double>[4*2*2];
      complex<double> *xx2;
      xx2=new complex<double>[4*2*2];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[2][2]={(0.,0.)};

      ppart[0][0]=pow(M2,0.5)/2.;
      ppart[1][0]=pow(M2,0.5)/2.;
      ppart[0][3]=pow(M2,0.5)/2.;
      ppart[1][3]=-pow(M2,0.5)/2.;

      ppart[2][0]=pow(M2,0.5)/2.;
      ppart[3][0]=pow(M2,0.5)/2.;

      ppart[2][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][3]=pow(M2,0.5)/2.*Xbet;
      ppart[3][3]=-pow(M2,0.5)/2.*Xbet;
      
      //ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);

      res=std::imag(TraceHSDiff(x,chan,M2,xx,tu)*ColinearDiff(x,chan,M2)*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));
      
      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      //ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);
      tu=1;
      
      res+=std::imag(TraceHSDiff(x,chan,M2,xx2,tu)*ColinearDiff(x,chan,M2)*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));
      
    }

    else if(mel==3) //Exp
    {
      int tu=0;
      complex<double> *xx;
      xx=new complex<double>[4*2*2];
      complex<double> *xx2;
      xx2=new complex<double>[4*2*2];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[2][2]={(0.,0.)};

      ppart[0][0]=pow(M2,0.5)/2.;
      ppart[1][0]=pow(M2,0.5)/2.;
      ppart[0][3]=pow(M2,0.5)/2.;
      ppart[1][3]=-pow(M2,0.5)/2.;

      ppart[2][0]=pow(M2,0.5)/2.;
      ppart[3][0]=pow(M2,0.5)/2.;

      ppart[2][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][3]=pow(M2,0.5)/2.*Xbet;
      ppart[3][3]=-pow(M2,0.5)/2.*Xbet;
      
      //ml5_0_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);

      res=std::imag((TraceHSExpDiff(x,chan,M2,xx,tu)+(1.+ColinearExpDiff(x,chan,M2))*TraceBornDiff(x,chan,M2))*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));

      ppart[3][2]=pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[2][2]=-pow(M2,0.5)/2.*bet*Xsin(x[1]);
      ppart[3][3]=pow(M2,0.5)/2.*Xbet;
      ppart[2][3]=-pow(M2,0.5)/2.*Xbet;
      
      //ml5_0_sloopmatrix_thres_(ppart,xx2,prec_ask,prec_f,ret_code);
      tu=1;
      
      res+=std::imag((TraceHSExpDiff(x,chan,M2,xx2,tu)+(1.+ColinearExpDiff(x,chan,M2))*TraceBornDiff(x,chan,M2))*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));
    }

  
  if(isfinite(res)==0){
    std::cout << "!! Result non finite !!" << std::endl;
    exit(1);
  }
  return res;
}
