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
std::complex<double> MellinPDFPhase(double,int,double);
std::complex<double> GlobalPhase(double,double,double);
std::complex<double> TracePhase(double*,int,double,complex<double>*,double);
std::complex<double> TraceBornPhase(double*,int,double);
std::complex<double> TraceExpPhase(double*, int , double , complex<double> *,double);
std::complex<double> ColinearPhase(double, int , double );
std::complex<double> ColinearExpPhase(double, int , double);

struct my_f_params { double sc; int mel; double M2; };  

extern "C"{
  void __genps_nbody_MOD_generate_mom(int&,double&,double*,double*,double*,double&);
  void ml5_2_sloopmatrix_thres_(double*,complex<double>(*),double&,double[3][3],int&);
}

// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi, xmax, xmin;
extern int chan, npart;

//------------------------------------------------------------//
//--------------------Total Integrand-------------------------//
//------------------------------------------------------------//
double TotPhaseSpace(double *x, size_t dim, void *params)
{
  double res=0., bet, t, u, jac=0.;
  std::complex<double> temp=0, i(0.0,1.0), Nb;
  const int Len=chan+2;
  
  struct my_f_params * fparams = (struct my_f_params *)params;  
  double M2=fparams->M2, sc=fparams->sc, Q=pow(M2,0.5);
  int mel=fparams->mel;
  int d=dim-1;
  double *mpart;
  size_t siz=npart;
  size_t sizp=npart+2;  
  mpart=new double[siz];
  for(int k=0; k<siz; k++) mpart[k]=Mt;
  
  double *ppart;
  ppart= new double [4*sizp];

  double *y;
  y=new double[d];
  for(int j=0; j < d; j++)
    {
      y[j]=x[j];
    }
  
  __genps_nbody_MOD_generate_mom(npart,Q,mpart,y,ppart,jac);

  // std::cout << " jac = " << jac << std::endl;
  // for(int i=0; i<npart+2; i++)
  //   {
  //     for(int j=0; j<4; j++) std::cout << " p^t("<< i << "," << j << ") = "  << ppart[i+(npart+2)*j] << std::endl;
  //   }

  double *ppar;
  ppar= new double [4*sizp]; // transpose for MadLoop input 
  for(int i=0; i<npart+2; i++)
    {
      for(int j=0; j<4; j++) ppar[j+i*4]=ppart[i+j*(npart+2)];
    }

  // for(int i=0; i<npart+2; i++)
  //   {
  //     for(int j=0; j<4; j++) std::cout << " p("<< i << "," << j << ") = "  << ppar[j+i*4] << std::endl;
  //   }

  //Nb=N(x[d])*exp(-Psi(1.));
  //std::cout << " Nb = " << Nb << std::endl;
  //bet=pow(1.-4.*Mt*Mt/M2,0.5);
  //t=-M2/2.*(1-Xbet);
  complex<double> *xx;
  xx=new complex<double>[4*Len*Len];
  int ret_code=0;
  double prec_ask=-1;
  double prec_f[3][3]={(0.,0.)};

  if(mel==0)
    {
      std::complex<double> tmp=0.;      
      tmp=TraceBornPhase(ppart,chan,M2);
      res=std::imag(tmp*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);
    }

  else if(mel==1) //Diff
    {
      
      ml5_2_sloopmatrix_thres_(ppar,xx,prec_ask,prec_f,ret_code);

      res+=std::imag(TracePhase(ppart,chan,M2,xx,x[d])*ColinearPhase(x[d],chan,M2)*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);
      //res+=std::imag(TracePhase(ppart,chan,M2,xx,x[d])*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);
      //res+=std::imag(TraceBornPhase(ppart,chan,M2)*ColinearPhase(x[d],chan,M2)*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);

      res-=std::imag((TraceExpPhase(ppart,chan,M2,xx,x[d])+(1.+ColinearExpPhase(x[d],chan,M2))*TraceBornPhase(ppart,chan,M2))*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);
      //res-=std::imag((TraceExpPhase(ppart,chan,M2,xx,x[d])+TraceBornPhase(ppart,chan,M2))*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);
      //res-=std::imag((1.+ColinearExpPhase(x[d],chan,M2))*TraceBornPhase(ppart,chan,M2)*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);
      
    }

    else if(mel==2) //NLL 
    {
      ml5_2_sloopmatrix_thres_(ppar,xx,prec_ask,prec_f,ret_code);

      res=std::imag(TracePhase(ppart,chan,M2,xx,x[d])*ColinearPhase(x[d],chan,M2)*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);
      
    }

    else if(mel==3) //Exp
    {
      ml5_2_sloopmatrix_thres_(ppar,xx,prec_ask,prec_f,ret_code);

      res=std::imag((TraceExpPhase(ppart,chan,M2,xx,x[d])+(1.+ColinearExpPhase(x[d],chan,M2))*TraceBornPhase(ppart,chan,M2))*MellinPDFPhase(x[d],chan,M2)*GlobalPhase(x[d],sc,M2)*jac);

    }

  delete [] ppart;
  delete [] ppar;
  delete [] xx;
  
  if(isfinite(res)==0){
    std::cout << "!! Result non finite !!" << std::endl;
    exit(1);
  }
  return res;
}
