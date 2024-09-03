// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <iostream>
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
// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi;
extern int Nc, npart;

extern "C"{
  //void gg2ttx_shfmatrix_(int&, int&, double&, double&, double&, double&,complex<double>(*));
  //void ml5_0_sloopmatrix_thres_(double[4][4],complex<double>(*),double&,double[2][2],int&);
  //void update_mur_(double&);
}

double Xcos(double);
double Xsin(double);

double Sqrt(double x)
{
  return pow(x,0.5);
}

std::complex<double> C0(double a,double b, double c, double m1, double m2, double m3) // C0 =-i*M_PI^2/2./Mt^2* lim_{eps->0} 3F2({1,1,1},{3/2,2},(M2+i*eps)/(4.*Mt^2))=-i*M_PI^2/2./Mt^2*lim_{eps->0} arcsin(sqrt(M2+i*eps)/2./Mt)^2/(M2/4./Mt^2) !!! Careful z=M2/4./Mt^2 > 1 here so no log formula to express arcin -> use gsl_complex function to evaluate. Careful of branch cut, here M2 as an infinitesimal imaginary part which changes the sign of imaginary part. See 0712.1851 for definition of the integral: C0({x})=i*M_PI^2*I_3({x})
{
  gsl_complex resgsl=gsl_complex_pow_real(gsl_complex_arcsin_real(Sqrt(c/4./m1)),2.);
  std::complex<double> i(0,1);
  resgsl=gsl_complex_mul_real(resgsl,4.*m1/c);
  std::complex<double> res=GSL_REAL(resgsl)-i*GSL_IMAG(resgsl);
  res*=-i*pow(M_PI,2.)/2./pow(Mt,2.);
  return res;
}


// *********************************************************************** //
//  Hard coefficients //
// *********************************************************************** //


double H0(double *x, int m, int n, int chan, double M2)
{
  double res=0.;
  std::complex<double> i(0.0,1.0);
  double Xbet=0, bet=0., t=0., u=0., dR=0.;
  int tmp=0, Len=chan+2, Nsym=1;

  //H0 symmetric
  if(m>n)
    {
      tmp=m;
      m=n;
      n=tmp;
    }
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;

  t=-M2/2.*(1-Xbet);
  u=-M2/2.*(1+Xbet);


  if(chan==0)
    {
      dR=double(Nc)*2.; //         
    }
  else if(chan==1)
    {
      dR=2.*(pow(Nc,2)-1.);// CDR convention          
      Nsym=1;
    }


  if(chan==0)
    {
      if(m*n==1)
        {
          res=2.*((pow(t,2)+pow(u,2))/pow(M2,2)+2.*Mt*Mt/M2); //see 0805.1885 or SCET 1003.5827
	}
    }

  else if(chan==1)
    {
      double hgg=0.;
      hgg=pow(M2,2)/2./t/u*(pow(t,2)/pow(M2,2)+pow(u,2)/pow(M2,2)+4.*pow(Mt,2)/M2-4.*pow(Mt,4)/t/u);
      res=hgg;
      if(m==0)
        {
          res*=1./double(Nc);
          if(n==0) res*=1./double(Nc);
          else if(n==1) res*=(t-u)/M2;
	}

      else if(m==1)
        {
          res*=(t-u)/M2;
          if(n==1) res*=(t-u)/M2;;
        }
    }
  
  return res/2./M2/dR/dR/Nsym*pow(alpha_s*4.*M_PI,2)*4.;  //Global factor with respect to SCET because of conventions M_ren -> M_ren^(0) and H^(0)_SCET = H^(0)/4, no 3/8 nor 8/3 here         
}

double H0Phase(double *ppart, int m, int n, int chan, double M2)
{
  double res=0.;
  std::complex<double> i(0.0,1.0);
  double Xbet=0, bet=0., t=0., u=0., dR=0.;
  int tmp=0, Len=chan+2, Nsym=1;
  
  //H0 symmetric
  if(m>n)
    {
      tmp=m;
      m=n;
      n=tmp;
    }
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  int sign=1;
  for(int l=0; l<4; l++)
    {
      if(l>0) sign=-1;
      t+=sign*pow(ppart[0+(npart+2)*l]-ppart[2+(npart+2)*l],2);
      u+=sign*pow(ppart[0+(npart+2)*l]-ppart[3+(npart+2)*l],2);
    }

  //std::cout << " t = " << t << std::endl;
  
  t-=pow(Mt,2); //t -> t1
  u-=pow(Mt,2); //u -> u1
  
  Xbet=1.+2.*t/M2;

  if(chan==0)
    {
      dR=double(Nc)*2.; //         
    }
  else if(chan==1)
    {
      dR=2.*(pow(Nc,2)-1.);// CDR convention          
      Nsym=1;
    }


  if(chan==0)
    {
      if(m*n==1)
        {
          res=2.*((pow(t,2)+pow(u,2))/pow(M2,2)+2.*Mt*Mt/M2); //see 0805.1885 or SCET 1003.5827
	}
    }

  else if(chan==1)
    {
      double hgg=0.;
      hgg=pow(M2,2)/2./t/u*(pow(t,2)/pow(M2,2)+pow(u,2)/pow(M2,2)+4.*pow(Mt,2)/M2-4.*pow(Mt,4)/t/u);
      res=hgg;
      if(m==0)
        {
          res*=1./double(Nc);
          if(n==0) res*=1./double(Nc);
          else if(n==1) res*=(t-u)/M2;
	}

      else if(m==1)
        {
          res*=(t-u)/M2;
          if(n==1) res*=(t-u)/M2;;
        }
    }
  //std::cout << " H0 = " << res/2./M2 << " m = " << m << " n = " << n << " alpha_s = " << alpha_s << " M2 = " << M2 << std::endl;
  
  return res/2./M2/dR/dR/Nsym*pow(alpha_s*4.*M_PI,2)*4.;  //Global factor with respect to SCET because of conventions M_ren -> M_ren^(0) and H^(0)_SCET = H^(0)/4, no 3/8 nor 8/3 here         
}
