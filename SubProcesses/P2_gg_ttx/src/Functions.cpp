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
// -------- Classes ----------------------------------- //
#include "process.h"	// Process//
#include "Constants.h"
// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi;
extern double A[8][8];
extern int Nc, Nf;

double Ci(int chan) //Color coefficient in soft functions, channel 0 is qqbar -> CF, channel 1 is gg -> CA=Nc
{
  double res=0;
  if(chan==0) res=CF;
  else if(chan==1) res=double(Nc);
  else std::cout <<  " Improper channel, select 0 for qqbar and 1 for gg " << std::endl; 
  return res;
}

double B1(int chan) //Color coefficient in soft functions, channel 0 is qqbar -> CF, channel 1 is gg -> CA=Nc
{
  double res=0;
  if(chan==0) res=-3.*CF;
  else
    {
      if(chan==1) res=-2.*beta0;
      else std::cout <<  " Improper channel, select 0 for qqbar and 1 for gg " << std::endl; 
    }

  return res;
}


std::complex<double> N(double x)
{
  std::complex<double> res=0., i(0.0,1.0);
  //res=C+(cos(Phi)+i*sin(Phi))*x/(1.-x);
  //res=C+(cos(Phi)+i*sin(Phi))*pow(x/(1.-x),2.);
  res=C+(cos(Phi)+i*sin(Phi))*(-log(x));
  res+=1.;// N+1
  return res;

}

double Xcos(double x)
{
  double res=0.;
  res=(tanh(20.*x-5.)+tanh(20.*x-15.))/2.;
  return res;
}

double Xsin(double x) // integration from 0 to pi so sin(x) >= 0
{
  double res=0.;
  res=pow(1.-pow(Xcos(x),2.),0.5);
  return res;
}

