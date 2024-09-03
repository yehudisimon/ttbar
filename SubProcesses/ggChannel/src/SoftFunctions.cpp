// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <iostream>
#include <iomanip>
#include <cstddef>
#include <math.h>
#include <complex>
#include <boost/math/special_functions/zeta.hpp>
#include <gsl/gsl_sf_dilog.h>
#include <Eigen/Eigenvalues> 

// -------- Classes ----------------------------------- //
#include "process.h"	// Process//
// -------- Functions --------------------------------- //
// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, CF, beta0;
extern int Nc, Nf;

double Sqrt(double);
double Xcos(double);
double Ci(int);
std::complex<double> N(double);
std::complex<double> Psi(std::complex<double>);

std::complex<double> GammGlobal(double *x, int m, int n, int chan, double M2)
{
  // Gamma not symmetric, only complex in diagonal
  
  std::complex<double> res=0;
  const int Len=chan+2;
  std::complex<double> i(0.0,1.0), Lbet=0;
  double Logb=0.;
  double Xbet=0, bet=0.;

  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))-i*M_PI); //HS
  Xbet=Xcos(x[1])*bet;
  Logb=log((1.-Xbet)/(1.+Xbet));
  
  if(chan==0)
    {
      exit(5);
    }
  
  
  if(chan==1)
    {
      if(n==0 and m==0) res=-CF*(1.+Lbet); 
      else if(m==1 and n==0) res=Logb;
      else if(m==0 and n==1) res=2.*Logb;
      else if(m==1 and n== 1) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.-i*M_PI)+double(Nc)/4.*2.*log(M2/4./pow(Mt,2)*(1.-pow(Xbet,2)));
      else if(m==2 and n==1) res= Logb/2./double(Nc)*(double(Nc*Nc)-4.);
      else if(m==1 and n==2) res= Logb*double(Nc)/2.;
      else if(m==2 and n==2) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.-i*M_PI)+double(Nc)/4.*2.*log(M2/4./pow(Mt,2)*(1.-pow(Xbet,2)));
    }

  
  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in Gamma !!" << std::endl;
    exit(1);
  }

  return res*(-2.); //Should be multiplied by alpha_s/2./M_PI to recover total anomalous dimension (not like in 2109.15039)
}


void EigenvGlobal(double *x, double M2, int chan, std::complex<double> **Rdag, std::complex<double> **Rmdag, std::complex<double> **R, std::complex<double> **Rm, std::complex<double> *Lamb)
{
  const int Len=chan+2;
  Eigen::MatrixXcf Ag(Len,Len);
  for(int m=0; m < Len; m++)
    {
      for(int j=0; j < Len; j++)
	{
	  Ag(m,j)=GammGlobal(x,m,j,chan,M2);
	}
    }
  
  
  Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
  ces.compute(Ag);

  for(int n=0; n < Len; n++)
    {
      Lamb[n]=ces.eigenvalues()(n);
      for(int i=0; i < Len; i++)
	{
	  R[n][i]=ces.eigenvectors().inverse()(n,i);
	  Rm[n][i]=ces.eigenvectors()(n,i);
	  Rdag[n][i]=ces.eigenvectors().inverse().transpose().conjugate()(n,i);
	  Rmdag[n][i]=ces.eigenvectors().transpose().conjugate()(n,i);
	}  
    }
}

double S0(int m, int n, int chan)
{
  double res=0.;

  if(m==n)
    {
      if(chan==0)
	{
	  //Color basis = {delta_(q,qx)*delta_(t,tx),T^i_(qb,q)*T^i_(t,tx)}
	  if(n==0) res= double(Nc)*double(Nc);
	  else if(n==1) res=CF/2.*double(Nc);
	}
    

      else if (chan==1)
	{
	  //Color basis = {delta^(g1,g2)*delta_(t,tx),i*f^(g1,g2,i)*T^i_(t,tx),d^(g1,g2,i)*T^i_(t,tx)}
	  if(n==0) res= double(Nc)*(pow(double(Nc),2)-1.);
	  else if(n==1) res=double(Nc)*(pow(double(Nc),2)-1.)/2.;
	  else if(n==2) res=(pow(double(Nc),2)-1.)*(double(Nc*Nc)-4.)/2./double(Nc);
	}
    }
  
  return res;
}

double g3Nindep(int chan, double M2)
{
  double res=0.;

  res=Ci(chan)*4.; //Gamma_cusp^0                     

  return res*alpha_s/2./M_PI*pow(log(muR*muR/M2),2)/8.*2.; //*2 for each incomming particule  
}


double Beam(int chan,double M2)
{
  double res=0.;
  double gamma_col=0.;

  if(chan==0)
    {
      gamma_col=3.*CF; // 2*gamma_coll^(1) of resummino 1304.0790
    }
  else
    {
      gamma_col=beta0*2.; // 2*gamma_coll^(1) of resummino 1304.0790
    }
  
  res=-2.*gamma_col; // *2 is for 2 incomming particles -1. from HS's notes  
  
  return res*alpha_s/4./M_PI*log(muF/muR)*2.; //*2 for mu^2 
}

std::complex<double> Gamma_sNindep(double *x,int m, int n, int chan, double M2)
{
  std::complex<double> res=0.;
  int Len=chan+2;
  for(int k=0; k < Len; k++)
    {
      res+=S0(m,k,chan)*std::conj(GammGlobal(x,n,k,chan,M2))+GammGlobal(x,m,k,chan,M2)*S0(k,n,chan); //HS
    }
                                                                                      
  return res*alpha_s/2./M_PI*log(muR*muR/M2)/2.; //HS Convention Gamma=alpha_s/2.*M_PI*Gamma^1
}


std::complex<double> SG(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0.;
  int Len=chan+2;
  std::complex<double> Nb=(N(x[0]))*exp(-Psi(1.));

  for(int k=0; k < Len; k++)
    {   
      res+=S0(m,k,chan)*std::conj(GammGlobal(x,n,k,chan,M2))+GammGlobal(x,m,k,chan,M2)*S0(k,n,chan); //HS
    }
  
  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in SG !!" << std::endl;
    exit(1);
  }
  return res*alpha_s/2./M_PI*log(Nb); // With covention of literature log(1-2*lambda)*Gamma^1(without alpha)/2*beta0 factor 2 in addition by convention alpha_s/M_PI HS
}

std::complex<double> S1old(double *x, int m, int n, int chan, double M2) // S1
{
  std::complex<double> res=0.;
  double t=0, u=0;
  std::complex<double> Lbet=0., i(0.0,1.0);
  double Xbet=0, bet=0.;
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet)));
  t=-M2/2.*(1.-Xbet); //t1
  u=-M2/2.*(1.+Xbet); //u1
  
    if(chan==0)
    {
      exit(1);
    }

    else if (chan==1)
    {
      if(n==0 && m==0) res=M_PI*M_PI/3.*double(Nc)+2.*CF*log((1.+bet)/(1.-bet))/bet+CF*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)));

      else if(n==1 && m==0) res=4.*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==1 && n==1) res=2.*CF*log((1.+bet)/(1.-bet))/bet+double(Nc)*(M_PI*M_PI/3.-0.5*pow(log((1.+bet)/(1.-bet)),2))-double(Nc)*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet)))-(double(Nc)/2.-CF)*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)));

      else if(m==1 && n==2) res=double(Nc)*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==1 && n==0) res=2.*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==2 && n==1) res=(double(Nc*Nc)-4.)/double(Nc)*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==2 && n==2) res=2.*CF*log((1.+bet)/(1.-bet))/bet+double(Nc)*(M_PI*M_PI/3.-0.5*pow(log((1.+bet)/(1.-bet)),2))-double(Nc)*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet)))-(double(Nc)/2.-CF)*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)));

    }

    return res*alpha_s/2./M_PI;
}


std::complex<double> S1t(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0., res2=0., un(1.,0.), g1=0.,g2=0.;
  int Len=chan+2;

  for(int k=0; k < Len; k++)
    {
      res+=S1old(x,m,k,chan,M2)*S0(k,n,chan); // S1(muS)
    }
  
  res+=Gamma_sNindep(x,m,n,chan,M2);
  res+=S0(m,n,chan)*g3Nindep(chan,M2);
  res+=S0(m,n,chan)*Beam(chan,M2);
  return res;
}
