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
#include <Eigen/Eigenvalues> 

// -------- Classes ----------------------------------- //
#include "process.h"	// Process//
// -------- Functions --------------------------------- //
// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, CF;
extern int Nc, Nf;

double Sqrt(double);
double Xcos(double);
std::complex<double> N(double);
std::complex<double> Psi(std::complex<double>);

std::complex<double> GammGlobal(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0;
  const int Len=chan+2;
  std::complex<double> i(0.0,1.0), Lbet=0, G[Len][Len]={0.};
  double Xbet=0, bet=0.;

  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI);
  Xbet=Xcos(x[1])*bet;

  // Carefull transpose compared to mathematica code to match the convention in arXiv:0805.1885  
  if(chan==0)
    {
      G[0][0]=-CF*(1.+Lbet);
      G[1][0]=2.*log((1.-Xbet)/(1.+Xbet)); 
      G[0][1]=G[1][0]*CF/2./double(Nc);
      G[1][1]=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.-i*M_PI+log(M2*pow(1.-Xbet,2)/4./pow(Mt,2)))+2./double(Nc)*log((1.+Xbet)/(1.-Xbet));
      
    }
  
  
  if(chan==1)
    {
      G[0][0]=-CF*(1.+Lbet);
      G[0][2]=log((1.-Xbet)/(1.+Xbet));
      G[1][1]=(1.+Lbet+pow(Nc,2)*(-1.-i*M_PI))/2./double(Nc)+double(Nc)*log(M2*(1.-pow(Xbet,2))/4./pow(Mt,2))/2.;
      G[1][2]=double(Nc)/2.*G[0][2];
      G[2][0]=2.*G[0][2];
      G[2][1]=G[0][2]*(pow(Nc,2)-4.)/double(Nc)/2.;
      G[2][2]=G[1][1];
    }

  res=G[m][n];
  return res;  
}


void EigenvGlobal(double *x, double M2, int chan, std::complex<double> **Rdag, std::complex<double> **Rmdag, std::complex<double> **R, std::complex<double> **Rm, std::complex<double> *Lamb)
{
  const int Len=chan+2;
  Eigen::MatrixXcf Ag(Len,Len);
  for(int m=0; m < Len; m++)
    {
      for(int j=0; j < Len; j++) Ag(m,j)=GammGlobal(x,m,j,chan,M2);
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
	  else if(n==2) res=double(pow(Nc,2)-1)*(double(Nc*Nc)-4.)/2./double(Nc);
	}
    }
  
  return res;
}

double S1(double *x, int m, int n, int chan, double M2)
{
  double res=0.;

  double t=0, u=0;
  std::complex<double> Lbet=0, i(0.0,1.0);
  double Xbet=0, bet=0.;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);

  Xbet=Xcos(x[1])*bet;
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI);
  t=-M2/2.*(1-Xbet);
  u=-M2/2.*(1+Xbet);

  
  if(chan==0)
    {
      if(n==0 && m==0) res=std::real(log((1.+bet)/(1.-bet))/bet+M_PI*M_PI/6.+log(muR*muR/M2)+0.5*pow(log(muR*muR/M2),2)+log(muR*muR/M2)*Lbet)*2.*CF;
 
      else if(n==1 && m==1) res=std::real(2.*CF*(log((1.+bet)/(1.-bet))/bet+log(muR*muR/M2))+CF*(M_PI*M_PI/3.+pow(log(muR*muR/M2),2))-Lbet*log(muR*muR/M2)*(double(Nc)-2.*CF)-double(Nc)/2.*pow(log((1.+bet)/(1.-bet)),2)+2.*log(muR*muR/M2)*(log(t/u)*(2.*double(Nc)-4.*CF)-double(Nc)*log(-t/Mt/pow(M2,0.5)))+4.*(2.*CF-double(Nc))*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet)))+2.*(double(Nc)-4.*CF)*(gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet))));

      else if(m==0 && n==1) res=4.*log(u/t)*log(muR*muR/M2)+4.*std::real(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==1 && n==0) res=CF/double(Nc)*2.*log(u/t)*log(muR*muR/M2)+CF/double(Nc)*2.*std::real(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));
      
    }
  

  else if (chan==1)
    {
      if(n==0 && m==0) res=std::real((M_PI*M_PI/3.+pow(log(muR*muR/M2),2))*double(Nc)+2.*CF*(log((1.+bet)/(1.-bet))/bet+log(muR*muR/M2)*(Lbet+1.)));

      else if(n==2 && m==0) res=4.*log(u/t)*log(muR*muR/M2)+4.*std::real(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==1 && n==1) res=2.*CF*log(muR*muR/M2)+std::real(-Lbet*log(muR*muR/M2)*(double(Nc)-2.*CF)+2.*CF*log((1.+bet)/(1.-bet))/bet+double(Nc)*(M_PI*M_PI/3.+pow(log(muR*muR/M2),2)-0.5*pow(log((1.+bet)/(1.-bet)),2)-log(muR*muR/M2)*log(t*u/M2/Mt/Mt)-(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet)))));
      
      else if(m==1 && n==2) res=(3./2.*double(Nc)-4.*CF)*2.*std::real((log(muR*muR/M2)*log(t/u)-gsl_sf_dilog(1.+M2/2./u*(1.+bet))-gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet))));
  
      else if(m==2 && n==0) res=2*log(u/t)*log(muR*muR/M2)+2.*std::real(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==2 && n==1) res=-double(Nc)*(log(muR*muR/M2)*log(t/u)+std::real(gsl_sf_dilog(1.+M2*(1.+bet)/2./t)-gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2*(1.-bet)/t/2.)-gsl_sf_dilog(1.+M2/2.*(1.-bet)/u)));

      else if(m==2 && n==2) res=2.*CF*log(muR*muR/M2)+std::real(-Lbet*log(muR*muR/M2)*(double(Nc)-2.*CF)+2.*CF*log((1.+bet)/(1.-bet))/bet+double(Nc)*(M_PI*M_PI/3.+pow(log(muR*muR/M2),2)-0.5*pow(log((1.+bet)/(1.-bet)),2)-log(muR*muR/M2)*log(t*u/M2/Mt/Mt)-(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet)))));
    
    
    }
  return res/double(Nc)*alpha_s/2./M_PI;
}

double S1t(double *x, int m, int n, int chan, double M2)
{
  double res=0.;
  int Len=chan+1;

  for(int k=0; k < Len; k++)
    {
      res+=S0(m,k,chan)*S1(x,k,n,chan,M2)*double(Nc);
    }
  
  return res;
}


std::complex<double> SG(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0.;
  int Len=chan+1;
  std::complex<double> Nb=N(x[0])*exp(-Psi(1.));

  for(int k=0; k < Len; k++)
    {
      res+=S0(m,k,chan)*GammGlobal(x,k,n,chan,M2)+std::conj(GammGlobal(x,k,m,chan,M2))*S0(k,n,chan);
    }

  return res*(-alpha_s/2./M_PI*log(Nb));
}


void DisplayS1(double *x, int chan, double M2)
{

  int Len=chan+1;
  double t=0, u=0;
  std::complex<double> Lbet=0, i(0.0,1.0);
  double Xbet=0, bet=0.;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);

  Xbet=Xcos(x[1])*bet;
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI);
  t=-M2/2.*(1-Xbet);
  u=-M2/2.*(1+Xbet);

  for(int j=0; j < Len+1; j++)
    {
      for(int k=0; k< Len+1; k++)
	{
	  std::cout << " S1 qq " << j << k << " = " << S1(x,j,k,chan,M2) << " Lbet = " << Lbet << " bet = " << bet << " t1 = " << t << " u1 = " << u <<" muR = " << muR << " M2 = " << M2 << std::endl;

	}
    }
  
}


void DisplayG1(double *x, int chan, double M2)
{

  int Len=chan+1;
  std::complex<double> Lbet=0, i(0.0,1.0);
  double Xbet=0, bet=0.;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);

  Xbet=Xcos(x[1])*bet;
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI);
  
  for(int j=0; j < Len+1; j++)
    {
      for(int k=0; k< Len+1; k++)
	{
	  std::cout << " G1 qq " << j << k << " = " << GammGlobal(x,j,k,chan,M2) << " Lbet = " << Lbet << " Xbet = " << Xbet  << " M2 = " << M2 << std::endl;
	}
    }
}
