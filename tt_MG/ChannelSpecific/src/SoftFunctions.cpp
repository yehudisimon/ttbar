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
  // Gamma not symmetric but real
  
  std::complex<double> res=0;
  const int Len=chan+2;
  std::complex<double> i(0.0,1.0), Lbet=0;
  double Logb=0.;
  double Xbet=0, bet=0.;

  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI);
  //Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet)));
  Xbet=Xcos(x[1])*bet;
  Logb=log((1.-Xbet)/(1.+Xbet));

  if(chan==0)
    {
      if(n==0 and m==0) res=-CF*(1.+Lbet); 
      else if(m==0 and n==1) res=Logb*CF/double(Nc);
      else if(m==1 and n==0) res=2.*Logb;
      else if(m==1 and n== 1) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.+i*M_PI)-2.*Logb/double(Nc)+double(Nc)/2.*log(M2/4./pow(Mt,2)*pow(1.-Xbet,2));
    }
  
  
  if(chan==1)
    {
      if(n==0 and m==0) res=-CF*(1.+Lbet); 
      else if(m==0 and n==1) res=-Logb;
      else if(m==1 and n==0) res=-2.*Logb;
      else if(m==1 and n== 1) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.+i*M_PI)+double(Nc)/4.*2.*log(M2/4./pow(Mt,2)*(1.-pow(Xbet,2)));
      else if(m==1 and n==2) res= -Logb/2./double(Nc)*(double(Nc*Nc)-4.);
      else if(m==2 and n==1) res= -Logb*double(Nc)/2.;
      else if(m==2 and n==2) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.+i*M_PI)+double(Nc)/4.*2.*log(M2/4./pow(Mt,2)*(1.-pow(Xbet,2)));
    }

  
  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in Gamma !!" << std::endl;
    exit(1);
  }
  
  return res; //Should be multiplied by alpha_s/M_PI to recover total anomalous dimension like in 2109.15039
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
	  //std::cout << " Chan = " << chan << " M2 = " << M2 << " muR = " << muR << " muF = " << muF << " Xbet = " << Xcos(x[1])*pow(1.-4.*Mt*Mt/M2,0.5) << "Gamma = " << GammGlobal(x,m,j,chan,M2) << " m = " << m << " n = " << j << std::endl;
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
	  //std::cout << " Rdag = " << Rdag[n][i] << " R = " << R[n][i] << std::endl;
	}
      //std::cout << "Lamb = " << Lamb[n] << " n = " << n << std::endl;  
    }
  //exit(3);
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

std::complex<double> S1old(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0.;
  double t=0, u=0;
  std::complex<double> Lbet=0., i(0.0,1.0);
  double Xbet=0, bet=0.;
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet)));
  t=-M2/2.*(1.-Xbet);
  u=-M2/2.*(1.+Xbet);
  
    if(chan==0)
    {
      if(n==0)
	{
	  if(m==0) res=(log((1.+bet)/(1.-bet))/bet+M_PI*M_PI/6.+log(muR*muR/M2)*(1.+Lbet)+0.5*pow(log(muR*muR/M2),2))*2.*CF+CF*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)));
	  
	  else
	    {
	      res=4.*(log(u/t)*log(muR*muR/M2)+gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));
	    }
	}
      else
	{
	  if(m==1) res=2.*CF*(log((1.+bet)/(1.-bet))/bet+log(muR*muR/M2))+CF*(M_PI*M_PI/3.+pow(log(muR*muR/M2),2))-Lbet*log(muR*muR/M2)*(double(Nc)-2.*CF)-double(Nc)/2.*pow(log((1.+bet)/(1.-bet)),2)+2.*log(muR*muR/M2)*(log(t/u)*(2.*double(Nc)-4.*CF)-double(Nc)*log(-t/Mt/pow(M2,0.5)))+4.*(2.*CF-double(Nc))*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet)))+2.*(double(Nc)-4.*CF)*(gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet)))+(CF-double(Nc)/2.)*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)));

	  else
	    {
	      res=CF/double(Nc)*2.*(log(u/t)*log(muR*muR/M2)+gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));
	    }
	}	  
    }

    else if (chan==1)
    {
      if(n==0 && m==0) res=(M_PI*M_PI/3.+pow(log(muR*muR/M2),2))*double(Nc)+2.*CF*(log((1.+bet)/(1.-bet))/bet+log(muR*muR/M2)*(Lbet+1.))+CF*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)));

      else if(n==0 && m==1) res=4.*log(u/t)*log(muR*muR/M2)+4.*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==1 && n==1) res=2.*CF*log(muR*muR/M2)-Lbet*log(muR*muR/M2)*(double(Nc)-2.*CF)+2.*CF*log((1.+bet)/(1.-bet))/bet+double(Nc)*(M_PI*M_PI/3.+pow(log(muR*muR/M2),2)-0.5*pow(log((1.+bet)/(1.-bet)),2)-log(muR*muR/M2)*log(t*u/M2/Mt/Mt)-(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet))))-(double(Nc)/2.-CF)*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)));

      else if(m==2 && n==1) res=double(Nc)*(log(u/t)*log(muR*muR/M2)+gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==0 && n==1) res=2*log(u/t)*log(muR*muR/M2)+2.*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==1 && n==2) res=(double(Nc*Nc)-4.)/double(Nc)*(log(u/t)*log(muR*muR/M2)+gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));

      else if(m==2 && n==2) res=2.*CF*log(muR*muR/M2)-Lbet*log(muR*muR/M2)*(double(Nc)-2.*CF)+2.*CF*log((1.+bet)/(1.-bet))/bet+double(Nc)*(M_PI*M_PI/3.+pow(log(muR*muR/M2),2)-0.5*pow(log((1.+bet)/(1.-bet)),2)-log(muR*muR/M2)*log(t*u/M2/Mt/Mt)-(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet))))-(double(Nc)/2.-CF)*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)));

    }
    return res*alpha_s/2./M_PI;       
}

std::complex<double> S1(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0.;
    int tmp; // S1 "hermitian" (symmetric in fact)
  if(m>n)
    {
      tmp=m;
      m=n;
      n=tmp;
    }
  
  double t=0, u=0;
  std::complex<double> Lbet=0., i(0.0,1.0), Logb=0.;
  double Xbet=0, bet=0.;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);

  Xbet=Xcos(x[1])*bet;
  Lbet=(1.+bet*bet)/2./bet*log((1.-bet)/(1.+bet));
  //Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI);
  Logb=log((1.-Xbet)/(1.+Xbet));
  t=-M2/2.*(1.-Xbet);
  u=-M2/2.*(1.+Xbet);
  
  if(chan==0)
    {
      if(n==0) res=double(Nc)*(double(Nc*Nc)-1.)*(1./bet*log((1.-bet)/(1.+bet))+pow(M_PI,2)/6.+(Lbet+1.)*log(muR*muR/M2)+0.5*pow(log(muR*muR/M2),2)); //m=0
      
      else
	{
	  if(m==0) res=(double(Nc*Nc)-1.)*(log(muR*muR/M2)*log(u/t)+gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));
	  else res=CF*(CF*double(Nc)*(pow(log(muR*muR/M2),2)/2.+pow(M_PI,2)/6.+1./bet*log((1.-bet)/(1.+bet)))-(1.+Lbet)*log(muR*muR/M2)/2.+2.*log(muR*muR/M2)*log(t/u)+double(Nc*Nc)*(log(muR*muR/M2)/2.-pow(log((1.+bet)/(1.-bet)),2)/4.-log(muR*muR/M2)*log(t/Mt/pow(M2,0.5))-gsl_sf_dilog(1.+M2*(1.-bet)/2./t)-gsl_sf_dilog(1.+M2/2.*(1.+bet)))-2.*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)))); //m=1
	}
    }

    
  

  else if (chan==1)
    {
      if(n==0) res=(double(Nc*Nc)-1.)*((pow(M_PI,2)/3.+pow(log(muR*muR/M2),2))*pow(Nc,2)+(double(Nc*Nc)-1.)*(log(muR*muR/M2)*(1.+Lbet)+log((1.+bet)/(1.-bet))/bet)); //m=0

      else if(n==1)
	{
	  if(m==0) res=2.*double(Nc)*(-1.+double(Nc*Nc))*(log(u/t)*log(muR*muR/M2)+gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));
	  else res=(double(Nc*Nc)-1.)/2.*((double(Nc*Nc)-1.)*log((1.+bet)/(1.-bet))/bet-log(muR*muR/M2)*(1.+Lbet)+double(Nc*Nc)*(pow(M_PI,2)/3.+log(muR*muR/M2)+pow(log(muR*muR/M2),2)-1./2.*pow(log((1.+bet)/(1.-bet)),2)-log(t*u/M2/Mt/Mt)*log(muR*muR/M2)-gsl_sf_dilog(1.+M2/2./u*(1.+bet))-gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet))));
	}
      else if(n==2)
	{
	  if(m==0) res=0.;
	  else if(m==1) res=(double(Nc*Nc)-4.)/2.*(-1.+double(Nc*Nc))*(log(u/t)*log(muR*muR/M2)+gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)));
	  else res=(double(Nc*Nc)-4.)/double(Nc*Nc)*(double(Nc*Nc)-1.)/2.*((double(Nc*Nc)-1.)*log((1.+bet)/(1.-bet))/bet-log(muR*muR/M2)*(1.+Lbet)+double(Nc*Nc)*(pow(M_PI,2)/3.+log(muR*muR/M2)+pow(log(muR*muR/M2),2)-1./2.*pow(log((1.+bet)/(1.-bet)),2)-log(t*u/M2/Mt/Mt)*log(muR*muR/M2)-gsl_sf_dilog(1.+M2/2./u*(1.+bet))-gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet))));
	}    
    }

  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in S1 !!" << std::endl;
    std::cout << "chan = " << chan << " i = " << m << " j = " << n << std::endl;
    exit(1);
  }

  return res*alpha_s/2./M_PI;
}

std::complex<double> S1t(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0., un(1.,0.);
  int Len=chan+2;

  for(int k=0; k < Len; k++)
    {
      res+=S0(m,k,chan)*S1old(x,k,n,chan,M2); //*Nc in old implementation
    }
  return res;
}


std::complex<double> SG(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0.;
  int Len=chan+2;
  std::complex<double> Nb=N(x[0])*exp(-Psi(1.));

  for(int k=0; k < Len; k++)
    {
      res+=S0(m,k,chan)*GammGlobal(x,k,n,chan,M2)+std::conj(GammGlobal(x,k,m,chan,M2))*S0(k,n,chan);
      //res+=GammGlobal(x,m,k,chan,M2)*S0(k,n,chan)+std::conj(GammGlobal(x,n,k,chan,M2))*S0(m,k,chan);
      //std::cout << " Gamma_km = " << GammGlobal(x,k,m,chan,M2) << " Gamma_mk ="  << GammGlobal(x,m,k,chan,M2) << std::endl;
    }
  
  //exit(5);
  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in SG !!" << std::endl;
    exit(1);
  }
  return res*(-alpha_s/M_PI*log(Nb)); // With covention of literature log(1-2*lambda)*Gamma(without alpha)/beta0
}


void DisplayS1(double *x, int chan, double M2)
{

  int Len=chan+2;
  double t=0, u=0;
  std::complex<double> Lbet=0, i(0.0,1.0);
  double Xbet=0, bet=0.;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);

  Xbet=Xcos(x[1])*bet;
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI);
  t=-M2/2.*(1-Xbet);
  u=-M2/2.*(1+Xbet);

  for(int j=0; j < Len; j++)
    {
      for(int k=0; k< Len; k++)
	{
	  std::cout << " S1 qq " << j << k << " = " << S1(x,j,k,chan,M2) << " Lbet = " << Lbet << " bet = " << bet << " t1 = " << t << " u1 = " << u <<" muR = " << muR << " M2 = " << M2 << std::endl;

	}
    }
  
}


void DisplayG1(double *x, int chan, double M2)
{

  int Len=chan+2;
  std::complex<double> Lbet=0, i(0.0,1.0);
  double Xbet=0, bet=0.;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);

  Xbet=Xcos(x[1])*bet;
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI);
  
  for(int j=0; j < Len; j++)
    {
      for(int k=0; k< Len; k++)
	{
	  std::cout << " G1 qq " << j << k << " = " << GammGlobal(x,j,k,chan,M2) << " Lbet = " << Lbet << " Xbet = " << Xbet  << " M2 = " << M2 << std::endl;
	}
    }
}
