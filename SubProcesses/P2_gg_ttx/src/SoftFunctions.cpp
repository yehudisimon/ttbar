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
extern int Nc, Nf, npart;

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
  //Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI); //Me
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))-i*M_PI); //HS
  Xbet=Xcos(x[1])*bet;
  Logb=log((1.-Xbet)/(1.+Xbet));
  
  if(chan==0)
    {
      exit(5);
    }
  
  
  // if(chan==1) //Me -> *2 at the end
  //   {
  //     if(n==0 and m==0) res=-CF*(1.+Lbet); 
  //     else if(m==0 and n==1) res=Logb;
  //     else if(m==1 and n==0) res=2.*Logb;
  //     else if(m==1 and n== 1) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.+i*M_PI)+double(Nc)/4.*2.*log(M2/4./pow(Mt,2)*(1.-pow(Xbet,2)));
  //     else if(m==1 and n==2) res= Logb/2./double(Nc)*(double(Nc*Nc)-4.);
  //     else if(m==2 and n==1) res= Logb*double(Nc)/2.;
  //     else if(m==2 and n==2) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.+i*M_PI)+double(Nc)/4.*2.*log(M2/4./pow(Mt,2)*(1.-pow(Xbet,2)));
  //   }

  if(chan==1) // HS -> *-2 at the end
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

std::complex<double> GammPhase(double *ppart, int m, int n, int chan, double M2)
{
  // Gamma not symmetric, only complex in diagonal
  
  std::complex<double> res=0;
  const int Len=chan+2;
  std::complex<double> i(0.0,1.0), Lbet=0;
  double Logb=0.;
  double Xbet=0, bet=0., t=0.;

  int sign=1;
  for(int l=0; l<4; l++)
    {
      if(l>0) sign=-1;
      t+=sign*pow(ppart[0+(npart+2)*l]-ppart[2+(npart+2)*l],2);
    }
  t-=pow(Mt,2); //t -> t1
  Xbet=1.+2./M2*t;
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  //Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))+i*M_PI); //Me
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet))-i*M_PI); //HS
  Logb=log((1.-Xbet)/(1.+Xbet));
  
  if(chan==0)
    {
      exit(5);
    }
  
  
  // if(chan==1) //Me -> *2 at the end
  //   {
  //     if(n==0 and m==0) res=-CF*(1.+Lbet); 
  //     else if(m==0 and n==1) res=Logb;
  //     else if(m==1 and n==0) res=2.*Logb;
  //     else if(m==1 and n== 1) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.+i*M_PI)+double(Nc)/4.*2.*log(M2/4./pow(Mt,2)*(1.-pow(Xbet,2)));
  //     else if(m==1 and n==2) res= Logb/2./double(Nc)*(double(Nc*Nc)-4.);
  //     else if(m==2 and n==1) res= Logb*double(Nc)/2.;
  //     else if(m==2 and n==2) res=(1.+Lbet)/2./double(Nc)+double(Nc)/2.*(-1.+i*M_PI)+double(Nc)/4.*2.*log(M2/4./pow(Mt,2)*(1.-pow(Xbet,2)));
  //   }

  if(chan==1) // HS -> *-2 at the end
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

void EigenvPhase(double *ppart, double M2, int chan, std::complex<double> **Rdag, std::complex<double> **Rmdag, std::complex<double> **R, std::complex<double> **Rm, std::complex<double> *Lamb)
{
  const int Len=chan+2;
  Eigen::MatrixXcf Ag(Len,Len);
  for(int m=0; m < Len; m++)
    {
      for(int j=0; j < Len; j++)
	{
	  Ag(m,j)=GammPhase(ppart,m,j,chan,M2);
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
      //res+=S0(m,k,chan)*GammGlobal(x,k,n,chan,M2)+std::conj(GammGlobal(x,k,m,chan,M2))*S0(k,n,chan); //Me
    }
                                                                                      
  return res*alpha_s/2./M_PI*log(muR*muR/M2)/2.; //HS Convention Gamma=alpha_s/2.*M_PI*Gamma^1
  //return -res*alpha_s/2./M_PI*log(muR*muR/M2)/2.; //Me Convention Gamma=alpha_s/2.*M_PI*Gamma^1
}

std::complex<double> Gamma_sNindepPhase(double *ppart,int m, int n, int chan, double M2)
{
  std::complex<double> res=0.;
  int Len=chan+2;
  for(int k=0; k < Len; k++)
    {
      res+=S0(m,k,chan)*std::conj(GammPhase(ppart,n,k,chan,M2))+GammPhase(ppart,m,k,chan,M2)*S0(k,n,chan); //HS
      //res+=S0(m,k,chan)*GammGlobal(x,k,n,chan,M2)+std::conj(GammGlobal(x,k,m,chan,M2))*S0(k,n,chan); //Me
    }
                                                                                      
  return res*alpha_s/2./M_PI*log(muR*muR/M2)/2.; //HS Convention Gamma=alpha_s/2.*M_PI*Gamma^1
  //return -res*alpha_s/2./M_PI*log(muR*muR/M2)/2.; //Me Convention Gamma=alpha_s/2.*M_PI*Gamma^1
}


std::complex<double> SG(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0.;
  int Len=chan+2;
  std::complex<double> Nb=N(x[0])*exp(-Psi(1.));

  for(int k=0; k < Len; k++)
    {   
      res+=S0(m,k,chan)*std::conj(GammGlobal(x,n,k,chan,M2))+GammGlobal(x,m,k,chan,M2)*S0(k,n,chan); //HS
      //res+=S0(m,k,chan)*GammGlobal(x,k,n,chan,M2)+std::conj(GammGlobal(x,k,m,chan,M2))*S0(k,n,chan); // Me
      
    }
  
  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in SG !!" << std::endl;
    exit(1);
  }
  return res*alpha_s/2./M_PI*log(Nb); // With covention of literature log(1-2*lambda)*Gamma^1(without alpha)/2*beta0 factor 2 in addition by convention alpha_s/M_PI HS
  //return res*(-alpha_s/2./M_PI*log(Nb)); // With covention of literature log(1-2*lambda)*Gamma^1(without alpha)/2*beta0 factor 2 in addition by convention alpha_s/M_PI Me
}

std::complex<double> SGPhase(double *ppart, int m, int n, int chan, double M2, double x)
{
  std::complex<double> res=0.;
  int Len=chan+2;
  std::complex<double> Nb=N(x)*exp(-Psi(1.));

  for(int k=0; k < Len; k++)
    {   
      res+=S0(m,k,chan)*std::conj(GammPhase(ppart,n,k,chan,M2))+GammPhase(ppart,m,k,chan,M2)*S0(k,n,chan); //HS
      //res+=S0(m,k,chan)*GammGlobal(x,k,n,chan,M2)+std::conj(GammGlobal(x,k,m,chan,M2))*S0(k,n,chan); // Me
      
    }
  
  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in SG !!" << std::endl;
    exit(1);
  }
  return res*alpha_s/2./M_PI*log(Nb); // With covention of literature log(1-2*lambda)*Gamma^1(without alpha)/2*beta0 factor 2 in addition by convention alpha_s/M_PI HS
  //return res*(-alpha_s/2./M_PI*log(Nb)); // With covention of literature log(1-2*lambda)*Gamma^1(without alpha)/2*beta0 factor 2 in addition by convention alpha_s/M_PI Me
}

// std::complex<double> SuperG(double *x, int m, int n, int chan, double M2)
// {
//   std::complex<double> res=0., res2=0.;
//   int Len=chan+2;
//   std::complex<double> Nb=(N(x[0]))*exp(-Psi(1.)), lambdN=log(Nb)*alpha_s/2./M_PI*beta0;

//   for(int k=0; k < Len; k++)
//     {   
//       res+=(S0(m,k,chan)*std::conj(GammGlobal(x,n,k,chan,M2))+GammGlobal(x,m,k,chan,M2)*S0(k,n,chan)); //HS
//       //res+=S0(m,k,chan)*GammGlobal(x,k,n,chan,M2)+std::conj(GammGlobal(x,k,m,chan,M2))*S0(k,n,chan); // Me

//       for(int l=0; l< Len; l++)
// 	{
// 	  res2+=(S0(m,k,chan)*std::conj(GammGlobal(x,l,k,chan,M2))*std::conj(GammGlobal(x,n,l,chan,M2))+GammGlobal(x,m,k,chan,M2)*GammGlobal(x,k,l,chan,M2)*S0(l,n,chan))/2.; //HS
// 	  res2+=GammGlobal(x,m,k,chan,M2)*S0(k,l,chan)*std::conj(GammGlobal(x,n,l,chan,M2));
// 	}
      
//     }
//   if(isfinite(abs(res))==0){
//     std::cout << "!! Result non finite in SG !!" << std::endl;
//     exit(1);
//   }
//   return (res/beta0+res2/pow(beta0,2))*lambdN*lambdN; // With covention of literature log(1-2*lambda)*Gamma^1(without alpha)/2*beta0 factor 2 in addition by convention alpha_s/M_PI HS
// }


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

std::complex<double> S1Phase(double *ppart, int m, int n, int chan, double M2) // S1
{
  std::complex<double> res=0.;
  double t=0, u=0;
  std::complex<double> Lbet=0., i(0.0,1.0);
  double Xbet=0, bet=0.;
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet)));
  
  int sign=1;
  for(int l=0; l<4; l++)
    {
      if(l>0) sign=-1;
      t+=sign*pow(ppart[0+(npart+2)*l]-ppart[2+(npart+2)*l],2);
      u+=sign*pow(ppart[0+(npart+2)*l]-ppart[3+(npart+2)*l],2);
    }
  t-=pow(Mt,2); //t -> t1                                                                     
  u-=pow(Mt,2); //u -> u1                                                                      
  Xbet=1.+2./M2*t;

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


std::complex<double> S1oldMu(double *x, int m, int n, int chan, double M2) // S1
{
  std::complex<double> res=0.;
  double t=0, u=0;
  std::complex<double> Lbet=0., i(0.0,1.0), Nb=N(x[0])*exp(-Psi(1.));
  double Xbet=0, bet=0.;
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;
  //Xbet=Xcos(x[0])*bet;
  Lbet=(1.+bet*bet)/2./bet*(log((1.-bet)/(1.+bet)));
  t=-M2/2.*(1.-Xbet); //t1
  u=-M2/2.*(1.+Xbet); //u1
  
    if(chan==0)
    {
      exit(1);
    }

    else if (chan==1)
    {
      if(n==0 && m==0) res=M_PI*M_PI/3.*double(Nc)+2.*CF*log((1.+bet)/(1.-bet))/bet+CF*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)))+2.*CF*(1.+Lbet)*log(muR*muR*Nb*Nb/M2)+double(Nc)*pow(log(muR*muR*Nb*Nb/M2),2);

      else if(n==1 && m==0) res=4.*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)))-4.*log(muR*muR*Nb*Nb/M2)*log(t/u);

      else if(m==1 && n==1) res=2.*CF*log((1.+bet)/(1.-bet))/bet+double(Nc)*(M_PI*M_PI/3.-0.5*pow(log((1.+bet)/(1.-bet)),2))-double(Nc)*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet)))-(double(Nc)/2.-CF)*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)))+double(Nc)*pow(log(muR*muR*Nb*Nb/M2),2)+log(muR*muR*Nb*Nb/M2)*(2.*CF-double(Nc)*log(u*t/M2/pow(Mt,2))+Lbet*(2.*CF-double(Nc)));

      else if(m==1 && n==2) res=double(Nc)*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)))-double(Nc)*log(muR*muR*Nb*Nb/M2)*log(t/u);

      else if(m==1 && n==0) res=2.*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)))-2.*log(muR*muR*Nb*Nb/M2)*log(t/u);

      else if(m==2 && n==1) res=(double(Nc*Nc)-4.)/double(Nc)*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))-gsl_sf_dilog(1.+M2/2./t*(1.+bet))-gsl_sf_dilog(1.+M2/2./t*(1.-bet)))-(double(Nc*Nc)-4.)/double(Nc)*log(muR*muR*Nb*Nb/M2)*log(t/u);

      else if(m==2 && n==2) res=2.*CF*log((1.+bet)/(1.-bet))/bet+double(Nc)*(M_PI*M_PI/3.-0.5*pow(log((1.+bet)/(1.-bet)),2))-double(Nc)*(gsl_sf_dilog(1.+M2/2./u*(1.+bet))+gsl_sf_dilog(1.+M2/2./u*(1.-bet))+gsl_sf_dilog(1.+M2/2./t*(1.+bet))+gsl_sf_dilog(1.+M2/2./t*(1.-bet)))-(double(Nc)/2.-CF)*(1.+pow(bet,2))/bet*(gsl_sf_dilog(2.*bet/(bet-1.))-gsl_sf_dilog(2.*bet/(1.+bet)))+double(Nc)*pow(log(muR*muR*Nb*Nb/M2),2)+log(muR*muR*Nb*Nb/M2)*(2.*CF-double(Nc)*log(u*t/M2/pow(Mt,2))+Lbet*(2.*CF-double(Nc)));

    }

    // std::cout.precision(10);                                                                  
    // std::cout << " S1 =  " << res << " m = " << m << " n = " << n << " beta = " << bet << " M2 = " << M2 << " t = " << t << " muR*muR/M2 = " << muR*muR/M2 << std::endl;                  

    return res*alpha_s/2./M_PI;
}


std::complex<double> S1t(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0., res2=0., un(1.,0.), g1=0.,g2=0.;
  int Len=chan+2;
  //std::complex<double> Nb=N(x[0]);
  
  for(int k=0; k < Len; k++)
    {
      //res+=S0(m,k,chan)*S1old(x,k,n,chan,M2); //*Nc in old implementation
      //res+=S0(m,k,chan)*(S1old(x,k,n,chan,M2)+g3Nindep(x,k,n,chan,M2)+Gamma_sNindep(x,k,n,chan,M2)); // S1(mu_s)-> initial value, g3 is the N independant part of the colienar factor after integration
      //res+=S0(m,k,chan)*(S1old(x,k,n,chan,M2)+Gamma_sNindep(x,k,n,chan,M2)); // S1(mu_s)-> initial value
      //res2+=S1oldMu(x,m,k,chan,M2)*S0(k,n,chan); // S1(muR)
      res+=S1old(x,m,k,chan,M2)*S0(k,n,chan); // S1(muS)
    }

  // std::complex<double> Nb=N(x[0])*exp(-Psi(1.)), lambdN=alpha_s/2./M_PI*beta0*log(Nb);
  // g1=2.*Ci(chan)*lambdN/beta0;

  // g2=-2.*Ci(chan)/beta0*lambdN*log(M2/muR/muR);
  // g2+=2.*Ci(chan)*log(muF/muR)*2.*lambdN/beta0;

  // g1*=2.;
  // g2*=2.;

  // std::cout << " S1(mu_S) =  " << res << " S1(muR) = " << res2  << " S1 evolved = " << res  + Gamma_sNindep(x,m,n,chan,M2) + S0(m,n,chan)*(g3Nindep(chan,M2)+g1*log(Nb)+g2) + SG(x,m,n,chan,M2)  << " m = " << m << " n = " << n << std::endl;
  //std::cout << " S1(mu_S) =  " << res << " S1(muR) = " << res2  << " S1 evolved = " << res  + Gamma_sNindep(x,n,m,chan,M2) + S0(m,n,chan)*(g3Nindep(chan,M2)+g1*log(Nb)+g2) + SG(x,n,m,chan,M2)  << " m = " << m << " n = " << n << std::endl;  

  //if(m==2) exit(1);
  
  res+=Gamma_sNindep(x,m,n,chan,M2);
  res+=S0(m,n,chan)*g3Nindep(chan,M2);
  res+=S0(m,n,chan)*Beam(chan,M2);

  // res+=S0(m,n,chan)*alpha_s/2./M_PI*(-Psi(1.))*Ci(chan)*log(muF*muF/M2)*4.;
  // res+=S0(m,n,chan)*alpha_s/2./M_PI*pow(-Psi(1.),2.)*Ci(chan)*4.;
  // res+=SG(x,m,n,chan,M2)*(-Psi(1.))/log(Nb);
  
  return res;
}

std::complex<double> S1tPhase(double *ppart, int m, int n, int chan, double M2)
{
  std::complex<double> res=0., res2=0., un(1.,0.), g1=0.,g2=0.;
  int Len=chan+2;
  //std::complex<double> Nb=N(x[0]);
  
  for(int k=0; k < Len; k++)
    {
      res+=S1Phase(ppart,m,k,chan,M2)*S0(k,n,chan); // S1(muS)
    }

  
  res+=Gamma_sNindepPhase(ppart,m,n,chan,M2);
  res+=S0(m,n,chan)*g3Nindep(chan,M2);
  res+=S0(m,n,chan)*Beam(chan,M2);

  // res+=S0(m,n,chan)*alpha_s/2./M_PI*(-Psi(1.))*Ci(chan)*log(muF*muF/M2)*4.;
  // res+=S0(m,n,chan)*alpha_s/2./M_PI*pow(-Psi(1.),2.)*Ci(chan)*4.;
  // res+=SG(x,m,n,chan,M2)*(-Psi(1.))/log(Nb);
  
  return res;
}


std::complex<double> S1G(double *x, int m, int n, int chan, double M2)
{
  std::complex<double> res=0., res2=0.;
  int Len=chan+2;
  std::complex<double> Nb=N(x[0])*exp(-Psi(1.));
  
  
  for(int k=0; k < Len; k++)
    {    
      for(int l=0; l < Len; l++)
	{   
	  res+=S1old(x,m,k,chan,M2)*S0(k,l,chan)*std::conj(GammGlobal(x,n,l,chan,M2))+GammGlobal(x,m,k,chan,M2)*S1old(x,k,l,chan,M2)*S0(l,n,chan); //HS
	  res2+=(S0(m,k,chan)*std::conj(GammGlobal(x,l,k,chan,M2))*std::conj(GammGlobal(x,n,l,chan,M2))+2.*GammGlobal(x,m,k,chan,M2)*S0(k,l,chan)*std::conj(GammGlobal(x,n,l,chan,M2))+GammGlobal(x,m,k,chan,M2)*GammGlobal(x,k,l,chan,M2)*S0(l,n,chan))*alpha_s/2./M_PI*log(Nb)*alpha_s/2./M_PI*log(muR*muR/M2)/2.; //HS

	}
    }

  res2+=(g3Nindep(chan,M2)+Beam(chan,M2))*SG(x,m,n,chan,M2);
  
  return res*alpha_s/2./M_PI*log(Nb)+res2; // With covention of literature log(1-2*lambda)*Gamma^1(without alpha)/2*beta0 factor 2 in addition by convention alpha_s/M_PI HS
}
