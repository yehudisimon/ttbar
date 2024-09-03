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
// -------- Functions --------------------------------- //
std::complex<double> H1(double*,int,int,int,double);
void EigenvGlobal(double*,double,int,std::complex<double>**,std::complex<double>**,std::complex<double>**,std::complex<double>**,std::complex<double>*);
double Sqrt(double);
std::complex<double> Psi(std::complex<double>);
std::complex<double> Gamma(const std::complex<double>);
std::complex<double> B(std::complex<double>&,double,double);
void SetPDFN(std::complex<double>&,std::complex<double>*,std::complex<double>*,std::complex<double>&,double (*)[8]);
void EvolvePDF(std::complex<double>&, std::complex<double>*, std::complex<double>*, std::complex<double>&, double);
double Ci(int);
double B1(int);
std::complex<double> N(double);
double Xcos(double);
std::complex<double> GammGlobal(double*, int, int, int , double );
double H0(double*,int,int,int,double);
double S0(int,int,int);
std::complex<double> S1(double*,int,int,int,double);
std::complex<double> S1t(double*,int,int,int,double);
std::complex<double> SG(double*,int,int,int,double);
void DisplayS1(double*,int,double);
void DisplayG1(double*,int,double);

// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi;
extern double A[8][8];
extern int Nc, Nf;

std::complex<double> TraceBorn(double *x, int chan, double M2)
{
    // rep is representation index (0 for singlet, 1 for octet (asymmetric if relevant), 2 for symmetric octet
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  // Hard part
  
  std::complex<double> res=1.;
    
  res=0.; // Compute the trace of Hard and Soft
  for(int ka=0; ka < Len; ka++)
     {
       for(int ji=0; ji < Len; ji++)
	 {
	   res+=S0(ka,ji,chan)*H0(x,ji,ka,chan,M2);
	 }
     }
   return res;
}

std::complex<double> TraceHSGlobal(double *x, int chan, double M2)
{
  // rep is representation index (0 for singlet, 1 for octet (asymmetric if relevant), 2 for symmetric octet
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  // Hard part
  
  std::complex<double> res=1.;
  std::complex<double> **R, **Rm, **Rdag, **Rdagm1, *Lamb;
  std::complex<double> Ht[Len][Len]={0}, i(0.0,1.0);

  std::complex<double> Nb, DeltaProc;
  std::complex<double> St[Len][Len]={0};

  R=new std::complex<double> *[Len];
  Rm=new std::complex<double> *[Len];
  Rdag=new std::complex<double> *[Len];
  Rdagm1=new std::complex<double> *[Len];
  Lamb=new std::complex<double> [Len];
  
  for(int h=0; h < Len; h++)
    {
      R[h]=new std::complex<double>[Len];
      Rdag[h]=new std::complex<double>[Len];
      Rm[h]=new std::complex<double>[Len];
      Rdagm1[h]=new std::complex<double>[Len];
    }

  
  EigenvGlobal(x,M2,chan,Rdag,Rdagm1,R,Rm,Lamb); // Computes rotation matrices and eigenvalues from soft anomalous dimension matrix
  Nb=N(x[0])*exp(-Psi(1.));
   
  for(int k=0; k < Len; k++)
    {
      for(int l=0; l < Len; l++)
	{
	  for(int mm=0; mm < Len; mm++)
	    {
	      for(int nn=0; nn < Len; nn++)
		{
		  Ht[k][l]+=R[k][mm]*(H0(x,mm,nn,chan,M2)+H1(x,mm,nn,chan,M2))*Rdag[nn][l];
		  St[k][l]+=Rdagm1[k][mm]*(S0(mm,nn,chan)+S1t(x,mm,nn,chan,M2))*Rm[nn][l];
		  //std::cout << " Chan = " << chan << " M2 = " << M2 << " muR = " << muR << " muF = " << muF << " Xbet = " << Xcos(x[1])*pow(1.-4.*Mt*Mt/M2,0.5) << " S1 = " << S1t(x,mm,nn,chan,M2)*2.*M_PI/alpha_s << " m = " << mm << " n = " << nn << std::endl;

		  //St[k][l]+=Rdagm1[k][mm]*(S0(mm,nn,chan)+S1(x,mm,nn,chan,M2))*Rm[nn][l];
		}
	    }
	  DeltaProc=exp((Lamb[l]+std::conj(Lamb[k]))/beta0*log(1.-alpha_s/M_PI*beta0*log(Nb)));
	  St[k][l]*=DeltaProc;
	}
    }
  
  res=0.; // Compute the trace of Hard and Soft
   for(int ka=0; ka < Len; ka++)
     {
       for(int ji=0; ji < Len; ji++) res+=St[ka][ji]*Ht[ji][ka];
     }

   for (int i=0; i<Len; i++)
     {
       delete [] R[i];
       delete [] Rm[i];
       delete [] Rdagm1[i];
       delete [] Rdag[i];
     }
    
   delete [] R;
   delete [] Rm;
   delete [] Rdag;
   delete [] Rdagm1;
   delete [] Lamb;

   return res;
 
}

std::complex<double> TraceHSExpGlobal(double *x, int chan, double M2)
{
  // rep is representation index (0 for singlet, 1 for octet (asymmetric if relevant), 2 for symmetric octet
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  // Hard part
  
  std::complex<double> res=1.;

  res=0.; // Compute the trace of Hard and Soft
  for(int ka=0; ka < Len; ka++)
     {
       for(int ji=0; ji < Len; ji++)
	 {
	   res+=(S0(ka,ji,chan)+S1t(x,ka,ji,chan,M2)+SG(x,ka,ji,chan,M2))*H0(x,ji,ka,chan,M2);
	   //std::cout << " S1t = " << S1t(x,ka,ji,chan,M2) << " for M2 = " << M2 << " chan = " << chan << " ka  = " << ka << " ji = " << ji << std::endl; 
	   //res+=(S0(ka,ji,chan)+S1(x,ka,ji,chan,M2)+SG(x,ka,ji,chan,M2))*H0(x,ji,ka,chan,M2);
	   res+=S0(ka,ji,chan)*H1(x,ji,ka,chan,M2); //Virtual contribution
	 }
     }
  //exit(2);
  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in Tr expanded !!" << std::endl;
    exit(1);
  }
  
   return res;
}


// *************************************************** //
// Soft colinear radiations //
// *************************************************** //


std::complex<double> ColinearGlobal(double *x, int chan, double M2) 
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  // Computation of the kinematical variables
  std::complex<double> i(0.0,1.0), Nb=0, G=0, g1=0, g2=0;
  double Q=pow(M2,0.5);
  //Change of variable
  Nb=N(x[0])*exp(-Psi(1.));
  std::complex<double> lambdN=alpha_s/2./M_PI*beta0*log(Nb);
  // g1, g2 and corresponding G functions of the N variable, independant of the correction considered (global) cf. Nucl. Phys. B. 529 (1998) Bociani, Catani, Mangano, Nason or 1005.2909 Debove, Fuks, Klasen or 0805.1885 

  g1=Ci(chan)/beta0/lambdN*(2.*lambdN+(1.-2.*lambdN)*log(1.-2.*lambdN));
  
  g2=Ci(chan)/beta0*log(1.-2.*lambdN)*2.*log(Q/muR);
  g2+=2.*Ci(chan)*log(muF/muR)*2.*lambdN/beta0;
  g2+=-(2.*lambdN+log(1.-2.*lambdN))/pow(beta0,2)*Ci(chan)*(67.*Nc/18.-Nc/6.*pow(M_PI,2)-5.*Nf/9.);
  g2+=beta1/pow(beta0,3)*Ci(chan)*(2.*lambdN+log(1.-2.*lambdN)+0.5*pow(log(1.-2.*lambdN),2));

  //Sudakov Colinear= Delta_i^2 for qq > ttbar or gg > ttbar
  g1*=2.;
  g2*=2.;
    
  G=std::exp(g1*log(Nb)+g2);
  // result
  return G;
}

std::complex<double> ColinearExpGlobal(double *x, int chan, double M2) 
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  // Computation of the kinematical variables
  std::complex<double> i(0.0,1.0), Nb=0, G=0, g1=0, g2=0;
  double Q=pow(M2,0.5);
  
  // Change of variable
  Nb=N(x[0])*exp(-Psi(1.));
  
  std::complex<double> lambdN=alpha_s/2./M_PI*beta0*log(Nb);
  g1=2.*Ci(chan)*lambdN/beta0;
  
  g2=-2.*Ci(chan)/beta0*lambdN*2.*log(Q/muR);
  g2+=2.*Ci(chan)*log(muF/muR)*2.*lambdN/beta0;
  
  g1*=2.;
  g2*=2.;
  
  
  G=1.+g1*log(Nb)+g2;
  
  if(isfinite(abs(G))==0){
    std::cout << "!! Result non finite in Colinear expanded !!" << std::endl;
    exit(1);
  }
  // result
  return G-1.;
}

std::complex<double> MellinPDFGlobal(double *x, int chan, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform           
{
  std::complex<double> i(0.0,1.0), fAB=0, res=0;

  std::complex<double> q[5], g, qbar[5], Npdf;
  Npdf=N(x[0])+1.;
  SetPDFN(Npdf,q,qbar,g,A); // PDFs with N

  // Symmetrized PDF depend on channel, 0 is qqbar and 1 is gg                                     
  if(chan==0)
    {
      for(int j=0; j < Nf; j++) fAB+=2.*q[j]*qbar[j];        
    }
  else if(chan==1) fAB=g*g;//*2.; symmetrized because same particle        

  else std::cout << "Error in chan ID" << std::endl;

  //Global factors                                              
  res=fAB;
  return res;
}


std::complex<double> GlobalGlobal(double *x, double sc, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform
{
  
  std::complex<double> i(0.0,1.0), res=0, jac=0;
  double jac2=0, bet=0;
  
  //Tau computation
  double tau = M2/sc;
  
  // Change of variable
  jac = (cos(Phi)+i*sin(Phi))/pow(1.-x[0],2);
  
  // Change of variable 2
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  jac2=20.-10.*pow(tanh(20.*x[1]-5.),2)-10.*pow(tanh(20.*x[1]-15.),2);

  //Global factors
  res=pow(tau,-N(x[0]))*2./2./M_PI; //*2 from Im part, /2*M_PI from inverse Mellin transform, N-1 because the functions are supposed to be evaluated at N+1 when taking the Mellin inverse in N
  res*=bet/16./M_PI/M2; 
  
  // Result
  res=res*jac*jac2*0.38937966e9;//*jac3;

  return res;
}

//------------------------------------------------------------//
//--------------------Total Integrand-------------------------//
//------------------------------------------------------------//


double TotGlobal(double *x, double& sc, int& mel, double& M2)
{
  double res=0., bet, t, hgg,  u, Xbet;
  std::complex<double> temp=0, i(0.0,1.0), Nb;
  
  Nb=N(x[0])*exp(-Psi(1.));
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;
  t=-M2/2.*(1-Xbet);
  u=-M2/2.*(1+Xbet);

  //hgg=pow(M2,2)/2./t/u*(pow(t,2)/pow(M2,2)+pow(u,2)/pow(M2,2)+4.*pow(Mt,2)/M2-4.*pow(Mt,4)/t/u);

  if(mel==2)
    {
      for(int chan=1; chan < 2; chan++)
	{
	  res+=std::imag(TraceHSGlobal(x,chan,M2)*ColinearGlobal(x,chan,M2)*MellinPDFGlobal(x,chan,M2)*GlobalGlobal(x,sc,M2));
	}
    }

  else if(mel==3)
    {
      for(int chan=1; chan < 2; chan++)
	{
	  //if(chan==0) temp=CF/2.*pow(alpha_s,2)*3./8./(double)Nc*(1.+pow(Xbet,2)+4.*Mt*Mt/M2);
	  //else if(chan==1) temp=pow(alpha_s,2)*3./8./(double(Nc*Nc)-1.)*hgg*(1./double(Nc)+double(Nc)/2.*pow(t-u,2)/pow(M2,2)+(-4.+pow(Nc,2))/2./double(Nc));
	  res+=std::imag((TraceHSExpGlobal(x,chan,M2)+ColinearExpGlobal(x,chan,M2)*TraceBorn(x,chan,M2))*MellinPDFGlobal(x,chan,M2)*GlobalGlobal(x,sc,M2));
	}
    }

  else if(mel==0)
    {
      for(int chan=0; chan < 2; chan++)
	{
	 int Len=chan+2;
	 std::complex<double> tmp=0.;
	 tmp=TraceBorn(x,chan,M2);
	 // for(int ka=0; ka < Len; ka++)
	 //    {
	 //      for(int ji=0; ji < Len; ji++)
	 // 	{
	 // 	 tmp+=S0(ka,ji,chan)*H0(x,ji,ka,chan,M2);
	 // 	}
	 //    }
	 res+=std::imag(tmp*MellinPDFGlobal(x,chan,M2)*GlobalGlobal(x,sc,M2));
	}
    }
  
  if(isfinite(res)==0){
    std::cout << "!! Result non finite !!" << std::endl;
    exit(1);
  }
  return res;
}
