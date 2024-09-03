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
void EigenvPhase(double*,double,int,std::complex<double>**,std::complex<double>**,std::complex<double>**,std::complex<double>**,std::complex<double>*);
double Sqrt(double);
std::complex<double> Psi(std::complex<double>);
std::complex<double> Gamma(const std::complex<double>);
std::complex<double> B(std::complex<double>&,double,double);
void SetPDFN(std::complex<double>&,std::complex<double>*,std::complex<double>*,std::complex<double>&,double (*)[8]);
double Ci(int);
double B1(int);
std::complex<double> N(double);
std::complex<double> GammPhase(double*, int, int, int, double);
double H0Phase(double*,int,int,int,double);
double S0(int,int,int);
//std::complex<double> S1Phase(double*,int,int,int,double);
std::complex<double> S1tPhase(double*,int,int,int,double);
std::complex<double> SGPhase(double*,int,int,int,double,double);
void pdfFit(double&, double (*)[8], double, double, double);

// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi, xmax, xmin;
extern double A[8][8], A1min;
extern int Nc, Nf, chan;

std::complex<double> TraceBornPhase(double* ppart, int chan, double M2)
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  std::complex<double> res=0.;    
  res=0.; // Compute the trace of Hard and Soft
  for(int ka=0; ka < Len; ka++)
     {
       for(int ji=0; ji < Len; ji++)
	 {
	   res+=S0(ka,ji,chan)*H0Phase(ppart,ji,ka,chan,M2);
	 }
     }

  return res;
}

std::complex<double> TracePhase(double *ppart, int chan, double M2, complex<double> *xx, double xN)
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  std::complex<double> res=1.;
  std::complex<double> **R, **Rm, **Rdag, **Rdagm1, *Lamb;
  std::complex<double> Ht[Len][Len]={0}, i(0.0,1.0);

  std::complex<double> Nb, DeltaProc;
  std::complex<double> St[Len][Len]={0};
  std::complex<double> S0t[Len][Len]={0};
  std::complex<double> H1t[Len][Len]={0};

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
  
  EigenvPhase(ppart,M2,chan,Rdag,Rdagm1,R,Rm,Lamb); // Computes rotation matrices and eigenvalues from soft anomalous dimension matrix
  Nb=N(xN)*exp(-Psi(1.));
  
  double alpha_p=0.118;
  
  for(int k=0; k < Len; k++)
    {
      for(int l=0; l < Len; l++)
	{
	  for(int mm=0; mm < Len; mm++)
	    {
	      for(int nn=0; nn < Len; nn++)
		{
		  H1t[k][l]+=Rdagm1[k][mm]*xx[1+nn*4+mm*Len*4]*pow(alpha_s/alpha_p,3)/2./M2*Rm[nn][l]; //HS                                                                        
                  Ht[k][l]+=Rdagm1[k][mm]*H0Phase(ppart,mm,nn,chan,M2)*Rm[nn][l]; // HS                
 
                  St[k][l]+=R[k][mm]*(S0(mm,nn,chan)+S1tPhase(ppart,mm,nn,chan,M2))*Rdag[nn][l]; // HS
                  S0t[k][l]+=R[k][mm]*S0(mm,nn,chan)*Rdag[nn][l]; //HS              
		  
		  //H1t[k][l]+=R[k][mm]*xx[1+nn*4+mm*Len*4]*pow(alpha_s/alpha_p,3)/2./M2*Rdag[nn][l]; //Me
		  //Ht[k][l]+=R[k][mm]*H0Phase(ppart,mm,nn,chan,M2)*Rdag[nn][l]; //Me
		  
		  //S0t[k][l]+=Rdagm1[k][mm]*S0(mm,nn,chan)*Rm[nn][l]; //Me
		  //St[k][l]+=Rdagm1[k][mm]*(S0(mm,nn,chan)+S1Phase(ppart,mm,nn,chan,M2))*Rm[nn][l]; //Me

		}
	    }
	  DeltaProc=exp(-(Lamb[k]+std::conj(Lamb[l]))/2./beta0*log(1.-alpha_s/M_PI*beta0*log(Nb))); //HS convention
	  //DeltaProc=exp((Lamb[l]+std::conj(Lamb[k]))/2./beta0*log(1.-alpha_s/M_PI*beta0*log(Nb))); //Me convention

	  St[k][l]*=DeltaProc;
	  S0t[k][l]*=DeltaProc;
	}
    }
  
  res=0.; // Compute the trace of Hard and Soft
   for(int ka=0; ka < Len; ka++)
     {
       for(int ji=0; ji < Len; ji++)
	 {
	   res+=St[ka][ji]*Ht[ji][ka];
	   res+=S0t[ka][ji]*H1t[ji][ka];
	 }
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
   
  
   if(isfinite(abs(res))==0){
     std::cout << "!! Result non finite in Res !!" << std::endl;
     exit(1);
   }
   
   return res;
 
}

std::complex<double> TraceExpPhase(double *ppart, int chan, double M2, complex<double> *xx, double xN)
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  double alpha_p=0.118;
  
  std::complex<double> res=1.;

  res=0.; // Compute the trace of Hard and Soft
  for(int ka=0; ka < Len; ka++)
     {
       for(int ji=0; ji < Len; ji++)
	 {
	   res+=(S1tPhase(ppart,ka,ji,chan,M2)+SGPhase(ppart,ka,ji,chan,M2,xN))*H0Phase(ppart,ji,ka,chan,M2);
	   //res+=(SGPhase(ppart,ka,ji,chan,M2,xN))*H0Phase(ppart,ji,ka,chan,M2);
	   //res+=(S1tPhase(ppart,ka,ji,chan,M2))*H0Phase(ppart,ji,ka,chan,M2);
	   res+=S0(ka,ji,chan)*xx[1+ka*4+ji*Len*4]*pow(alpha_s/alpha_p,3)/2./M2; //Virtual contribution
	 }
     }

  if(isfinite(abs(res))==0){
    std::cout << "!! Result non finite in Tr expanded !!" << std::endl;
    exit(1);
  }
  
   return res;
}

// *************************************************** //
// Soft colinear radiations //
// *************************************************** //

std::complex<double> ColinearPhase(double x, int chan, double M2) 
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  std::complex<double> i(0.0,1.0), Nb=0, G=0, g1=0, g2=0;
  double Q=pow(M2,0.5);
  //Change of variable
  Nb=N(x)*exp(-Psi(1.));

  std::complex<double> lambdN=alpha_s/2./M_PI*beta0*log(Nb);
  // g1, g2 and corresponding G functions of the N variable, independant of the correction considered (global) cf. Nucl. Phys. B. 529 (1998) Bociani, Catani, Mangano, Nason or 1005.2909 Debove, Fuks, Klasen or 0805.1885 

  g1=Ci(chan)/beta0/lambdN*(2.*lambdN+(1.-2.*lambdN)*log(1.-2.*lambdN));
  
  g2=Ci(chan)/beta0*log(1.-2.*lambdN)*2.*log(Q/muR);
  g2+=2.*Ci(chan)*log(muF/muR)*2.*lambdN/beta0;
  g2+=-(2.*lambdN+log(1.-2.*lambdN))/pow(beta0,2)*Ci(chan)*(67.*double(Nc)/18.-double(Nc)/6.*pow(M_PI,2)-5.*double(Nf)/9.);
  g2+=beta1/pow(beta0,3)*Ci(chan)*(2.*lambdN+log(1.-2.*lambdN)+0.5*pow(log(1.-2.*lambdN),2));

  //g2+=g1*(-Psi(1.)); // N resummation
  //g2+=Psi(1.)/beta0*Ci(chan)/lambdN*(2.*lambdN+log(1.-2.*lambdN));

  //Sudakov Colinear= Delta_i^2 for qq > ttbar or gg > ttbar
  g1*=2.;
  g2*=2.;

  G=std::exp(g1*log(Nb)+g2);
  
  // result
  return G;
}

std::complex<double> ColinearExpPhase(double x, int chan, double M2) 
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  std::complex<double> i(0.0,1.0), Nb=0, G=0, g1=0, g2=0;
  double Q=pow(M2,0.5);
  
  // Change of variable
  Nb=(N(x))*exp(-Psi(1.));
  
  std::complex<double> lambdN=alpha_s/2./M_PI*beta0*log(Nb);
  g1=2.*Ci(chan)*lambdN/beta0;
  
  g2=-2.*Ci(chan)/beta0*lambdN*2.*log(Q/muR);
  g2+=2.*Ci(chan)*log(muF/muR)*2.*lambdN/beta0;

  //g2+=g1*(-Psi(1.)); // N resummation
  //g2+=(-Psi(1.))*lambdN*2.*Ci(chan)/beta0;
  
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

std::complex<double> MellinPDFPhase(double x, int chan, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform           
{
  std::complex<double> i(0.0,1.0), fAB=0, res=0;

  std::complex<double> q[5], g, qbar[5], Npdf;
  Npdf=N(x);
  SetPDFN(Npdf,q,qbar,g,A); // PDFs with N

  //g=1./Npdf/(Npdf+1);
  // Symmetrized PDF depend on channel, 0 is qqbar and 1 is gg                                     
  if(chan==0)
    {
      exit(1);
    }
  else if(chan==1) fAB=g*g;       

  else std::cout << "Error in chan ID" << std::endl;

  //Global factors                                              
  res=fAB;
  return res;
}


std::complex<double> GlobalPhase(double x, double sc, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform
{  
  std::complex<double> i(0.0,1.0), res=0, jac=0;
  double jac2=0, bet=0;
  double tau = M2/sc;
  
  // Change of variable
  jac = (cos(Phi)+i*sin(Phi))/x;
  
  //Global factors
  res=pow(tau,-N(x)+1.)*2./2./M_PI; //*2 from Im part, /2*M_PI from inverse Mellin transform, N-1 because the functions are supposed to be evaluated at N+1 when taking the Mellin inverse in N
  res*=1./M2;
  
  // Result
  res=res*jac*0.38937966e9;//*jac3;

  return res;
}
