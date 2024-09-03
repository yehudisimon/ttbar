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
double Xsin(double);
std::complex<double> GammGlobal(double*, int, int, int , double );
double H0(double*,int,int,int,double);
double S0(int,int,int);
std::complex<double> S1(double*,int,int,int,double);
std::complex<double> S1t(double*,int,int,int,double);
std::complex<double> SG(double*,int,int,int,double);
std::complex<double> SSuperG(double*,int,int,int,double);
std::complex<double> GammGlobal(double*,int,int,int,double);
void pdfFit(double&, double (*)[8], double, double, double);
double Beam(int);
double g3Nindep(int,double);

extern "C"{
  void update_mur_(double&);
  void ml5_2_sloopmatrix_thres_(double[4][4],complex<double>(*),double&,double[3][3],int&);
}

// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi, xmax, xmin;
extern double A[8][8], A1min;
extern int Nc, Nf, chan;

std::complex<double> TraceBornDiff(double *x, int chan, double M2)
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;
  std::complex<double> res=0.;    
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


std::complex<double> TraceHSDiff(double *x, int chan, double M2, complex<double> *xx)
{
  // chan is channel index, 0 for qqbar and 1 for gluon gluon
  const int Len=chan+2;  
  std::complex<double> res=0.;
  std::complex<double> **R, **Rm, **Rdag, **Rdagm1, *Lamb;
  std::complex<double> Ht[Len][Len]={0}, i(0.0,1.0);
  std::complex<double> H1t[Len][Len]={0};

  std::complex<double> Nb, DeltaProc;
  std::complex<double> St[Len][Len]={0};
  std::complex<double> S0t[Len][Len]={0};

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

  
  double alpha_p=0.118; // Fixed value from MadLoop -> need rescaling for H
  
  for(int k=0; k < Len; k++)
    {
      for(int l=0; l < Len; l++)
	{
	  for(int mm=0; mm < Len; mm++)
	    {
	      for(int nn=0; nn < Len; nn++)
		{
		  Ht[k][l]+=Rdagm1[k][mm]*xx[0+nn*4+mm*Len*4]*pow(alpha_s/alpha_p,2)/2./M2*Rm[nn][l]; 
		  H1t[k][l]+=Rdagm1[k][mm]*xx[1+nn*4+mm*Len*4]*pow(alpha_s/alpha_p,3)/2./M2*Rm[nn][l];
		  St[k][l]+=R[k][mm]*(S0(mm,nn,chan)+S1t(x,mm,nn,chan,M2))*Rdag[nn][l]; 
		  S0t[k][l]+=R[k][mm]*S0(mm,nn,chan)*Rdag[nn][l]; 
		}
	    }
	  DeltaProc=exp(-(Lamb[k]+std::conj(Lamb[l]))/2./beta0*log(1.-alpha_s/M_PI*beta0*log(Nb)));
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
   
   return res;
 
}

// *************************************************** //
// Soft colinear radiations //
// *************************************************** //


std::complex<double> ColinearDiff(double *x, int chan, double M2) 
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
  g2+=-(2.*lambdN+log(1.-2.*lambdN))/pow(beta0,2)*Ci(chan)*(67.*double(Nc)/18.-double(Nc)/6.*pow(M_PI,2)-5.*double(Nf)/9.);
  g2+=beta1/pow(beta0,3)*Ci(chan)*(2.*lambdN+log(1.-2.*lambdN)+0.5*pow(log(1.-2.*lambdN),2));

  //Sudakov Colinear= Delta_i^2 for qq > ttbar or gg > ttbar
  g1*=2.;
  g2*=2.;

  G=std::exp(g1*log(Nb)+g2);

  // result
  return G;
}

std::complex<double> MellinPDFDiff(double *x, int chan, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform           
{
  std::complex<double> i(0.0,1.0), fAB=0, res=0;

  std::complex<double> q[5], g, qbar[5], Npdf;
  Npdf=N(x[0]);
  SetPDFN(Npdf,q,qbar,g,A); // PDFs with N
  
  // Symmetrized PDF depend on channel, 0 is qqbar and 1 is gg                                     
  if(chan==0)
    {
      for(int j=0; j < Nf; j++) fAB+=2.*q[j]*qbar[j];        
    }
  else if(chan==1) fAB=g*g;

  else std::cout << "Error in chan ID" << std::endl;

  //Global factors                                              
  res=fAB;
  return res;
}


std::complex<double> GlobalDiff(double *x, double sc, double M2) // Resummed formula, the numericall integration does the inversed Mellin transform
{
  
  std::complex<double> i(0.0,1.0), res=0, jac=0;
  double jac2=0, bet=0;
  
  //Tau computation
  double tau = M2/sc;
  
  // Change of variable N
  jac = (cos(Phi)+i*sin(Phi))/x[0];
  
  // Change of variable cos
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  jac2=20.-10.*pow(tanh(20.*x[1]-5.),2)-10.*pow(tanh(20.*x[1]-15.),2);

  //Global factors
  res=pow(tau,-N(x[0])+1.)*2./2./M_PI; //*2 from Im part, /2*M_PI from inverse Mellin transform, N-1 because the functions are supposed to be evaluated at N+1 when taking the Mellin inverse in N
  res*=bet/16./M_PI/M2; 
  
  // Result
  res=res*jac*jac2*0.38937966e9;

  return res;
}

//------------------------------------------------------------//
//--------------------Total Integrand-------------------------//
//------------------------------------------------------------//
double TotDiff(double *x, double& sc, int& mel, double& M2)
{
  double res=0., bet, t, hgg,  u, Xbet;
  std::complex<double> temp=0, i(0.0,1.0), Nb;
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  Xbet=Xcos(x[1])*bet;
  t=-M2/2.*(1-Xbet);

  if(mel==0) // Born
    {
      int Len=chan+2;
      std::complex<double> tmp=0.;
      tmp=TraceBornDiff(x,chan,M2);
      res+=std::imag(tmp*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));
    }

  else if(mel==1) //NLL
    {
      int Len=chan+2;

      complex<double> *xx;
      xx=new complex<double>[4*3*3];
      double ppart[4][4]={0.};
      int ret_code=0;
      double prec_ask=-1;
      double prec_f[3][3]={(0.,0.)};
      //kinematics
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
      
      ml5_2_sloopmatrix_thres_(ppart,xx,prec_ask,prec_f,ret_code);
      
      res+=std::imag(TraceHSDiff(x,chan,M2,xx)*ColinearDiff(x,chan,M2)*MellinPDFDiff(x,chan,M2)*GlobalDiff(x,sc,M2));
      
    }
  
  if(isfinite(res)==0){
    std::cout << "!! Result non finite !!" << std::endl;
    exit(1);
  }
  
  return res;
}
