// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <sstream>    	// String streams               //
#include <string>	// Strings                      //
#include <iostream>
#include <fstream>
#include <complex>
#include <stdio.h>
#include <math.h>
// -------- Classes ----------------------------------- //
#include "process.h"      	// Process definition   //
#include <gsl/gsl_monte_vegas.h>// GSL                  //
#include "LHAPDF/LHAPDF.h"
// -------- Functions --------------------------------- //

using namespace LHAPDF;
double TotDiff(double*,double&,int&,double&);
double TotPhaseSpace(double*, size_t, void*);
double TotVegas(double*, size_t, void*);
void DisplayXsec(const double&,const double&,const std::string&); //

extern "C"{
  void __nintlib_MOD_romberg_nd(double(*)(int&, double *, double&, int&, double&) ,double*,double*,int&,int*,int&,double&,double&,int&,int&,double&,int&,double&);
}

extern int npart;

double NewTot(int& dim_num, double* x, double& sc, int& mel, double& M2){
  return TotDiff(x,sc,mel,M2);
}

double Integrate(double& res, double& acc, double& sc, int& mel, double& M2)
{
  double init=1e-5, end=1.; //Not 0 to 1 for finite evaluations only
  double resbis=0.;
  int dim_num=2, sum_num[dim_num]={30,30};
  int it_max=1, ind=0, eval_num;
  double a[dim_num], b[dim_num], abis[dim_num], bbis[dim_num];
  
  for(int i=0; i < dim_num; i++)
    {
      a[i]=init;
      abis[i]=init/10.;
      b[i]=end;
      bbis[i]=init*9./10.+end;
    }
  
  acc=1e-5;

  __nintlib_MOD_romberg_nd(NewTot,a,b,dim_num,sum_num,it_max,acc,res,ind,eval_num,sc,mel,M2);
  //__nintlib_MOD_romberg_nd(NewTot,abis,bbis,dim_num,sum_num,it_max,acc,resbis,ind,eval_num,sc,mel,M2);
  //if(abs(res-resbis)>abs(res)*5.*acc) std::cout << " WARNING: INTEGRATION SENSITIVE TO INTEGRAL RANGES" << std::endl;
  
  return 0.;
}


double IntegrateVegas(double& res, double& err, double& chi, double& sc, int& mel, double& M2)
{
  // Initializing the random number generator                                                     
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  time_t seed = time(NULL);
  gsl_rng_set(r,seed);
  
  // Integrand                                                                               
  gsl_monte_function I; I.dim=0;
  size_t calls=0;
  
  // Number of integrations and calls                                                       
  calls=2000;
  I.dim=2;
  
  struct my_f_params { double sc; int mel; double M2; };
  struct my_f_params params = { sc, mel, M2 };
  
  
  // Selecting the good integrand + includes in the variable factor the                       
  // constant pieces of the squared matrix element
  
  I.f= &TotVegas;
  I.params=&params;
  // Integration bounds                                                                        
  double xmin[I.dim], xmax[I.dim];
  for(size_t i=0; i<I.dim; i++) { xmin[i]=1e-7, xmax[i]=1.; } // should be called once            
  
  // Integration                                                                                
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(I.dim);
  
  // Warm-up                                                                                   
  s->stage=0;  s->iterations=5;
  s->verbose=0;
  gsl_monte_vegas_integrate(&I, xmin, xmax, I.dim, calls/50, r, s, &res, &err);
  chi = s->chisq;
  double prec=100.;
  
  // Real things: stops if the precision reaches 1% or if one oscillates                      
  int counter=1; s->iterations=3;
  while(prec>5e-3 && counter<=10)
    {
      std::ostringstream ocnt; ocnt << counter; std::string cntstr= ocnt.str();
      s->stage=1;
      gsl_monte_vegas_integrate(&I, xmin, xmax, I.dim, calls, r, s, &res, &err);
      prec=std::abs(err/res);
      DisplayXsec(res,err,"Refine-"+cntstr);
      counter++;
      calls*=5;
    }
  chi = s->chisq;
  DisplayXsec(res,err,"final");
  // Cleaning the memory and closing the file                                                   
  gsl_monte_vegas_free(s);
  gsl_rng_free(r);

  return 0;
}

double IntegratePhase(double& res, double& err, double& chi, double& sc, int& mel, double& M2)
{
  // Initializing the random number generator                                                     
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  time_t seed = time(NULL);
  gsl_rng_set(r,seed);
  
  // Integrand                                                                               
  gsl_monte_function I; I.dim=0;
  size_t calls=0;
  
  // Number of integrations and calls                                                       
  calls=10000;
  I.dim=3*(npart-1)+1;
  
  struct my_f_params { double sc; int mel; double M2; };
  struct my_f_params params = { sc, mel, M2 };
  
  
  // Selecting the good integrand + includes in the variable factor the                       
  // constant pieces of the squared matrix element
  
  I.f= &TotPhaseSpace;
  I.params=&params;
  // Integration bounds                                                                        
  double xmin[I.dim], xmax[I.dim];
  for(size_t i=0; i<I.dim; i++) { xmin[i]=1e-7, xmax[i]=1.; } // should be called once            
  
  // Integration                                                                                
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(I.dim);
  
  // Warm-up                                                                                   
  s->stage=0;  s->iterations=5;
  s->verbose=0;
  gsl_monte_vegas_integrate(&I, xmin, xmax, I.dim, calls/500, r, s, &res, &err);
  chi = s->chisq;
  double prec=100.;
  
  // Real things: stops if the precision reaches 1% or if one oscillates                      
  int counter=1; s->iterations=3;
  while(prec>5e-2 && counter<=10)
    {
      std::ostringstream ocnt; ocnt << counter; std::string cntstr= ocnt.str();
      s->stage=1;
      gsl_monte_vegas_integrate(&I, xmin, xmax, I.dim, calls, r, s, &res, &err);
      prec=std::abs(err/res);
      DisplayXsec(res,err,"Refine-"+cntstr);
      counter++;
      calls*=5;
    }
  chi = s->chisq;
  DisplayXsec(res,err,"final");
  // Cleaning the memory and closing the file                                                   
  gsl_monte_vegas_free(s);
  gsl_rng_free(r);

  return 0;
}
