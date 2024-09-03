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

extern "C"{
  void __nintlib_MOD_romberg_nd(double(*)(int&, double *, double&, int&, double&) ,double*,double*,int&,int*,int&,double&,double&,int&,int&,double&,int&,double&);
}


double NewTot(int& dim_num, double* x, double& sc, int& mel, double& M2){
  return TotDiff(x,sc,mel,M2);
}

double Integrate(double& res, double& acc, double& sc, int& mel, double& M2)
{
  //double init=1e-10, end=1.-init; //Not 0 to 1 for finite evaluations only
  double init=1e-10, end=1.; //Not 0 to 1 for finite evaluations only
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
