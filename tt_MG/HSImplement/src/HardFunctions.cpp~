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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
// -------- Classes ----------------------------------- //
#include "process.h"	// Process//
// -------- Functions --------------------------------- //
// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi;
extern int Nc;

extern "C"{
  void gg2ttx_shfmatrix_(int&, int&, double&, double&, double&, double&,double(*)[2]);
}

double Xcos(double);

//std::complex<double> Sqrt(double x)
double Sqrt(double x)
{
  return pow(x,0.5);
}

std::complex<double> C0(double a,double b, double c, double m1, double m2, double m3) // C0 =-i*M_PI^2/2./Mt^2* lim_{eps->0} 3F2({1,1,1},{3/2,2},(M2+i*eps)/(4.*Mt^2))=-i*M_PI^2/2./Mt^2*lim_{eps->0} arcsin(sqrt(M2+i*eps)/2./Mt)^2/(M2/4./Mt^2) !!! Careful z=M2/4./Mt^2 > 1 here so no log formula to express arcin -> use gsl_complex function to evaluate. Careful of branch cut, here M2 as an infinitesimal imaginary part which changes the sign of imaginary part. See 0712.1851 for definition of the integral: C0({x})=i*M_PI^2*I_3({x})
{
  gsl_complex resgsl=gsl_complex_pow_real(gsl_complex_arcsin_real(Sqrt(c/4./m1)),2.);
  std::complex<double> i(0,1);
  resgsl=gsl_complex_mul_real(resgsl,4.*m1/c);
  std::complex<double> res=GSL_REAL(resgsl)-i*GSL_IMAG(resgsl);
  res*=-i*pow(M_PI,2.)/2./pow(Mt,2.);
  return res;
}


// *********************************************************************** //
//  Hard coefficients //
// *********************************************************************** //


double H0(double *x, int m, int n, int chan, double M2)
{
  double res=0.;
  double t=0, u=0, hgg=0, dR=0, Nsym=1.;
  std::complex<double> i(0.0,1.0);
  double Xbet=0, bet=0.;
  int tmp=0;

  //H0 symmetric
  if(m>n)
    {
      tmp=m;
      m=n;
      n=tmp;
    }
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);

  Xbet=Xcos(x[1])*bet;
  t=-M2/2.*(1-Xbet);
  u=-M2/2.*(1+Xbet);
 
  
  if(chan==0)
    {
      dR=double(Nc)*2.; //
    }
  else if(chan==1)
    {
      dR=2.*(pow(Nc,2)-1.);// CDR convention : (1-eps) ommited, to check
      Nsym=1.;
    }


  if(chan==0)
    {
      if(m*n==1)
	{
	  res=2.*((pow(t,2)+pow(u,2))/pow(M2,2)+2.*Mt*Mt/M2); //see 0805.1885 or SCET 1003.5827
	}
    }

  else if(chan==1)
    {
      hgg=pow(M2,2)/2./t/u*(pow(t,2)/pow(M2,2)+pow(u,2)/pow(M2,2)+4.*pow(Mt,2)/M2-4.*pow(Mt,4)/t/u);
      res=hgg;
      if(m==0)
	{
	  res*=1./double(Nc);
	  if(n==0) res*=1./double(Nc);
	  else if(n==1) res*=(t-u)/M2;
	}

      else if(m==1)
	{
	  res*=(t-u)/M2;
	  if(n==1) res*=(t-u)/M2;;
	}
    }
  
  
  return res/2./M2/dR/dR/Nsym*pow(alpha_s*4.*M_PI,2)*4.;  //Global factor with respect to SCET because of conventions M_ren -> M_ren^(0) and H^(0)_SCET = H^(0)/4, no 3/8 nor 8/3 here
}


std::complex<double> H1(double *x, int m, int n, int chan, double M2)
{
  const int Len=chan+2;

  //std::complex<double> H[Len][Len];
  std::complex<double> i(0.0,1.0), un(1.0,0.0), resi, res;
  double Xbet=0, bet=0, t=0,dR=0., Nsym=1.;
  int tmp;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  
  Xbet=Xcos(x[1])*bet;
  t=-M2/2.*(1.-Xbet);
  
  if(chan==0) dR=double(Nc)*2.;
  else if(chan==1) dR=(pow(Nc,2)-1.)*2.;

  
  if(chan==0) //H1qqb is symmetric (and hermitian) -> res_ij=2*Re(<C_i|A_0><A_1|C_j>)=res_ji
    {
      if(m>n)
	{
	  tmp=m;
	  m=n;
	  n=tmp;
	}
      if(n==0)
	{
	  resi=0; //n=0 -> m=0 -> res=0 because of symetrization 
	}

      else if (n==1)
	{
	  if(m==0)
	    {
	      resi=((-0.14814814814814814*(2.*pow(Mt,4) + M2*pow(Mt,2)*(-2. - 4.*Sqrt(1 \
- (4.*pow(Mt,2))/M2)) + pow(M2,2)*(0.25 + 1.*Sqrt(1 - \
(4.*pow(Mt,2))/M2)))*(M2 + 2.*t))/((1.*M2 - 4.*pow(Mt,2))*Sqrt(1 - \
(4.*pow(Mt,2))/M2)) - (0.1111111111111111*(pow(M2,2) - \
8.*M2*pow(Mt,2) + 8.*pow(Mt,4))*(M2 + 2.*t)*pow(std::log(un*(1. - \
1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - \
(4.*pow(Mt,2))/M2))),2))/(pow(M_PI,2)*(1.*M2 - 4.*pow(Mt,2))*Sqrt(1 - \
(4.*pow(Mt,2))/M2)) - (0.2222222222222222*pow(M2,2)*(1.*M2 + \
2.*t)*std::log(un*(-1.*M2)/pow(muR,2)))/(pow(M_PI,2)*(1.*M2 - \
4.*pow(Mt,2))) - (0.1111111111111111*M2*(M2 + \
2.*t)*pow(std::log(un*(-1.*M2)/pow(muR,2)),2))/pow(M_PI,2) - \
(0.2222222222222222*pow(M2,2)*(1.*M2 + \
2.*t)*std::log(un*pow(muR,2)/pow(Mt,2)))/(pow(M_PI,2)*(1.*M2 - \
4.*pow(Mt,2))) - (0.2222222222222222*M2*(M2 + \
2.*t)*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(un*pow(muR,2)/pow(Mt,2)\
))/pow(M_PI,2) + (0.1111111111111111*M2*(M2 + \
2.*t)*pow(std::log(un*pow(muR,2)/pow(Mt,2)),2))/pow(M_PI,2) - \
(0.33333333333333337*M2*(1.*M2 + \
0.6666666666666665*t)*std::log(un*(-1.*pow(Mt,2))/t))/pow(M_PI,2) - \
(0.2222222222222222*M2*(M2 + 2.*(pow(Mt,2) + \
t))*std::log(un*pow(muR,2)/pow(Mt,2))*std::log(un*(-1.*pow(Mt,2))/t))/\
pow(M_PI,2) - (0.1111111111111111*M2*(M2 + 2.*(pow(Mt,2) + \
t))*pow(std::log(un*(-1.*pow(Mt,2))/t),2))/pow(M_PI,2) - \
(0.3333333333333333*M2*(1.*M2*pow(Mt,2) + 1.6666666666666667*M2*t + \
0.6666666666666666*pow(Mt,2)*t)*std::log(un*(-1.*t)/pow(Mt,2)))/(pow(\
M_PI,2)*(pow(Mt,2) + t)) - (0.2222222222222222*(pow(M2,2) + \
2.*M2*pow(Mt,2) + 2.*M2*t + \
4.*pow(t,2))*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(un*(-1.*t)/pow(\
Mt,2)))/pow(M_PI,2) + (0.2222222222222222*(pow(M2,2) + \
2.*M2*pow(Mt,2) + 2.*M2*t + \
4.*pow(t,2))*std::log(un*pow(muR,2)/pow(Mt,2))*std::log(un*(-1.*t)/\
pow(Mt,2)))/pow(M_PI,2) + (0.1111111111111111*M2*(1.*M2 - \
2.*t)*std::log(un*pow(Mt,2)/(M2 + t)))/pow(M_PI,2) - \
(0.2222222222222222*M2*(M2 - 2.*pow(Mt,2) + \
2.*t)*std::log(un*pow(muR,2)/pow(Mt,2))*std::log(un*pow(Mt,2)/(M2 + \
t)))/pow(M_PI,2) - (0.1111111111111111*M2*(M2 - 2.*pow(Mt,2) + \
2.*t)*pow(std::log(un*pow(Mt,2)/(M2 + t)),2))/pow(M_PI,2) + \
(0.5555555555555556*M2*(1.*pow(M2,2) - 0.2*M2*pow(Mt,2) + \
1.0000000000000002*M2*t + \
0.3999999999999999*pow(Mt,2)*t)*std::log(un*(M2 + \
t)/pow(Mt,2)))/(pow(M_PI,2)*(M2 - 1.*pow(Mt,2) + t)) + \
(0.6666666666666666*(1.*pow(M2,2) + 0.6666666666666666*M2*pow(Mt,2) + \
2.*M2*t + \
1.3333333333333333*pow(t,2))*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(\
un*(M2 + t)/pow(Mt,2)))/pow(M_PI,2) - \
(0.6666666666666666*(1.*pow(M2,2) + 0.6666666666666666*M2*pow(Mt,2) + \
2.*M2*t + \
1.3333333333333333*pow(t,2))*std::log(un*pow(muR,2)/pow(Mt,2))*\
std::log(un*(M2 + t)/pow(Mt,2)))/pow(M_PI,2) - \
(0.4444444444444444*(pow(M2,2) - 8.*M2*pow(Mt,2) + 8.*pow(Mt,4))*(M2 \
+ 2.*t)*gsl_sf_dilog((-1. + 1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + \
Sqrt(1. - (4.*pow(Mt,2))/M2))))/(pow(M_PI,2)*(1.*M2 - \
4.*pow(Mt,2))*Sqrt(1 - (4.*pow(Mt,2))/M2)) + \
(0.2222222222222222*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*gsl_sf_dilog((M2 - \
1.*pow(Mt,2) + t)/(M2 + t)))/pow(M_PI,2) + (0.2222222222222222*M2*(M2 \
+ 2.*(pow(Mt,2) + t))*gsl_sf_dilog((pow(Mt,2) + \
t)/t))/pow(M_PI,2))/pow(M2,2);
	}	  
	  else if(m==1)
	    {
	      resi=((M2*(pow(M2,5)*(17.878510777812895*Sqrt(1 - (4.*pow(Mt,2))/M2) + \
pow(M_PI,2)*(0.4166666666666667 + 0.9444444444444444*Sqrt(1 - \
(4.*pow(Mt,2))/M2))) + pow(Mt,6)*(-3.5555555555555554*pow(M_PI,2) - \
85.33333333333334*Sqrt(1 - (4.*pow(Mt,2))/M2))*pow(t,2) + \
pow(M2,4)*(-123.27106466687738*Sqrt(1 - \
(4.*pow(Mt,2))/M2)*(1.*pow(Mt,2) - 0.1819055289486572*t) + \
pow(M_PI,2)*(pow(Mt,2)*(-2.2777777777777777 - \
6.555555555555555*Sqrt(1 - (4.*pow(Mt,2))/M2)) + (0.8333333333333334 \
+ 1.8888888888888888*Sqrt(1 - (4.*pow(Mt,2))/M2))*t)) + \
pow(M2,2)*pow(Mt,2)*(Sqrt(1 - \
(4.*pow(Mt,2))/M2)*(262.7790115566793*pow(Mt,4) + \
305.44567822334596*pow(Mt,2)*t - 160.722839111673*pow(t,2)) + \
pow(M_PI,2)*(pow(Mt,4)*(14.222222222222221 + 16.*Sqrt(1 - \
(4.*pow(Mt,2))/M2)) + pow(Mt,2)*(16.444444444444443 + \
30.22222222222222*Sqrt(1 - (4.*pow(Mt,2))/M2))*t + (-8.11111111111111 \
- 23.999999999999996*Sqrt(1 - (4.*pow(Mt,2))/M2))*pow(t,2))) + \
pow(M2,3)*(Sqrt(1 - (4.*pow(Mt,2))/M2)*(146.66666666666666*pow(Mt,4) \
- 160.722839111673*pow(Mt,2)*t + 22.423688222292455*pow(t,2)) + \
pow(M_PI,2)*(pow(Mt,4)*(2.111111111111111 + 7.111111111111111*Sqrt(1 \
- (4.*pow(Mt,2))/M2)) + pow(Mt,2)*(-4.777777777777778 - \
15.11111111111111*Sqrt(1 - (4.*pow(Mt,2))/M2))*t + \
(1.1111111111111112 + 2.9999999999999996*Sqrt(1 - \
(4.*pow(Mt,2))/M2))*pow(t,2))) + M2*pow(Mt,4)*(Sqrt(1 - \
(4.*pow(Mt,2))/M2)*(-85.33333333333334*pow(Mt,4) - \
85.33333333333334*pow(Mt,2)*t + 305.44567822334596*pow(t,2)) + \
pow(M_PI,2)*(-3.5555555555555554*pow(Mt,4) + \
5.333333333333334*pow(Mt,2)*t + (27.55555555555555 + \
47.99999999999999*Sqrt(1 - (4.*pow(Mt,2))/M2))*pow(t,2))) + M2*Sqrt(1 \
- (4.*pow(Mt,2))/M2)*(4.*pow(M2,4) + \
128.00000000000006*pow(Mt,4)*pow(t,2) + \
pow(M2,2)*t*(-64.00000000000003*pow(Mt,2) + 8.000000000000004*t) + \
pow(M2,3)*(-23.999999999999986*pow(Mt,2) + 8.000000000000004*t) + \
M2*(128.00000000000006*pow(Mt,6) + 128.00000000000006*pow(Mt,4)*t - \
64.00000000000003*pow(Mt,2)*pow(t,2)))*std::log(un*M_PI) + \
pow(M2,2)*pow(Mt,2)*Sqrt(1 - \
(4.*pow(Mt,2))/M2)*(1.4210854715202004e-14*M2 + \
1.4210854715202004e-14*t)*t*pow(std::log(un*M_PI),2)))/(pow(1.*M2 - \
4.*pow(Mt,2),2)*Sqrt(1 - (4.*pow(Mt,2))/M2)) - \
2.25*pow(M2,4)*std::log(un*(-1.*M2)/pow(Mt,2)) - \
(0.08333333333333333*(M2 + Sqrt(M2*(M2 - \
4.*pow(Mt,2))))*(1.*pow(M2,4) - 64.*pow(Mt,4)*pow(t,2) + \
pow(M2,3)*(6.*pow(Mt,2) + 2.*t) + M2*(-64.*pow(Mt,6) - \
64.*pow(Mt,4)*t) + pow(M2,2)*(-32.*pow(Mt,4) + \
2.*pow(t,2)))*std::log(un*(-1.*M2 + 1.*Sqrt(M2*(M2 - \
4.*pow(Mt,2))))/(M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2))))))/(1.*M2 - \
4.*pow(Mt,2)) + (0.08333333333333333*(-1.*M2 + Sqrt(M2*(M2 - \
4.*pow(Mt,2))))*(1.*pow(M2,4) - 64.*pow(Mt,4)*pow(t,2) + \
pow(M2,3)*(6.*pow(Mt,2) + 2.*t) + M2*(-64.*pow(Mt,6) - \
64.*pow(Mt,4)*t) + pow(M2,2)*(-32.*pow(Mt,4) + \
2.*pow(t,2)))*std::log(un*(1.*M2 + 1.*Sqrt(M2*(M2 - \
4.*pow(Mt,2))))/(-1.*M2 + 1.*Sqrt(M2*(M2 - 4.*pow(Mt,2))))))/(1.*M2 - \
4.*pow(Mt,2)) + (0.16666666666666666*M2*(1.*pow(M2,4) - \
64.*pow(Mt,4)*pow(t,2) + pow(M2,3)*(6.*pow(Mt,2) + 2.*t) + \
M2*(-64.*pow(Mt,6) - 64.*pow(Mt,4)*t) + pow(M2,2)*(-32.*pow(Mt,4) + \
						   2.*pow(t,2)))*std::log(un*-0.5 - (0.5*Sqrt(M2*(M2 - \
												  4.*pow(Mt,2))))/M2))/(1.*M2 - 4.*pow(Mt,2)) + \
(0.16666666666666666*M2*(1.*pow(M2,4) - 64.*pow(Mt,4)*pow(t,2) + \
pow(M2,3)*(6.*pow(Mt,2) + 2.*t) + M2*(-64.*pow(Mt,6) - \
64.*pow(Mt,4)*t) + pow(M2,2)*(-32.*pow(Mt,4) + \
			      2.*pow(t,2)))*std::log(un*-0.5 + (0.5*Sqrt(M2*(M2 - \
									     4.*pow(Mt,2))))/M2))/(1.*M2 - 4.*pow(Mt,2)) + \
(1.0833333333333333*pow(M2,2)*(1.*pow(M2,4) + \
pow(M2,3)*(-5.076923076923077*pow(Mt,2) + 2.*t) + \
pow(Mt,4)*t*(24.615384615384617*pow(Mt,2) + 66.46153846153847*t) + \
pow(M2,2)*(4.*pow(Mt,4) - 10.153846153846153*pow(Mt,2)*t + \
2.769230769230769*pow(t,2)) + M2*(34.46153846153846*pow(Mt,6) + \
35.69230769230769*pow(Mt,4)*t - \
19.384615384615383*pow(Mt,2)*pow(t,2)))*pow(std::log(un*(1. - \
1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - \
(4.*pow(Mt,2))/M2))),2))/(pow(1.*M2 - 4.*pow(Mt,2),2)*Sqrt(1 - \
(4.*pow(Mt,2))/M2)) - (0.3333333333333333*pow(M2,3)*(M2 - \
2.*pow(Mt,2))*std::log(un*(-1. + 1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. \
+ Sqrt(1. - (4.*pow(Mt,2))/M2))))/Sqrt(1 - (4.*pow(Mt,2))/M2) - \
(0.6666666666666666*M2*(M2 - 2.*pow(Mt,2))*(pow(M2,2) + \
2.*M2*pow(Mt,2) + 2.*M2*t + 2.*pow(t,2))*std::log(un*(4.*Sqrt(1. - \
(4.*pow(Mt,2))/M2))/pow(1. + Sqrt(1. - \
(4.*pow(Mt,2))/M2),2))*std::log(un*(-1. + 1.*Sqrt(1. - \
(4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2))))/Sqrt(1 - \
(4.*pow(Mt,2))/M2) + (0.16666666666666666*M2*(M2 - \
2.*pow(Mt,2))*(pow(M2,2) + 2.*M2*pow(Mt,2) + 2.*M2*t + \
2.*pow(t,2))*pow(std::log(un*(-1. + 1.*Sqrt(1. - \
(4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2))),2))/Sqrt(1 - \
(4.*pow(Mt,2))/M2) + (2.5*pow(M2,2)*(0.43333333333333335*pow(M2,4) - \
6.533333333333333*pow(M2,3)*pow(Mt,2) - \
8.533333333333339*pow(Mt,4)*pow(t,2) + \
pow(M2,2)*(6.933333333333334*pow(Mt,4) - \
12.26666666666667*pow(Mt,2)*t + 0.6666666666666665*pow(t,2)) + \
M2*(-8.533333333333331*pow(Mt,6) - 8.533333333333331*pow(Mt,4)*t - \
14.933333333333334*pow(Mt,2)*pow(t,2)))*std::log(un*(-1.*M2)/pow(muR,\
2)))/pow(1.*M2 - 4.*pow(Mt,2),2) + 1.25*pow(M2,2)*(1.*pow(M2,2) + \
1.4666666666666666*M2*pow(Mt,2) + 2.*M2*t + \
2.6666666666666665*pow(t,2))*pow(std::log(un*(-1.*M2)/pow(muR,2)),2) \
+ (0.16666666666666666*M2*(1.*pow(M2,4) - 64.*pow(Mt,4)*pow(t,2) + \
pow(M2,3)*(6.*pow(Mt,2) + 2.*t) + M2*(-64.*pow(Mt,6) - \
64.*pow(Mt,4)*t) + pow(M2,2)*(-32.*pow(Mt,4) + \
2.*pow(t,2)))*std::log(un*M2/pow(muR,2)))/(1.*M2 - 4.*pow(Mt,2)) - \
(0.3333333333333333*M2*(M2 - 2.*pow(Mt,2))*(pow(M2,2) + \
2.*M2*pow(Mt,2) + 2.*M2*t + 2.*pow(t,2))*std::log(un*(-1. + \
1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - \
(4.*pow(Mt,2))/M2)))*std::log(un*pow(Mt,2)/pow(muR,2)))/Sqrt(1 - \
(4.*pow(Mt,2))/M2) + \
(5.166666666666666*M2*(1.7258064516129032*pow(M2,5) + \
8.258064516129032*pow(Mt,6)*pow(t,2) + \
pow(M2,4)*(-13.032258064516132*pow(Mt,2) + 2.*t) + \
pow(M2,3)*(9.806451612903226*pow(Mt,4) - \
21.677419354838708*pow(Mt,2)*t + 2.3225806451612905*pow(t,2)) + \
pow(M2,2)*(28.903225806451612*pow(Mt,6) + \
24.774193548387096*pow(Mt,4)*t - \
22.967741935483872*pow(Mt,2)*pow(t,2)) + \
M2*(8.258064516129032*pow(Mt,8) + 8.258064516129032*pow(Mt,6)*t + \
24.774193548387096*pow(Mt,4)*pow(t,2)))*std::log(un*pow(muR,2)/pow(Mt,\
2)))/pow(1.*M2 - 4.*pow(Mt,2),2) + \
2.1666666666666665*pow(M2,2)*(1.*pow(M2,2) + \
1.3846153846153846*M2*pow(Mt,2) + 2.*M2*t + \
2.769230769230769*pow(t,2))*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(\
un*pow(muR,2)/pow(Mt,2)) - 2.583333333333333*pow(M2,2)*(1.*pow(M2,2) \
+ 1.7419354838709677*M2*pow(Mt,2) + 2.*M2*t + \
2.3225806451612905*pow(t,2))*pow(std::log(un*pow(muR,2)/pow(Mt,2)),2) \
- 4.*pow(M2,2)*(pow(M2,2) + 2.*M2*pow(Mt,2) + 2.*M2*t + \
2.*pow(t,2))*std::log(un*(4.*M_PI*pow(muR,2))/pow(Mt,2)) - \
3.500000000000001*pow(M2,3)*(1.*M2 + \
0.6666666666666665*t)*std::log(un*(-1.*pow(Mt,2))/t) - \
2.3333333333333335*pow(M2,3)*(M2 + 2.*(pow(Mt,2) + \
t))*std::log(un*pow(muR,2)/pow(Mt,2))*std::log(un*(-1.*pow(Mt,2))/t) \
- 1.1666666666666667*pow(M2,3)*(M2 + 2.*(pow(Mt,2) + \
t))*pow(std::log(un*(-1.*pow(Mt,2))/t),2) - \
(3.5*pow(M2,3)*(1.*M2*pow(Mt,2) + 1.6666666666666667*M2*t + \
0.6666666666666666*pow(Mt,2)*t)*std::log(un*(-1.*t)/pow(Mt,2)))/(pow(\
Mt,2) + t) - 2.3333333333333335*pow(M2,2)*(pow(M2,2) + \
2.*M2*pow(Mt,2) + 2.*M2*t + \
4.*pow(t,2))*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(un*(-1.*t)/pow(\
Mt,2)) + 2.3333333333333335*pow(M2,2)*(pow(M2,2) + 2.*M2*pow(Mt,2) + \
2.*M2*t + \
4.*pow(t,2))*std::log(un*pow(muR,2)/pow(Mt,2))*std::log(un*(-1.*t)/\
pow(Mt,2)) - 0.3333333333333333*pow(M2,3)*(1.*M2 - \
2.*t)*std::log(un*pow(Mt,2)/(M2 + t)) + \
0.6666666666666666*pow(M2,3)*(M2 - 2.*pow(Mt,2) + \
2.*t)*std::log(un*pow(muR,2)/pow(Mt,2))*std::log(un*pow(Mt,2)/(M2 + \
t)) + 0.3333333333333333*pow(M2,3)*(M2 - 2.*pow(Mt,2) + \
2.*t)*pow(std::log(un*pow(Mt,2)/(M2 + t)),2) - \
(1.6666666666666667*pow(M2,3)*(1.*pow(M2,2) - 0.2*M2*pow(Mt,2) + \
1.0000000000000002*M2*t + \
0.3999999999999999*pow(Mt,2)*t)*std::log(un*(M2 + t)/pow(Mt,2)))/(M2 \
- 1.*pow(Mt,2) + t) - 2.*pow(M2,2)*(1.*pow(M2,2) + \
0.6666666666666666*M2*pow(Mt,2) + 2.*M2*t + \
1.3333333333333333*pow(t,2))*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(\
un*(M2 + t)/pow(Mt,2)) + 2.*pow(M2,2)*(1.*pow(M2,2) + \
0.6666666666666666*M2*pow(Mt,2) + 2.*M2*t + \
1.3333333333333333*pow(t,2))*std::log(un*pow(muR,2)/pow(Mt,2))*\
std::log(un*(M2 + t)/pow(Mt,2)) - (0.3333333333333333*M2*(M2 - \
2.*pow(Mt,2))*(pow(M2,2) + 2.*M2*pow(Mt,2) + 2.*M2*t + \
2.*pow(t,2))*gsl_sf_dilog(pow(1. - 1.*Sqrt(1. - \
(4.*pow(Mt,2))/M2),2)/pow(1. + Sqrt(1. - \
(4.*pow(Mt,2))/M2),2)))/Sqrt(1 - (4.*pow(Mt,2))/M2) - \
(0.6666666666666666*M2*(M2 - 2.*pow(Mt,2))*(pow(M2,2) + \
2.*M2*pow(Mt,2) + 2.*M2*t + 2.*pow(t,2))*gsl_sf_dilog(2./(1. + \
Sqrt(1. - (4.*pow(Mt,2))/M2))))/Sqrt(1 - (4.*pow(Mt,2))/M2) + \
(4.333333333333333*pow(M2,2)*(1.*pow(M2,4) + \
pow(M2,3)*(-5.076923076923077*pow(Mt,2) + 2.*t) + \
pow(Mt,4)*t*(24.615384615384617*pow(Mt,2) + 66.46153846153847*t) + \
pow(M2,2)*(4.*pow(Mt,4) - 10.153846153846153*pow(Mt,2)*t + \
2.769230769230769*pow(t,2)) + M2*(34.46153846153846*pow(Mt,6) + \
35.69230769230769*pow(Mt,4)*t - \
19.384615384615383*pow(Mt,2)*pow(t,2)))*gsl_sf_dilog((-1. + \
1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - \
(4.*pow(Mt,2))/M2))))/(pow(1.*M2 - 4.*pow(Mt,2),2)*Sqrt(1 - \
(4.*pow(Mt,2))/M2)) - 0.6666666666666666*pow(M2,3)*(M2 - 2.*pow(Mt,2) \
+ 2.*t)*gsl_sf_dilog((M2 - 1.*pow(Mt,2) + t)/(M2 + t)) + \
2.3333333333333335*pow(M2,3)*(M2 + 2.*(pow(Mt,2) + \
t))*gsl_sf_dilog((pow(Mt,2) + t)/t))/(pow(M2,4)*pow(M_PI,2));
	    }
	}
      res=2.*std::real(resi)*pow(alpha_s*4.*M_PI,3)/dR/dR/2./M2/Nsym;
    }

  else if (chan==1)
    {
      double tMand=-(t+pow(Mt,2)+M2)+2.*pow(Mt,2), mu2=pow(muR,2), Mt2=pow(Mt,2.);
      double tmp[3][2];
      int nn=n+1, mm=m+1;
      if(abs(t+pow(Mt,2))<abs(tMand)) tMand=t+pow(Mt,2);
      gg2ttx_shfmatrix_(mm,nn,M2,tMand,Mt2,mu2,tmp); //here we have directly the H1_gg matrix with poles
      resi=(tmp[0][0]+i*tmp[0][1]); //We can take only the real part of H1 because we will trace it with S0 which is symmetric 

      //std::cout << std::setprecision(20);
      //std::cout << tmp[0][0] << " "  << tmp[0][1] << " "  << tmp[1][0] << " "  << tmp[1][1] << " "  << tmp[2][0] << " " << tmp[2][1]  << std::endl;
      //std::cout << " n = " << nn << " m = " << mm << " mu = " << muR << " t = "  << tMand << " s = " << M2 << " Mt = " << Mt <<  std::endl;
      //std::cout << tmp[0][0] << " + i * "  << tmp[0][1] << std::endl;
      res=resi*pow(alpha_s,3)*M_PI/dR/dR/2./M2/Nsym;
    }  
  return res;
}
