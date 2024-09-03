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
// -------- Classes ----------------------------------- //
#include "process.h"	// Process//
// -------- Functions --------------------------------- //
// ---------------------------------------------------- //
extern double MZ, Mt, alpha_s, muF, muR, beta0, beta1, CF, C, Phi;
extern int Nc;

double Xcos(double);

std::complex<double> Sqrt(double x)
{
  return pow(x,0.5);
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
      if(m*n==1) res=2.*((pow(t,2)+pow(u,2))/pow(M2,2)+2.*Mt*Mt/M2); //see 0805.1885 or SCET 1003.5827
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


double H1(double *x, int m, int n, int chan, double M2)
{
  const int Len=chan+2;

  std::complex<double> H[Len][Len];
  std::complex<double> i(0.0,1.0), un(1.0,0.0);
  double Xbet=0, Xcos=0, res=0, bet=0, t=0,dR=0.;
  
  bet=pow(1.-4.*Mt*Mt/M2,0.5);
  
  Xcos=2.*x[1]-1.;
  Xbet=Xcos*bet;
  t=-M2/2.*(1.-Xbet)+pow(Mt,2);
  
  if(chan==0) dR=double(Nc)*2.;
  else if(chan==1) dR=(pow(Nc,2)-1.)*2.;

  
  if(chan==0)
    {
  
      if(n==0)
	{
	  if(m==0) H[0][0]=0;
	  else if(m==1)  H[1][0]=1./pow(M2,2)*(-0.009259259259259259*(((8.*pow(Mt,4) - 8.*M2*pow(Mt,2)*(1. + 2.*Sqrt(1. - (4.*pow(Mt,2))/M2)) + pow(M2,2)*(1. + 4.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*          (M2 - 2.*pow(Mt,2) + 2.*t))/((M2 - 4.*pow(Mt,2))*Sqrt(1. - (4.*pow(Mt,2))/M2)) +					1.5386782351823112*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*std::log(un*(-1.*M2)/pow(Mt,2))       +        (0.3039635509270133*(pow(M2,2) - 8.*M2*pow(Mt,2) + 8.*pow(Mt,4))*(M2 - 2.*pow(Mt,2) + 2.*t)*          pow(std::log(un*(1. - 1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2))),2))/((M2 - 4.*pow(Mt,2))*Sqrt(1. - (4.*pow(Mt,2))/M2)) +        0.3039635509270133*M2*(8.062048493938581*M2 - 2.*pow(Mt,2) + 12.124096987877163*t)*std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t)) +        0.3039635509270133*M2*(M2 + 2.*t)*pow(std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t)),2) +        0.3039635509270133*M2*(4.0620484939385815*M2 - 22.248193975754326*pow(Mt,2) + 12.124096987877163*t)*std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t)) +        0.3039635509270133*M2*(M2 - 4.*pow(Mt,2) + 2.*t)*pow(std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t)),2) -        (0.3039635509270133*(-10.186145481815743*pow(M2,3) + pow(M2,2)*(44.620484939385804*pow(Mt,2) - 40.55843644544723*t) +             20.248193975754326*pow(pow(Mt,2) - 1.*t,2)*(2.*pow(Mt,2) - 1.*t) +             2.*M2*(-31.372290963631485*pow(Mt,4) + 61.74458192726297*pow(Mt,2)*t - 25.31024246969291*pow(t,2)))*std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2)))/        (M2 - 2.*pow(Mt,2) + t) - 0.6079271018540267*(3.*pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + M2*(-4.*pow(Mt,2) + 6.*t))*std::log(un*(-1.*M2)/pow(Mt,2))*				   std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2)) - (0.3039635509270133*          (pow(M2,2)*(2.*pow(Mt,2) + 0.06204849393858147*t) + 20.248193975754326*pow(pow(Mt,2) - 1.*t,2)*t +             2.*M2*(pow(Mt,4) - 1.*pow(Mt,2)*t + 5.0620484939385815*pow(t,2)))*std::log(1. - (1.*t)/pow(Mt,2)))/t +        0.6079271018540267*(pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(1. - (1.*t)/pow(Mt,2)) +        (0.6079271018540267*M2*(-1.5310242469692907*M2 + 10.124096987877163*pow(Mt,2))*(M2 - 2.*pow(Mt,2) + 2.*t)*std::log(un*(-12.566370614359172*M2)/pow(muR,2)))/        (M2 - 4.*pow(Mt,2)) + 0.3039635509270133*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*pow(std::log(un*(-12.566370614359172*M2)/pow(muR,2)),2) +        (0.6079271018540267*M2*(-1.5310242469692907*M2 + 10.124096987877163*pow(Mt,2))*(M2 - 2.*pow(Mt,2) + 2.*t)*          std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)))/(M2 - 4.*pow(Mt,2)) +        0.6079271018540267*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) +        0.6079271018540267*M2*(M2 + 2.*t)*std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) +        0.6079271018540267*M2*(M2 - 4.*pow(Mt,2) + 2.*t)*std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) +        0.6079271018540267*(3.*pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + M2*(-4.*pow(Mt,2) + 6.*t))*std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2))*        std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) - 0.6079271018540267*(pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*				   std::log(1. - (1.*t)/pow(Mt,2))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) - 				   0.3039635509270133*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*pow(std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)),2) +				   (1.2158542037080533*(pow(M2,2) - 8.*M2*pow(Mt,2) + 8.*pow(Mt,4))*(M2 - 2.*pow(Mt,2) + 2.*t)*           gsl_sf_dilog((-1. + sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + sqrt(1. - (4.*pow(Mt,2))/M2))))/((M2 - 4.*pow(Mt,2))*sqrt(1. - (4.*pow(Mt,2))/M2)) -        0.6079271018540267*M2*(M2 + 2.*t)* gsl_sf_dilog((-1.*t)/(pow(Mt,2) - 1.*t)) -        0.6079271018540267*M2*(M2 - 4.*pow(Mt,2) + 2.*t)* gsl_sf_dilog((M2 - 2.*pow(Mt,2) + t)/(M2 - 1.*pow(Mt,2) + t))));

	}

      else if (n==1)
	{
	  if(m==0) H[0][1]=(-0.009259259259259259*(((8.*pow(Mt,4) - 8.*M2*pow(Mt,2)*(1. + 2.*Sqrt(1. - (4.*pow(Mt,2))/M2)) + pow(M2,2)*(1. + 4.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*
          (M2 - 2.*pow(Mt,2) + 2.*t))/((M2 - 4.*pow(Mt,2))*Sqrt(1. - (4.*pow(Mt,2))/M2)) + 
       1.5386782351823112*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*std::log(un*(-1.*M2)/pow(Mt,2)) + 
       (0.3039635509270133*(pow(M2,2) - 8.*M2*pow(Mt,2) + 8.*pow(Mt,4))*(M2 - 2.*pow(Mt,2) + 2.*t)*
          pow(std::log(un*(1. - 1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2))),2))/((M2 - 4.*pow(Mt,2))*Sqrt(1. - (4.*pow(Mt,2))/M2)) + 
       0.3039635509270133*M2*(8.062048493938581*M2 - 2.*pow(Mt,2) + 12.124096987877163*t)*std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t)) + 
       0.3039635509270133*M2*(M2 + 2.*t)*pow(std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t)),2) + 
       0.3039635509270133*M2*(4.0620484939385815*M2 - 22.248193975754326*pow(Mt,2) + 12.124096987877163*t)*std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t)) + 
       0.3039635509270133*M2*(M2 - 4.*pow(Mt,2) + 2.*t)*pow(std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t)),2) - 
       (0.3039635509270133*(-10.186145481815743*pow(M2,3) + pow(M2,2)*(44.620484939385804*pow(Mt,2) - 40.55843644544723*t) + 
            20.248193975754326*pow(pow(Mt,2) - 1.*t,2)*(2.*pow(Mt,2) - 1.*t) + 
            2.*M2*(-31.372290963631485*pow(Mt,4) + 61.74458192726297*pow(Mt,2)*t - 25.31024246969291*pow(t,2)))*std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2)))/
        (M2 - 2.*pow(Mt,2) + t) - 0.6079271018540267*(3.*pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + M2*(-4.*pow(Mt,2) + 6.*t))*std::log(un*(-1.*M2)/pow(Mt,2))*
        std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2)) - (0.3039635509270133*
          (pow(M2,2)*(2.*pow(Mt,2) + 0.06204849393858147*t) + 20.248193975754326*pow(pow(Mt,2) - 1.*t,2)*t + 
            2.*M2*(pow(Mt,4) - 1.*pow(Mt,2)*t + 5.0620484939385815*pow(t,2)))*std::log(1. - (1.*t)/pow(Mt,2)))/t + 
       0.6079271018540267*(pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(1. - (1.*t)/pow(Mt,2)) + 
       (0.6079271018540267*M2*(-1.5310242469692907*M2 + 10.124096987877163*pow(Mt,2))*(M2 - 2.*pow(Mt,2) + 2.*t)*std::log(un*(-12.566370614359172*M2)/pow(muR,2)))/
        (M2 - 4.*pow(Mt,2)) + 0.3039635509270133*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*pow(std::log(un*(-12.566370614359172*M2)/pow(muR,2)),2) + 
       (0.6079271018540267*M2*(-1.5310242469692907*M2 + 10.124096987877163*pow(Mt,2))*(M2 - 2.*pow(Mt,2) + 2.*t)*
          std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)))/(M2 - 4.*pow(Mt,2)) + 
       0.6079271018540267*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) + 
       0.6079271018540267*M2*(M2 + 2.*t)*std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) + 
       0.6079271018540267*M2*(M2 - 4.*pow(Mt,2) + 2.*t)*std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) + 
       0.6079271018540267*(3.*pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + M2*(-4.*pow(Mt,2) + 6.*t))*std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2))*
        std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) - 0.6079271018540267*(pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*
        std::log(1. - (1.*t)/pow(Mt,2))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) - 
       0.3039635509270133*M2*(M2 - 2.*pow(Mt,2) + 2.*t)*pow(std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)),2) + 
       (1.2158542037080533*(pow(M2,2) - 8.*M2*pow(Mt,2) + 8.*pow(Mt,4))*(M2 - 2.*pow(Mt,2) + 2.*t)*
           gsl_sf_dilog((-1. + sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + sqrt(1. - (4.*pow(Mt,2))/M2))))/((M2 - 4.*pow(Mt,2))*sqrt(1. - (4.*pow(Mt,2))/M2)) - 
       0.6079271018540267*M2*(M2 + 2.*t)* gsl_sf_dilog((-1.*t)/(pow(Mt,2) - 1.*t)) - 
       0.6079271018540267*M2*(M2 - 4.*pow(Mt,2) + 2.*t)* gsl_sf_dilog((M2 - 2.*pow(Mt,2) + t)/(M2 - 1.*pow(Mt,2) + t))))/pow(M2,2);

	  else if(m==1)  H[1][1]=(0.0007036193308495678*((pow(M2,6)*(1138.7935823400087*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(15. + 42.*Sqrt(1. - (4.*pow(Mt,2))/M2))) - 
          128.*M2*pow(Mt,6)*(9.869604401089358 + 24.*Sqrt(1. - (4.*pow(Mt,2))/M2))*pow(pow(Mt,2) - 1.*t,2) + 
          pow(M2,5)*(-16.*pow(Mt,2)*(575.3967911700043*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(7. + 23.*Sqrt(1. - (4.*pow(Mt,2))/M2))) + 
             6.*(218.60508487698564*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(5. + 14.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*t) + 
          pow(M2,4)*(pow(Mt,4)*(20300.327826702058*Sqrt(1. - (4.*pow(Mt,2))/M2) + 39.47841760435743*(72. + 263.*Sqrt(1. - (4.*pow(Mt,2))/M2))) - 
             4.*pow(Mt,2)*(3111.0762731547843*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(63. + 230.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*t + 
             4.*(327.90762731547846*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(10. + 31.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*pow(t,2)) + 
          32.*pow(M2,2)*(pow(Mt,8)*(595.8152546309569*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(21. + 62.*Sqrt(1. - (4.*pow(Mt,2))/M2))) - 
             4.*pow(Mt,6)*(321.90762731547846*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(14. + 31.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*t + 
             pow(Mt,4)*(595.8152546309569*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(31. + 62.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*pow(t,2)) - 
          4.*pow(M2,3)*pow(Mt,2)*(pow(Mt,4)*(2839.2610185238277*Sqrt(1. - (4.*pow(Mt,2))/M2) + 
                9.869604401089358*(93. + 376.*Sqrt(1. - (4.*pow(Mt,2))/M2))) - 
             2.*pow(Mt,2)*(4838.522037047655*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(147. + 416.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*t + 
             (2455.2610185238277*Sqrt(1. - (4.*pow(Mt,2))/M2) + 9.869604401089358*(73. + 248.*Sqrt(1. - (4.*pow(Mt,2))/M2)))*pow(t,2)))/
        (pow(M2 - 4.*pow(Mt,2),2)*Sqrt(1. - (4.*pow(Mt,2))/M2)) - 
       3.*pow(M2,2)*(-38.80663042120155*pow(M2,2) + 10.124096987877163*M2*(4.*pow(Mt,2) - 13.*t) - 182.23374578178894*pow(pow(Mt,2) - 1.*t,2))*
        std::log(un*(-1.*M2)/pow(Mt,2)) + (6.*M2*(pow(M2,4) - 64.*pow(Mt,4)*pow(pow(Mt,2) - 1.*t,2) - 64.*M2*pow(Mt,4)*t + 
            2.*pow(M2,3)*(2.*pow(Mt,2) + t) + pow(M2,2)*(-30.*pow(Mt,4) - 4.*pow(Mt,2)*t + 2.*pow(t,2)))*
          std::log(un*(0.5*(-1.*M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2)))))/M2))/(M2 - 4.*pow(Mt,2)) - 
       (3.*(M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2))))*(pow(M2,4) - 64.*pow(Mt,4)*pow(pow(Mt,2) - 1.*t,2) - 64.*M2*pow(Mt,4)*t + 
            2.*pow(M2,3)*(2.*pow(Mt,2) + t) + pow(M2,2)*(-30.*pow(Mt,4) - 4.*pow(Mt,2)*t + 2.*pow(t,2)))*
          std::log(un*(-1.*M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2))))/(M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2))))))/(M2 - 4.*pow(Mt,2)) + 
       (6.*M2*(pow(M2,4) - 64.*pow(Mt,4)*pow(pow(Mt,2) - 1.*t,2) - 64.*M2*pow(Mt,4)*t + 2.*pow(M2,3)*(2.*pow(Mt,2) + t) + 
            pow(M2,2)*(-30.*pow(Mt,4) - 4.*pow(Mt,2)*t + 2.*pow(t,2)))*std::log(un*(-0.5*(M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2)))))/M2))/(M2 - 4.*pow(Mt,2)) + 
       (3.*(-1.*M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2))))*(pow(M2,4) - 64.*pow(Mt,4)*pow(pow(Mt,2) - 1.*t,2) - 64.*M2*pow(Mt,4)*t + 
            2.*pow(M2,3)*(2.*pow(Mt,2) + t) + pow(M2,2)*(-30.*pow(Mt,4) - 4.*pow(Mt,2)*t + 2.*pow(t,2)))*
          std::log(un*(M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2))))/(-1.*M2 + Sqrt(M2*(M2 - 4.*pow(Mt,2))))))/(M2 - 4.*pow(Mt,2)) + 
       (3.*pow(M2,2)*(13.*pow(M2,4) + pow(M2,3)*(-92.*pow(Mt,2) + 26.*t) + 4.*pow(M2,2)*(55.*pow(Mt,4) - 51.*pow(Mt,2)*t + 9.*pow(t,2)) - 
            4.*M2*(67.*pow(Mt,6) - 242.*pow(Mt,4)*t + 63.*pow(Mt,2)*pow(t,2)) + 32.*(17.*pow(Mt,8) - 44.*pow(Mt,6)*t + 27.*pow(Mt,4)*pow(t,2)))*
          pow(std::log(un*(1. - 1.*Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2))),2))/
        (pow(M2 - 4.*pow(Mt,2),2)*Sqrt(1. - (4.*pow(Mt,2))/M2)) + 
       (12.*M2*(M2 - 2.*pow(Mt,2))*(1.5310242469692907*pow(M2,2) + 5.0620484939385815*pow(Mt,4) + 5.0620484939385815*M2*t - 
            10.124096987877163*pow(Mt,2)*t + 5.0620484939385815*pow(t,2))*std::log(un*(-1. + Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2))))/
        Sqrt(1. - (4.*pow(Mt,2))/M2) - (24.*M2*(M2 - 2.*pow(Mt,2))*(pow(M2,2) + 2.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*
          std::log(un*(4.*Sqrt(1. - (4.*pow(Mt,2))/M2))/pow(1. + Sqrt(1. - (4.*pow(Mt,2))/M2),2))*
          std::log(un*(-1. + Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2))))/Sqrt(1. - (4.*pow(Mt,2))/M2) + 
       (6.*M2*(M2 - 2.*pow(Mt,2))*(pow(M2,2) + 2.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*
          pow(std::log(un*(-1. + Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2))),2))/Sqrt(1. - (4.*pow(Mt,2))/M2) + 
       42.*pow(M2,3)*(-8.062048493938581*M2 + 2.*pow(Mt,2) - 12.124096987877163*t)*std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t)) - 
       42.*pow(M2,3)*(M2 + 2.*t)*pow(std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t)),2) - 
       12.*pow(M2,3)*(-4.0620484939385815*M2 + 22.248193975754326*pow(Mt,2) - 12.124096987877163*t)*std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t)) + 
       12.*pow(M2,3)*(M2 - 4.*pow(Mt,2) + 2.*t)*pow(std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t)),2) - 
       (12.*pow(M2,2)*(-10.186145481815743*pow(M2,3) + pow(M2,2)*(44.620484939385804*pow(Mt,2) - 40.55843644544723*t) + 
            20.248193975754326*pow(pow(Mt,2) - 1.*t,2)*(2.*pow(Mt,2) - 1.*t) + 
            2.*M2*(-31.372290963631485*pow(Mt,4) + 61.74458192726297*pow(Mt,2)*t - 25.31024246969291*pow(t,2)))*std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2)))/
        (M2 - 2.*pow(Mt,2) + t) - 24.*pow(M2,2)*(3.*pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + M2*(-4.*pow(Mt,2) + 6.*t))*std::log(un*(-1.*M2)/pow(Mt,2))*
        std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2)) + (42.*pow(M2,2)*
          (pow(M2,2)*(2.*pow(Mt,2) + 0.06204849393858147*t) + 20.248193975754326*pow(pow(Mt,2) - 1.*t,2)*t + 
            2.*M2*(pow(Mt,4) - 1.*pow(Mt,2)*t + 5.0620484939385815*pow(t,2)))*std::log(1. - (1.*t)/pow(Mt,2)))/t - 
       84.*pow(M2,2)*(pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*std::log(un*(-1.*M2)/pow(Mt,2))*std::log(1. - (1.*t)/pow(Mt,2)) + 
       (3.*pow(M2,2)*(-62.93072740907872*pow(M2,4) - 3495.7110361206915*pow(Mt,4)*pow(pow(Mt,2) - 1.*t,2) - 
            4.*pow(M2,3)*(-112.9855518060346*pow(Mt,2) + 37.965363704539364*t) + 
            4.*pow(M2,2)*(-286.336170478718*pow(Mt,4) + 302.96387951508643*pow(Mt,2)*t - 45.62048493938582*pow(t,2)) + 
            64.*M2*(28.43433945757007*pow(Mt,6) - 78.58584864392517*pow(Mt,4)*t + 18.31024246969291*pow(Mt,2)*pow(t,2)))*
          std::log(un*(-12.566370614359172*M2)/pow(muR,2)))/pow(M2 - 4.*pow(Mt,2),2) + 
       3.*pow(M2,2)*(15.*pow(M2,2) + 40.*pow(pow(Mt,2) - 1.*t,2) + M2*(-8.*pow(Mt,2) + 30.*t))*pow(std::log(un*(-12.566370614359172*M2)/pow(muR,2)),2) + 
       (6.*M2*(pow(M2,4) - 64.*pow(Mt,4)*pow(pow(Mt,2) - 1.*t,2) - 64.*M2*pow(Mt,4)*t + 2.*pow(M2,3)*(2.*pow(Mt,2) + t) + 
            pow(M2,2)*(-30.*pow(Mt,4) - 4.*pow(Mt,2)*t + 2.*pow(t,2)))*std::log(un*(12.566370614359172*M2)/pow(muR,2)))/(M2 - 4.*pow(Mt,2)) - 
       (12.*M2*(M2 - 2.*pow(Mt,2))*(pow(M2,2) + 2.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*
          std::log(un*(-1. + Sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + Sqrt(1. - (4.*pow(Mt,2))/M2)))*std::log(un*(12.566370614359172*pow(Mt,2))/pow(muR,2)))/
        Sqrt(1. - (4.*pow(Mt,2))/M2) + (3.*M2*(-57.92350331209599*pow(M2,5) + 512.*pow(Mt,6)*pow(pow(Mt,2) - 1.*t,2) - 
            4.*pow(M2,4)*(-106.97110361206921*pow(Mt,2) + 51.461751656048*t) + 
            16.*pow(M2,3)*(-77.95091551057396*pow(Mt,4) + 110.48193975754322*pow(Mt,2)*t - 14.779218222723614*pow(t,2)) - 
            16.*pow(M2,2)*(-157.73013373329758*pow(Mt,6) + 436.3144981877698*pow(Mt,4)*t - 101.23374578178894*pow(Mt,2)*pow(t,2)) + 
            256.*M2*(-17.779218222723614*pow(Mt,8) + 37.558436445447235*pow(Mt,6)*t - 17.779218222723614*pow(Mt,4)*pow(t,2)))*
          std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)))/pow(M2 - 4.*pow(Mt,2),2) + 
       6.*pow(M2,2)*(13.*pow(M2,2) + 36.*pow(pow(Mt,2) - 1.*t,2) + M2*(-8.*pow(Mt,2) + 26.*t))*std::log(un*(-1.*M2)/pow(Mt,2))*
        std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) - 84.*pow(M2,3)*(M2 + 2.*t)*std::log(pow(Mt,2)/(pow(Mt,2) - 1.*t))*
        std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) + 24.*pow(M2,3)*(M2 - 4.*pow(Mt,2) + 2.*t)*std::log(pow(Mt,2)/(M2 - 1.*pow(Mt,2) + t))*
        std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) + 24.*pow(M2,2)*(3.*pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + M2*(-4.*pow(Mt,2) + 6.*t))*
        std::log(un*(M2 - 1.*pow(Mt,2) + t)/pow(Mt,2))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) + 
       84.*pow(M2,2)*(pow(M2,2) + 4.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*std::log(1. - (1.*t)/pow(Mt,2))*std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)) - 
       3.*pow(M2,2)*(31.*pow(M2,2) + 72.*pow(pow(Mt,2) - 1.*t,2) + M2*(-8.*pow(Mt,2) + 62.*t))*
        pow(std::log(un*(0.07957747154594767*pow(muR,2))/pow(Mt,2)),2) + 
       24.*pow(M2,2)*(pow(M2,2) + 2.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*std::log(pow(muR,2)/pow(Mt,2)) - 
       144.*pow(M2,2)*(pow(M2,2) + 2.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*std::log(un*(12.566370614359172*pow(muR,2))/pow(Mt,2)) - 
       (12.*M2*(M2 - 2.*pow(Mt,2))*(pow(M2,2) + 2.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)*
           gsl_sf_dilog(pow(-1. + sqrt(1. - (4.*pow(Mt,2))/M2),2)/pow(1. + sqrt(1. - (4.*pow(Mt,2))/M2),2)))/sqrt(1. - (4.*pow(Mt,2))/M2) - 
       (24.*M2*(M2 - 2.*pow(Mt,2))*(pow(M2,2) + 2.*pow(pow(Mt,2) - 1.*t,2) + 2.*M2*t)* gsl_sf_dilog(2./(1. + sqrt(1. - (4.*pow(Mt,2))/M2))))/
        sqrt(1. - (4.*pow(Mt,2))/M2) + (12.*pow(M2,2)*(13.*pow(M2,4) + pow(M2,3)*(-92.*pow(Mt,2) + 26.*t) + 
            4.*pow(M2,2)*(55.*pow(Mt,4) - 51.*pow(Mt,2)*t + 9.*pow(t,2)) - 4.*M2*(67.*pow(Mt,6) - 242.*pow(Mt,4)*t + 63.*pow(Mt,2)*pow(t,2)) + 
            32.*(17.*pow(Mt,8) - 44.*pow(Mt,6)*t + 27.*pow(Mt,4)*pow(t,2)))*
					gsl_sf_dilog((-1. + sqrt(1. - (4.*pow(Mt,2))/M2))/(1. + sqrt(1. - (4.*pow(Mt,2))/M2))))/(pow(M2 - 4.*pow(Mt,2),2)*Sqrt(1. - (4.*pow(Mt,2))/M2)) + 
       84.*pow(M2,3)*(M2 + 2.*t)* gsl_sf_dilog((-1.*t)/(pow(Mt,2) - 1.*t)) - 
       24.*pow(M2,3)*(M2 - 4.*pow(Mt,2) + 2.*t)* gsl_sf_dilog((M2 - 2.*pow(Mt,2) + t)/(M2 - 1.*pow(Mt,2) + t))))/pow(M2,4);
	}
      
    }

  else if (chan==1)
    {
      
    }
  
  res=2.*std::real(H[m][n])*pow(alpha_s,2)/dR/dR/2/M2;
  
  return res;
}
