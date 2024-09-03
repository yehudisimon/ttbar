// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <iostream>   // In/Out streams                 //
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
// -------- Classes ----------------------------------- //
#include "process.h"      // Process definition         //
#include "messages.h"     // Message services           //                                                           
//#include "Constants.h"
#include "PDFfitAN.h"
using namespace std;
using namespace LHAPDF;

// -------- Functions --------------------------------- //
void Integrate(double&,double&,double&,int&,double&);  // Integration
void pdfFit(double&, double (*)[8], double, double, double); 
// ---------------------------------------------------- //
extern double MZ, Mt, muF, muR, A1min, alpha_s, xmin, xmax;
extern double A[8][8];
// ************************************************************************* //
//  Main code                                                                //
// ************************************************************************* //
int main(int argc, char* argv[])
{
  // Checking the arguments
  if(argc!=4)
  {
    std::cout << "This function requires 3 arguments\n"
	      << "  - PDF ID (1=CT18 NLO, 2=CT14LN, 3=CT14LL, 4=CTEQ6M)\n"
	      << "Initial point and final point \n";
    exit(0);
  }	
  

  const clock_t begin_time = clock();
  
  std::ofstream TestdM;
  std::ofstream ExpanddM;
  std::ofstream ResumdM;
  
  const int ntot=21;
  //const int ntot=51;
   
  double s0=14000.;
  double Mi=2.*Mt;
  double Mf=5000.;

  int Ncores, CurrentCore;
  std::istringstream ncores(argv[3]); 
  std::istringstream currentCore(argv[2]);
  ncores >> Ncores;
  currentCore >> CurrentCore;

  //int nstep[Ncores]={5,5,5,5,5,5,5,5,5,3,2};
  //int nstep[Ncores]={5,5,5,5,5,5,5,5,5,5};
  int nstep[Ncores]={2,2,2,2,2,2,2,2,2,2};
  //int nstep[Ncores]={10,10,10,10,10,10,10,10,10,10};
  //int nstep[Ncores]={20,20,20,20,5,5,4,2,2,2};
  
  const int nsteps=nstep[CurrentCore-1];
  int initstep=0;
  for(int k=1; k<CurrentCore; k++)
    {
      initstep+=nstep[k-1];
    }
  
  double Mic= Mi + (Mf-Mi)*double(initstep)/double(ntot-1)+(Mf-Mi)*0.5/double(ntot-1);
  //double Mic= Mi*pow(Mf/Mi,double(initstep)/double(ntot-1));
  
  double sc=s0, M=Mic;
  
  const std::string & pathstr="/home/yehudi/LHAPDF6/share/LHAPDF";

  setPaths(pathstr);
  int melo=0;
  Process *init = new Process(argv[1], s0, melo, M);

  init->Init(); //PDF initialisation
  
  ExpanddM.open("ExpanddM.txt");
  TestdM.open("TestdM.txt");
  ResumdM.open("ResumdM.txt");
  
  for(int ic=0; ic < nsteps; ic++) 
    {
      double M2=M*M;
      double s=sc*sc;
      xmax=max(0.8,M/sc*10.);
      xmin=M/sc/10.;
      
      for(int mu1=-1; mu1< 2; mu1++)
	{
	  muF=M*pow(2.,mu1);
	  //Doing the fit, carefull: only valid for CT18NLO (hard coded in src/pdfN.cpp)
	  pdfFit(A1min,A,-1.6,-1.6,-1.6); //Mandatory for each muF=M if dynamic scale
	  
	  for(int mu2=-1; mu2< 2; mu2++)
	    {
	      if(mu1*mu2>=0)
		{
		  muR=M*pow(2.,mu2);
		  alpha_s=init->F->alphasQ(muR);
	
		  int mel=2;
		  double res0=0., err0=0.;
		  Integrate(res0,err0,s,mel,M2);
		  std::cout << "Resummed done " << std::endl;
		  // Writing in output file
		  ResumdM << res0 << ",";
      
		  mel=0;
		  double res=0., err=0.;
		  Integrate(res,err,s,mel,M2);
		  std::cout << "Born done " << std::endl;
		  // Writing in output file
		  TestdM << res << ",";
      
		  mel=3;
		  double resfit=0., errfit=0.;
		  Integrate(resfit,errfit,s,mel,M2);
		  std::cout << "Expanded done " << std::endl;
		  // Writing in output file
		  ExpanddM << resfit << ",";
		}
	    }
	}

      ResumdM << M2 << std::endl;
      TestdM << M2 << std::endl;
      ExpanddM << M2 << std::endl;
      std::cout << " Step n ° " << ic+1 << " over " << nsteps << std::endl;
      std::cout << " Duration : " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
      //M*=pow(Mf/Mi,1./(ntot-1));
      M+=(Mf-Mi)/double(ntot-1);
    }

  TestdM.close();
  ExpanddM.close();
  ResumdM.close();
  // End of the program
  return 0;
}
