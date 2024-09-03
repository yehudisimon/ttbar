// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <iostream>   // In/Out streams                 //
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <complex>
// -------- Classes ----------------------------------- //
#include "process.h"      // Process definition         //
#include "messages.h"     // Message services           //
#include "LHAPDF/LHAPDF.h"
#include "PDFfitAN.h"
using namespace std;
using namespace LHAPDF;


// -------- Functions --------------------------------- //
void Integrate(double&,double&,double&,int&,double&);  // Integration
void pdfFit(double&, double (*)[8], double, double, double); 
// ---------------------------------------------------- //
extern double MZ, Mt, muF, muR, A1min, alpha_s, xmin, xmax;
extern double A[8][8];
extern int PDFID, chan, dim;
extern string namePDF;

extern "C" {
  void setmadlooppath_(char*);
  void setpara_(char*);
  void update_mur_(double&);
};


// ************************************************************************* //
//  Main code                                                                //
// ************************************************************************* //
void ConvertToFortran(char* fstring, std::size_t fstring_len, const char* cstring)
{
  std::size_t inlen = strlen(cstring);
  std::size_t cpylen = std::min(inlen, fstring_len);

  std::copy(cstring, cstring + cpylen, fstring);
  std::fill(fstring + cpylen, fstring + fstring_len, ' ');
}


int main(int argc, char* argv[])
{
  // Checking the arguments
  if(argc!=4)
  {
    std::cout << "This function requires 3 arguments\n"
	      << "  - PDF ID \n"
	      << "Current core and total number of cores \n";
    exit(0);
  }	

  
  const clock_t begin_time = clock();
  
  std::ofstream BorndM;
  std::ofstream NLLdM;
  
  const int ntot=31;
   
  double s0=7000.;

  double Mi=2.*Mt;
  double Mf=500.;

  int Ncores, CurrentCore;
  std::istringstream ncores(argv[3]); 
  std::istringstream currentCore(argv[2]);
  ncores >> Ncores;
  currentCore >> CurrentCore;

  int nstep[Ncores]={3,3,3,3,3,3,3,3,3,3};

  std::string name1("param_card.dat");
  std::string name2("../MadLoop5_resources");
  const size_t length = name1.length()+10;
  const size_t length2 = name2.length()+1;
  char cpp_string[length], stri[length];
  char cpp_string2[length2], stri2[length2];
  
  strcpy(cpp_string, name1.c_str());
  strcpy(cpp_string2, name2.c_str());
  ConvertToFortran(stri, sizeof stri, cpp_string);
  ConvertToFortran(stri2, sizeof stri2, cpp_string2);

  setpara_(stri);
  std::cout << " MadLoop parameter passed " << std::endl;
		  
		    
  const int nsteps=nstep[CurrentCore-1];
  int initstep=0;
  for(int k=1; k<CurrentCore; k++)
    {
      initstep+=nstep[k-1];
    }
  
  double Mic= Mi + (Mf-Mi)*double(initstep)/double(ntot-1)+(Mf-Mi)*0.5/double(ntot-1);
  double sc=s0, ML=Mic;
  
  const std::string & pathstr="/home/yehudi/LHAPDF6/share/LHAPDF";

  setPaths(pathstr);

  int melo=0;
  string name, tag;
  const PDF* F;  
  std::istringstream pdf(argv[1]);
  pdf >> PDFID;
  if(PDFID==14400) namePDF = "CT18NLO";
  else if(PDFID==14000) namePDF = "CT8NNLO";
  else if (PDFID==21100) namePDF = "MSTW2008nlo68cl";
  else if (PDFID==21200) namePDF = "MSTW2008nnlo68cl";
  F=mkPDF(namePDF,0);

  if(chan==0) tag="qqb";
  else if(chan==1) tag="gg";
  BorndM.open("BorndM_"+tag+".txt");
  NLLdM.open("ResumdM_"+tag+".txt");
  
  for(int ic=0; ic < nsteps; ic++) 
    {
      double M2=ML*ML;
      double s=sc*sc;

      for(int mu1=-1; mu1< 2; mu1++)
	{
	  muF=ML*pow(2.,mu1);
	  xmax=max(0.8,muF/sc*10.);
	  xmin=min(muF/sc/10.,0.01);

	  //Doing the fit, carefull: only valid for CT18NLO (hard coded in src/pdfN.cpp)
	  pdfFit(A1min,A,-1.6,-1.6,-1.6); //Mandatory for each muF=M if dynamic scale
	  
	  for(int mu2=-1; mu2< 2; mu2++)
	    {
	      if(mu1*mu2>=0)
		{
		  muR=ML*pow(2.,mu2);
		  update_mur_(muR);
		  alpha_s=F->alphasQ(muR);
		  
		  int mel=0;
		  double res=0., err=0.;
		  Integrate(res,err,s,mel,M2);
		  std::cout << "Born done " << std::endl;
		  // Writing in output file
		  BorndM << res << ",";
		  
		  mel=1;
		  double resD=0., errD=0.;
		  Integrate(resD,errD,s,mel,M2);
		  std::cout << "Diff done " << std::endl;
		  // Writing in output file
		  NLLdM << resD << ",";
		}
	    }
	}


      NLLdM << M2 << std::endl;
      BorndM << M2 << std::endl;

      std::cout << " Step n ° " << ic+1 << " over " << nsteps << std::endl;
      std::cout << " Duration : " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

      ML+=(Mf-Mi)/double(ntot-1);
    }

  BorndM.close();
  NLLdM.close();
  delete F;
  // End of the program
  return 0;
}
