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
double TotDiff(double*,double&,int&,double&);
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
	      << "  - PDF ID (1=CT18 NLO, 2=CT14LN, 3=CT14LL, 4=CTEQ6M)\n"
	      << "Initial point and final point \n";
    exit(0);
  }	

  
  const clock_t begin_time = clock();
  
  std::ofstream TestdM;
  // std::ofstream ExpanddM;
  // std::ofstream ResumdM;
  std::ofstream Diff;
  
  //const int ntot=21;
  const int ntot=51;
   
  double s0=14000.;//1400.;

  double Mi=2.*Mt;
  double Mf=5000.;

  int Ncores, CurrentCore;
  std::istringstream ncores(argv[3]); 
  std::istringstream currentCore(argv[2]);
  ncores >> Ncores;
  currentCore >> CurrentCore;

  //int nstep[Ncores]={1,1,1,1,1,1,1,1,1,1};
  //int nstep[Ncores]={2,2,2,2,2,2,2,2,2,2};
  int nstep[Ncores]={5,5,5,5,5,5,5,5,5,5};
  //int nstep[Ncores]={10,10,10,10,10,10,10,10,10,10};

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
  double sc=s0, M=Mic;
  
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

  double alpha_0=1e-7;

  alpha_s=F->alphasQ(M)*pow(alpha_0/F->alphasQ(M),(double(initstep)+.5)/double(ntot-1));
  M=4000.;
  
  if(chan==0) tag="qqb";
  else if(chan==1) tag="gg";

  double x[2]={0.3,0.3};
  
  // ExpanddM.open("ExpanddM_"+tag+".txt");
  TestdM.open("BorndM_"+tag+".txt");
  // ResumdM.open("ResumdM_"+tag+".txt");
  Diff.open("Diff_"+tag+".txt");
  
  for(int ic=0; ic < nsteps; ic++) 
    {
      double M2=M*M;
      double s=sc*sc;
      xmax=max(1.,M/sc*10.);
      xmin=min(M/sc/10.,0.01);

      //for(int mu1=-1; mu1< 2; mu1++)
      for(int mu1=0; mu1< 1; mu1++)
	{
	  muF=M*pow(2.,mu1);

	  //Doing the fit, carefull: only valid for CT18NLO (hard coded in src/pdfN.cpp)
	  pdfFit(A1min,A,-1.6,-1.6,-1.6); //Mandatory for each muF=M if dynamic scale

	  //for(int mu2=-1; mu2< 2; mu2++)
	  for(int mu2=1; mu2< 2; mu2++)
	    {
	      if(mu1*mu2>=0)
		{
		  muR=M*pow(2.,mu2);
		  update_mur_(muR);
		  //alpha_s=F->alphasQ(muR);

		  std::cout << " alpha_s main = " << alpha_s << std::endl;
		  
		  int mel=1;
		  double res0=0., err0=0.;
		  res0=TotDiff(x,sc,mel,M2);
		  std::cout << " Diff done " << std::endl;
		  // Writing in output file
		  Diff << res0 << ",";
      
		  // int mel=2;
		  // double res0=0., err0=0.;
		  // Integrate(res0,err0,s,mel,M2);
		  // std::cout << "Resummed done " << std::endl;
		  // // Writing in output file
		  // ResumdM << res0 << ",";
      
		  mel=0;
		  double res=0., err=0.;
		  //Integrate(res,err,s,mel,M2);
		  res=TotDiff(x,sc,mel,M2);
		  std::cout << "Born done " << std::endl;
		  // Writing in output file
		  TestdM << res << ",";
      
		  // mel=3;
		  // double resfit=0., errfit=0.;
		  // Integrate(resfit,errfit,s,mel,M2);
		  // std::cout << "Expanded done " << std::endl;
		  // // Writing in output file
		  // ExpanddM << resfit << ",";
		}
	    }
	}

      // ResumdM << M2 << std::endl;
      // TestdM << M2 << std::endl;
      // ExpanddM << M2 << std::endl;

      Diff << alpha_s << std::endl;
      TestdM << alpha_s << std::endl;
      
      std::cout << " Step n ° " << ic+1 << " over " << nsteps << std::endl;
      std::cout << " Duration : " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

      alpha_s*=pow(alpha_0/F->alphasQ(M),1./double(ntot-1));
      //M+=(Mf-Mi)/double(ntot-1);
    }

  TestdM.close();
  // ExpanddM.close();
  // ResumdM.close();
  Diff.close();
  
  // End of the program
  return 0;
}
