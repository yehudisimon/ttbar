#ifndef PROCESS_H
#define PROCESS_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <string>    // String streams                  //
#include "LHAPDF/LHAPDF.h"
// -------- Processor variables ----------------------- //
using namespace LHAPDF;
#define cst const std::string                           //
// -------- Functions --------------------------------- //
extern "C" { void setct18_(char[40]); void setctq6_(int&); };                    //
// ---------------------------------------------------- //

// ************************************************************************* //
// Definition of the class Process                                           //
// ************************************************************************* //
class Process
{
public:
  Process()  { };
  //Process(cst&,cst&,cst&);
  Process(cst &, double&, int, double&);
  ~Process() { };
  
  // Kinematical quantities and scales
  double sh;               // Hadronic center of mass energy
  //int NLO;              // Cross section ID 0=LO, 1=NLO, 2=Resummed, 3=Expanded 
  int pdf;              // PDF identifier
  //int rep;                // number of relevant representation

  double m2; 
  int mel;
  const PDF* F;
  
  // Functions
  void Init();          // Init of PDFs, etc...
};

//Process *p;
#endif
