Implementation of resummation is checked against Mathematica for different values of muR and muF (one phase space point -> M, Nb, cos(theta)) for qqb channel (no Jacobian of Mellin PDF)
Expansion and difference checked for few values of (muR,muF).

7 point scale variation of M_tt differential cross section.
2 integrators are available: Vegas (MC) or Ronberg (trapezoidal method), the choice is to be made in the main.cpp file for each computation.
The code supports computation of pp > ttb at Born level, NLL, NLL_SV and NLL-NLL_SV. To perform the matching, the later should be added to fixed order (NLO).
The optimal center scale choice seems to be muR=muF=M_tt/2.
Make sure that the parameters in inc/Constants.h are coherent with the param_card.dat in the ../../Cards/ folder.
If anything is modified from the MadLoop side, make -f makefile_loop to update the fortran routine.
To modify something in inc/, make clean and then make (make alone is not enough).
Choice of PDF is done in Runpar.sh (enter LHAPDF ID, if it's supported).
