import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import csv

HistNLO=np.loadtxt('../tt_VaryMG/MGfO/data/MinvNLO_ttx1_dXsec.txt',delimiter=" ")
HistLO=np.loadtxt('../tt_VaryMG/MGfO/data/MinvLO_ttx1_dXsec.txt',delimiter=" ")

M2hist=(np.transpose(HistLO)[0]+np.transpose(HistLO)[1])**2/4
NLOhist=-np.transpose(HistNLO)[2]/(np.transpose(HistNLO)[0]**2-np.transpose(HistNLO)[1]**2)
LOhist=-np.transpose(HistLO)[2]/(np.transpose(HistLO)[0]**2-np.transpose(HistLO)[1]**2)
Ratiohist=NLOhist/LOhist        


Megg=np.loadtxt("SubProcesses/P2_gg_ttx/Results/ExpanddM.txt",delimiter=",")
Borngg=np.loadtxt("SubProcesses/P2_gg_ttx/Results/BorndM.txt",delimiter=",")
Meqq=np.loadtxt("SubProcesses/P0_uux_ttx/Results/ExpanddM.txt",delimiter=",")
Bornqq=np.loadtxt("SubProcesses/P0_uux_ttx/Results/BorndM.txt",delimiter=",")


M2=np.transpose(Borngg)[-1]

M2=np.sqrt(M2)
M2hist=np.sqrt(M2hist)

#print(M2,M2hist)

sigmaBorndm=np.transpose(Borngg)[3]+np.transpose(Bornqq)[3]
RatioExpdm=np.transpose(Megg)[3]+np.transpose(Meqq)[3]
RatioExpdm/=sigmaBorndm 
RBorn=LOhist/sigmaBorndm

Rhist,= plt.plot(M2hist,Ratiohist,label=r'MG5@NLO')
R3dm, = plt.plot(M2,RatioExpdm, label=r'Resummed Exp@NLO')
#RBdm, = plt.plot(M2,RBorn, label=r'RBorn')
#plt.xscale("log")
plt.grid()
plt.xlabel(r'$M$ in ${\rm GeV}$')
plt.ylabel(r'$d\sigma/dM^2 $/LO')
#plt.ylim(2,10)
#plt.title(r'gg Ratio to LO, ansatz at $\sqrt{s}=1.5$TeV')
#plt.title(r'qqb Ratio to LO, ansatz at $\sqrt{s}=1.5$TeV')
plt.title(r'Total ratio to LO, ansatz at $\sqrt{s}=1.5$TeV')
plt.legend(handles=[R3dm,Rhist])
plt.savefig('XSecRatiodm.png')
plt.legend([])
plt.close()
