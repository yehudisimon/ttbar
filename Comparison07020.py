import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc
import pandas as pd

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams.update({'font.size': 14})

indixMG=[10,11,12,13,14,16,18]
indixMe=[3,4,2,5,6,1,0]
#indixMG=[10,11,12] #muR = 0,1,-1 & muF=0
#indixMe=[1,2,0]
#indixMe=[3,4,2]
#indixMG=[10,13,16] #muF = 0,1,-1 & muR=0
#indixMe=[3,5,1]

rangeL=np.array([345,400,470,550,650,800,1100])
rangeR=np.array([400,470,550,650,800,1100,1600])
#rangeL=np.arange(2*173.2,4000+2*173.2,200)
#rangeR=np.arange(2*173.2+200,4000+2*173.2+200,200)
print(rangeL,rangeR)
errange=(rangeR-rangeL)/2
rangeM=(rangeR+rangeL)/2

labels=["NLO \; MG"]+["NLO \oplus NLL"]
types=["Diff"]
#types=["Born","Resum"]
color=['g','k','g','magenta','orange','m','b','c']

# n1=[0]
# n2=[199]
# nm=[10]

# n1MG=[0]
# n2MG=[199]
# nmMG=[10]

n1=[151,91,61,41,25,11,0]
n2=[250,150,90,60,40,24,10]
nm=[100,60,30,20,16,14,11]

n1MG=[151,91,61,41,25,11,0]
n2MG=[250,150,90,60,40,24,10]
nmMG=[100,60,30,20,16,14,11]


def combine(xL,xR,init,end):
    if(init==end):
        return (xL, xR)
    newxL=np.zeros(len(xL)-end+init)
    newxR=np.zeros(len(xR)-end+init)
    c=0
    for k in range(len(newxR)):
        if(k<init):
            newxL[k]=xL[k]
            newxR[k]=xR[k]
        elif(c==0):
            newxL[k]=xL[k]
            newxR[k]=xR[k+end-init]
            c=1
        else:
            newxL[k]=xL[k+end-init]
            newxR[k]=xR[k+end-init]
    return (newxL, newxR)

def mergex(xL,xR,init,end,n):
    if((end-init+1)%n !=0):
        print("Unsuitable bin range for this separation: rest = ",(end-init+1)%n)
        exit()
    newxR=np.copy(xR)
    newxL=np.copy(xL)
    for i in range(int((end-init+1)/n)):
        newxL, newxR = combine(newxL,newxR,init+i,init+i+n-1)
    return (newxL, newxR)

def sumhist(hist,init,end):
    if(init==end):
        return hist
    newhist=np.zeros(len(hist)-end+init)
    c=0
    for k in range(len(newhist)):
        if(k<init):
            newhist[k]=hist[k]
        elif(c==0):
            newhist[k]=sum([hist[p] for p in range(init, end+1)])
            c=1
        else:
            newhist[k]=hist[k+end-init]
    return newhist


def merge(hist,init,end,n):
    if((end-init+1)%n !=0):
        print("Unsuitable bin range for this separation: rest = ",(end-init+1)%n)
        exit()
    newhist=np.copy(hist)
    for i in range(int((end-init+1)/n)):
        newhist=sumhist(newhist,init+i,init+i+n-1)
    return newhist


def surmerge(hist,xL,xR,init,end,n):
    newhist=np.copy(hist)
    newxL=np.copy(xL)
    newxR=np.copy(xR)
    for i in range(len(n)):
        newhist=merge(newhist,init[i],end[i],n[i])
        newxL, newxR= mergex(newxL,newxR,init[i],end[i],n[i])
    return (newhist, newxL, newxR)


def loadgg(order,res):
    #Mycode=np.loadtxt("SubProcesses/P2_gg_ttx/Results/"+order+"dM.txt",delimiter=",")
    #Mycode=np.loadtxt("SubProcesses/P2_gg_ttx/Results/Bkp13TeV.txt",delimiter=",")
    Mycode=np.loadtxt("SubProcesses/P2_gg_ttx/Results/Bkp8TeV.txt",delimiter=",")

    M=np.transpose(Mycode)[-1]
    M=np.sqrt(M)
    Mb=np.copy(M)
    a=1
    ML=M-(M[1]-M[0])/2
    MR=M+(M[1]-M[0])/2
    dres=[np.transpose(Mycode)[p]*2*M*(MR-ML) for p in range(len(np.transpose(Mycode))-1)]
    
    dres[0],ML,MR=surmerge(dres[0],ML,MR,n1,n2,nm)
    for k in range(1,len(dres)):
        dres[k],a,a=surmerge(dres[k],Mb,Mb,n1,n2,nm)

    res=np.zeros((len(dres),len(dres[0])))
    for k in range(len(dres)):
        res[k]=dres[k]/(MR-ML) # mu_R=mu_F=M -> k=3

    return(ML,MR,res)
        


def loadqqb(order,res):
    #Mycode=np.loadtxt("SubProcesses/P0_uux_ttx/Results/"+order+"dM.txt",delimiter=",")
    #Mycode=np.loadtxt("SubProcesses/P0_uux_ttx/Results/Bkp13TeVDiff.txt",delimiter=",")
    Mycode=np.loadtxt("SubProcesses/P0_uux_ttx/Results/Bkp8TeVDiff.txt",delimiter=",")
    M=np.transpose(Mycode)[-1]
    M=np.sqrt(M)
    Mb=np.copy(M)
    a=1
    ML=M-(M[1]-M[0])/2
    MR=M+(M[1]-M[0])/2
    
    dres=[np.transpose(Mycode)[p]*2*M*(MR-ML) for p in range(len(np.transpose(Mycode))-1)]
    dres[0],ML,MR=surmerge(dres[0],ML,MR,n1,n2,nm)
    for k in range(1,len(dres)):
        dres[k],a,a=surmerge(dres[k],Mb,Mb,n1,n2,nm)

    res=np.zeros((len(dres),len(dres[0])))
    for k in range(len(dres)):
        res[k]=dres[k]/(MR-ML) # mu_R=mu_F=M -> k=3

    return(ML,MR,res)
        
def loadMGNLO(res,resExt,ntot):
    
    for i in range(1,ntot+1):
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/MinvNLO_ttx'+str(i)+'_dXsec.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/BkpTotNLO13TeV.txt',delimiter=" ")
        MadG=np.loadtxt('../tt_VaryMG/MGfO/data/BkpTotNLO8TeV.txt',delimiter=" ")

        resc=np.array([np.transpose(MadG)[p] for p in indixMG])
        MG=np.transpose(MadG)[0]
        MR=np.transpose(MadG)[1]

        if(i==1):
            MGf=MG
            MRf=MR
            resf=resc
        else:
            MGf=np.concatenate((MGf,MG),axis=None)
            MRf=np.concatenate((MRf,MR),axis=None)
            resf=np.concatenate((resf,resc),axis=1)
        
    resExt=np.zeros((2,len(MGf)))
    for k in range(len(MGf)):
        resExt[0][k]=min(resf[p][k] for p in range(len(resf)))
        resExt[1][k]=max(resf[p][k] for p in range(len(resf)))
        
    res,MG,MR=surmerge(resf[0],MGf,MRf,n1MG,n2MG,nmMG)

    resExtf=np.zeros((2,len(MG)))
    resExtf[0],MG,MR=surmerge(resExt[0],MGf,MRf,n1MG,n2MG,nmMG)
    resExtf[1],MG,MR=surmerge(resExt[1],MGf,MRf,n1MG,n2MG,nmMG)

    res=res/(MR-MG)
    resExtf[0]=resExtf[0]/(MR-MG)
    resExtf[1]=resExtf[1]/(MR-MG)
    
    resMu=np.zeros((len(resf),len(MG)))
    for p in range(len(resf)):
        resMu[p],MG,MR=surmerge(resf[p],MGf,MRf,n1MG,n2MG,nmMG)
        resMu[p]/=(MR-MG)

    return(MG,MR,res,resExtf,resMu)



def BuildPlot():
    figure, axis=plt.subplots(2,gridspec_kw={'height_ratios': [1.67, 1]})    

    NLO07= np.loadtxt('NLO.csv',delimiter=';')
    #NLO07= np.loadtxt('NLO13.csv',delimiter=';')
    NLO0 = np.copy(np.transpose(NLO07)[1])
    # NLOm=np.copy(NLO0[::2])
    # NLOp=np.copy(NLO0[1::2])
    NLOp=np.copy(NLO0[::2])
    NLOm=np.copy(NLO0[1::2])

    Matched07= np.loadtxt('Matched.csv',delimiter=';')
    #Matched07= np.loadtxt('Matched13.csv',delimiter=';')
    Matched = np.copy(np.transpose(Matched07)[1])
    # Matchedm=np.copy(Matched[::2])
    # Matchedp=np.copy(Matched[1::2])
    Matchedp=np.copy(Matched[::2])
    Matchedm=np.copy(Matched[1::2])

    Atl = np.loadtxt('Atlas8TeV.csv',delimiter=';')
    Atlas = np.copy(np.transpose(Atl)[1])
    AtlasM=np.copy(Atlas[:7])
    Atlasp=np.copy(Atlas[7:14])
    Atlasm=np.copy(Atlas[14:])
    Atlaserr=np.zeros((2,len(AtlasM)))
    Atlaserr[0]=np.copy(AtlasM-Atlasm)
    Atlaserr[1]=np.copy(Atlasp-AtlasM)

    print(AtlasM,Atlasp,Atlasm)
    #print(AtlasM,Atlaserr[0],Atlaserr[1])

    #print("range=",rangeL)
    #print("NLO=",NLOm)
    
    axis[0].fill_between(rangeL,NLOm,NLOp,alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+4],label=r'${\rm \mathbf{NLO \, PRL}}$ ',step='post')
    axis[0].fill_between(rangeL,Matchedm,Matchedp,alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+5],label=r'${\rm \mathbf{NLO + NNLLp \, PRL}}$ ',step='post')
    axis[0].fill_between(rangeR[-2:],NLOm[-2:],NLOp[-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+4],step='pre')
    axis[0].fill_between(rangeR[-2:],Matchedm[-2:],Matchedp[-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+5],step='pre')

    # axis[0].scatter(rangeM,NLO0,alpha=0.5,linewidth=1,antialiased=True,color=color[len(types)+4],label=r'${\rm \mathbf{NLO \, PRL}}$ ',marker='^')
    # axis[0].scatter(rangeM,Matched,alpha=0.5,linewidth=1,antialiased=True,color=color[len(types)+5],label=r'${\rm \mathbf{NLO + NNLLp \, PRL}}$ ',marker='+')
    
    axis[0].errorbar(rangeM,AtlasM,yerr=Atlaserr,xerr=errange,fmt=',r',ecolor='r',elinewidth=0.5,capsize=1.5,barsabove=True,label=r'${\rm \mathbf{Atlas}}$')
    
    ntot=1
    resMGNLO=resMGExtNLO=1
    resMu=1
    aNLO=bNLO=1
    
    aNLO, bNLO, resMGNLO, resMGExtNLO, resMu = loadMGNLO(resMGNLO,resMGExtNLO,ntot) 
    
    axis[0].step(aNLO,resMGNLO, c=color[0], linewidth=0.5,where='post')
    axis[0].fill_between(aNLO,resMGExtNLO[0],resMGExtNLO[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[0],label=r'${\rm \mathbf{'+labels[0]+'}}$ ',step='post')
    axis[0].step(aNLO,resMGExtNLO[0], c=color[0], linewidth=0.5, alpha=0.3,where='post')
    axis[0].step(aNLO,resMGExtNLO[1], c=color[0], linewidth=0.5, alpha=0.3,where='post')

    axis[0].step(bNLO[-2:],resMGNLO[-2:], c=color[0], linewidth=0.5,where='pre')
    axis[0].fill_between(bNLO[-2:],resMGExtNLO[0][-2:],resMGExtNLO[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[0],step='pre')
    axis[0].step(bNLO[-2:],resMGExtNLO[0][-2:], c=color[0], linewidth=0.5, alpha=0.3,where='pre')
    axis[0].step(bNLO[-2:],resMGExtNLO[1][-2:], c=color[0], linewidth=0.5, alpha=0.3,where='pre')
   
    i=0    
    resgg=resqqb=1
    
    ML, MR, resgg = loadgg("Diff",resgg)
    ML, MR, resqqb = loadqqb("Diff",resqqb)

    resTot=resgg+resqqb
    combine=[resMu[p]+resTot[indixMe[p]] for p in range(len(indixMe))]
    resMatched=np.copy(combine[0])
    
    combineExt=np.zeros((2,len(combine[0])))
    for k in range(len(combine[0])):
        combineExt[0][k]=min(combine[p][k] for p in range(len(indixMe)))
        combineExt[1][k]=max(combine[p][k] for p in range(len(indixMe)))
        
    axis[0].step(ML,combine[0], c=color[-1], linewidth=0.5,where='post')
    axis[0].fill_between(ML,combineExt[0],combineExt[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[-1],label=r'${\rm \mathbf{NLO \oplus NLL}}$ ',step='post')
    axis[0].step(ML,combineExt[0], c=color[-1], linewidth=0.5, alpha=0.3,where='post')
    axis[0].step(ML,combineExt[1], c=color[-1], linewidth=0.5, alpha=0.3,where='post')

    axis[0].step(MR[-2:],combine[0][-2:], c=color[-1], linewidth=0.5,where='pre')
    axis[0].fill_between(MR[-2:],combineExt[0][-2:],combineExt[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[-1],step='pre')
    axis[0].step(MR[-2:],combineExt[0][-2:], c=color[-1], linewidth=0.5, alpha=0.3,where='pre')
    axis[0].step(MR[-2:],combineExt[1][-2:], c=color[-1], linewidth=0.5, alpha=0.3,where='pre')


    for k in range(len(combineExt)):
        combineExt[k]/=resMatched
       
    for p in range(len(combine)):
        combine[p]/=resMatched

    for k in range(len(resMGExtNLO)):
        resMGExtNLO[k]/=resMatched
        
    resMGNLO/=resMatched    
    AtlasM/=resMatched
    Atlaserr[0]/=resMatched
    Atlaserr[1]/=resMatched
    NLOm/=resMatched
    Matchedm/=resMatched
    NLOp/=resMatched
    Matchedp/=resMatched
    
    axis[1].errorbar(rangeM,AtlasM,yerr=Atlaserr,xerr=errange,fmt=',r',ecolor='r',elinewidth=0.5,capsize=1.5,barsabove=True)
    #axis[1].scatter(rangeM,NLO0,alpha=0.5,linewidth=1,antialiased=True,color=color[len(types)+4],marker='^')
    #axis[1].scatter(rangeM,Matched,alpha=0.5,linewidth=1,antialiased=True,color=color[len(types)+5],marker='+')

    axis[1].fill_between(rangeL,NLOm,NLOp,alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+4],label=r'${\rm \mathbf{NLO \, PRL}}$ ',step='post')
    axis[1].fill_between(rangeL,Matchedm,Matchedp,alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+5],label=r'${\rm \mathbf{NLO + NNLLp \, PRL}}$ ',step='post')
    axis[1].fill_between(rangeR[-2:],NLOm[-2:],NLOp[-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+4],step='pre')
    axis[1].fill_between(rangeR[-2:],Matchedm[-2:],Matchedp[-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+5],step='pre')

    
    axis[1].step(aNLO,resMGNLO, c=color[i+len(types)+1], linewidth=0.5, where='post')
    axis[1].step(aNLO,resMGExtNLO[0], c=color[i+len(types)+1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].step(aNLO,resMGExtNLO[1], c=color[i+len(types)+1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].fill_between(aNLO,resMGExtNLO[0],resMGExtNLO[1],alpha=0.5,facecolor=color[i+len(types)+1],linewidth=1,antialiased=True,step='post')

    axis[1].step(bNLO[-2:],resMGNLO[-2:], c=color[i+len(types)+1], linewidth=0.5, where='pre')
    axis[1].step(bNLO[-2:],resMGExtNLO[0][-2:], c=color[i+len(types)+1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].step(bNLO[-2:],resMGExtNLO[1][-2:], c=color[i+len(types)+1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].fill_between(bNLO[-2:],resMGExtNLO[0][-2:],resMGExtNLO[1][-2:],alpha=0.5,facecolor=color[i+len(types)+1],linewidth=1,antialiased=True,step='pre')

    
    axis[1].step(ML,combine[0],label=r''+labels[-1]+' ', c=color[-1], linewidth=0.5,where='post')
    axis[1].step(ML,combineExt[0], c=color[-1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].step(ML,combineExt[1], c=color[-1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].fill_between(ML,combineExt[0],combineExt[1],alpha=0.5,facecolor=color[-1],linewidth=1,antialiased=True,step='post')

    axis[1].step(MR[-2:],combine[0][-2:], c=color[-1], linewidth=0.5,where='pre')
    axis[1].step(MR[-2:],combineExt[0][-2:], c=color[-1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].step(MR[-2:],combineExt[1][-2:], c=color[-1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].fill_between(MR[-2:],combineExt[0][-2:],combineExt[1][-2:],alpha=0.5,facecolor=color[-1],linewidth=1,antialiased=True,step='pre')


    axis[0].legend(loc='best', fontsize = 12)
    axis[0].set_yscale("log")
    axis[0].set_xticklabels([])

    #axis[1].set_ylim([0.2,1.4])
    #axis[0].set_ylim([0.001,10])

    #axis[1].set_xlim([2*173.2,4000+2*173.2])
    #axis[0].set_xlim([173.2*2,4000+2*173.2])

    axis[0].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='y',which='both',direction='in')
    axis[0].tick_params(axis='y',which='both',direction='in')
   
    axis[1].yaxis.set_minor_locator(MultipleLocator(0.1))
    axis[1].xaxis.set_minor_locator(MultipleLocator(100))
    axis[1].xaxis.set_major_locator(MultipleLocator(500))

    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1\
),numticks=12)
    axis[0].yaxis.set_minor_locator(locmin)
    axis[0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    
    axis[0].grid(linestyle='dashed',linewidth=0.5)
    axis[1].grid(linestyle='dashed',linewidth=0.5)
    
    
    axis[0].set_ylabel(r'$\dfrac{{\rm d} \sigma}{{\rm d} M}$ [pb/GeV]',fontsize=16)
    axis[1].set_ylabel(r'$\sigma/\sigma_{{\rm ref}}$',fontsize=16)
    axis[1].set_xlabel(r'${\rm M}$ [GeV]',fontsize=16)
        
    figure.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('Comparison07020.svg')
    plt.savefig('Comparison07020.png')
    plt.legend([])
    plt.close()
 
BuildPlot()

    
