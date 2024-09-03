import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams.update({'font.size': 14})

indixMG=[10,11,12,13,14,16,18]
indixMe=[3,4,2,5,6,1,0]
#indixMG=[10,11,12] #muR = 0,1,-1 & muF=0
#indixMe=[1,2,0]
#indixMe=[3,4,2]
#indixMG=[10,13,16] #muF = 0,1,-1 & muR=0
#indixMe=[3,5,1]

labels=["LO"]+["NLL SCET"]+["LO \, MG","NLO \; MG"]+["NLO \oplus NLL"]
types=["Born","Resum"]
color=['b','g','k','r','orange','m','y','c']

n1=[0]
n2=[29]
nm=[1]

n1MG=[0]
n2MG=[29]
nmMG=[1]

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
    Mycode=np.loadtxt("SubProcesses/P2_gg_ttx/Results/"+order+"dM.txt",delimiter=",")
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
    Mycode=np.loadtxt("SubProcesses/P0_uux_ttx/Results/"+order+"dM.txt",delimiter=",")
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
        

def loadMGLO(res,resExt,ntot):
    
    for i in range(1,ntot+1):
        MadG=np.loadtxt('../tt_VaryMG/MGfO/data/MinvLO_ttx'+str(i)+'_dXsec.txt',delimiter=" ")

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

    return(MG,MR,res,resExtf)

def loadMGNLO(res,resExt,ntot):
    
    for i in range(1,ntot+1):
        MadG=np.loadtxt('../tt_VaryMG/MGfO/data/MinvNLO_ttx'+str(i)+'_dXsec.txt',delimiter=" ")

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
    figure, axis=plt.subplots(2,gridspec_kw={'height_ratios': [0.001,1]})    
    ntot=1
    resgg=1
    resqqb=1
    ML=1
    MR=1
    i=0
    resTot=1
    
    ML, MR, resgg = loadgg("Resum",resgg)
    ML, MR, resqqb = loadqqb("Resum",resqqb)
    
    #print((ML+MR)/2,(ML[30]+MR[30])/2,(ML[32]+MR[32])/2,(ML[33]+MR[33])/2)
    
    resTot=resgg+resqqb
    #resTot=resqqb
    #resTot=resgg
    
    resExt=np.zeros((2,len(resTot[0])))
    for k in range(len(resTot[0])):
        resExt[0][k]=min(resTot[p][k] for p in indixMe)
        resExt[1][k]=max(resTot[p][k] for p in indixMe)
        
    print((ML[10]+MR[10])/2,resTot[3][10]*1000,(resTot[3][10]-resExt[0][10])*1000,(resTot[3][10]-resExt[1][10])*1000)
    #print((ML[35]+MR[35])/2,resTot[3][35]*1000,(resTot[3][35]-resExt[0][35])*1000,(resTot[3][35]-resExt[1][35])*1000)
        
    #axis[0].step(ML,resTot[1], c=color[i+1], linewidth=0.5,where='post')
    axis[1].step(ML,resTot[3], c=color[i+1], linewidth=0.5,where='post')
    axis[1].fill_between(ML,resExt[0],resExt[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i+1],label=r'${\rm \mathbf{'+labels[i+1]+'}}$ ',step='post')
    axis[1].step(ML,resExt[0], c=color[i+1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].step(ML,resExt[1], c=color[i+1], linewidth=0.5, alpha=0.3,where='post')

    #axis[0].step(MR[-2:],resTot[1][-2:], c=color[i], linewidth=0.5,where='pre')
    axis[1].step(MR[-2:],resTot[3][-2:], c=color[i+1], linewidth=0.5,where='pre')
    axis[1].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i+1],step='pre')
    axis[1].step(MR[-2:],resExt[0][-2:], c=color[i+1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].step(MR[-2:],resExt[1][-2:], c=color[i+1], linewidth=0.5, alpha=0.3,where='pre')
        
    axis[1].legend(loc='best', fontsize = 12)
    #axis[0].set_yscale("log")
    #axis[0].set_xticklabels([])
    axis[1].set_ylim([0,0.9])
    
    
    axis[1].tick_params(axis='x',which='both',direction='in',bottom=True)
    #axis[0].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='y',which='both',direction='in')
   
    axis[1].yaxis.set_minor_locator(MultipleLocator(0.025))    
    axis[1].grid(linestyle='dashed',linewidth=0.5)
        
    axis[1].set_ylabel(r'$\dfrac{{\rm d} \sigma}{{\rm d} M}$ [pb/GeV]',fontsize=16)
    axis[1].set_xlabel(r'${\rm M}$ [GeV]',fontsize=16)
        
    figure.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('SCET.svg')
    plt.savefig('SCET.png')
    plt.legend([])
    plt.close()
 
BuildPlot()

    
