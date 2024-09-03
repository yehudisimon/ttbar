import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams.update({'font.size': 14})

# indixMG=[10,11,12,13,14,16,18]
# indixMe=[3,4,2,5,6,0,1]
indixMG=[10,11,12] #muR = 0,1,-1 & muF=0                                                           
#indixMe=[1,2,0]                                                                                   
indixMe=[3,4,2]

labels=["LO"]+["NLL"]+["NLL_1"]+["NLL_2"]
types=["Born","Resum","Expand","Diff"]
color=['b','g','k','r','orange','m','y','c']

n1=[0]
n2=[29]
nm=[1]

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
        

def BuildPlot():
    figure, axis=plt.subplots(2,gridspec_kw={'height_ratios': [1.67, 1]})    
    for i in range(len(types)):
        resgg=1
        resqqb=1
        resTot=1
        resExt=1
        resMG=1
        resMGExt=1
        resMGNLO=1
        resMGExtNLO=1
        resMu=1
        ML=1
        MR=1
        a=1
        b=1
        aNLO=1
        bNLO=1
        
        ML, MR, resgg = loadgg(types[i],resgg)
        #ML, MR, resqqb = loadqqb(types[i],resqqb)
        #resTot=resqqb+resgg
        #resTot=resqqb
        resTot=resgg
        
        if(types[i]=="Resum"):
            resNLL=np.copy(resTot)
        elif(types[i]=="Expand"):
            resExp=np.copy(resTot)
        
        if(types[i]=="Born"):
            resBorn=np.copy(resTot[3])    
            
        resExt=np.zeros((2,len(resTot[0])))
        for k in range(len(resTot[0])):
            resExt[0][k]=min(resTot[p][k] for p in indixMe)
            resExt[1][k]=max(resTot[p][k] for p in indixMe)
        

            
        axis[0].step(ML,resTot[3], c=color[i], linewidth=0.5,where='post')
        axis[0].fill_between(ML,resExt[0],resExt[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],label=r'${\rm \mathbf{'+labels[i]+'}}$ ',step='post')
        axis[0].step(ML,resExt[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
        axis[0].step(ML,resExt[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')

        axis[0].step(MR[-2:],resTot[3][-2:], c=color[i], linewidth=0.5,where='pre')
        axis[0].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],step='pre')
        axis[0].step(MR[-2:],resExt[0][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
        axis[0].step(MR[-2:],resExt[1][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
            
    
        for k in range(len(resExt)):
            resExt[k]/=resBorn
                        
        for p in range(len(resTot)):
            resTot[p]/=resBorn


        axis[1].step(ML,resTot[3],label=r''+labels[i]+' ', c=color[i], linewidth=0.5,where='post')
        # axis[1].step(ML,resExt[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
        # axis[1].step(ML,resExt[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')
        # axis[1].fill_between(ML,resExt[0],resExt[1],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='post')

        axis[1].step(MR[-2:],resTot[3][-2:], c=color[i], linewidth=0.5,where='pre')
        # axis[1].step(MR[-2:],resExt[0][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
        # axis[1].step(MR[-2:],resExt[1][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
        # axis[1].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='pre')


            
    axis[0].legend(loc='best', fontsize = 12)
    axis[0].set_yscale("log")
    axis[0].set_xticklabels([])
    
    #axis[1].set_ylim([0,10])
    #axis[1].set_ylim([0,5])
    
    axis[0].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='y',which='both',direction='in')
    axis[0].tick_params(axis='y',which='both',direction='in')
   
    
    axis[0].grid(linestyle='dashed',linewidth=0.5)
    axis[1].grid(linestyle='dashed',linewidth=0.5)
    
    
    axis[0].set_ylabel(r'$\dfrac{{\rm d} \sigma}{{\rm d} M}$ [pb/GeV]',fontsize=16)
    axis[1].set_ylabel(r'$\sigma/\sigma_{{\rm ref}}$',fontsize=16)
    axis[1].set_xlabel(r'${\rm M}$ [GeV]',fontsize=16)
        
    figure.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('ComparisonExp.svg')
    plt.savefig('ComparisonExp.png')
    plt.legend([])
    plt.close()
 
BuildPlot()

    
