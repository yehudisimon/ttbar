import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams.update({'font.size': 14})

labels=["LO","LO \; MG","NLO \; MG"]
color=['b','g','k','r','orange','m','y','c']

n1=[22]
n2=[49]
nm=[2]


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




def loadataMycode(order,res,resExt):
    Mycode=np.loadtxt("MyCodeResum/Results/"+order+"dM.txt",delimiter=",")
    M=np.transpose(Mycode)[1]
    M=np.sqrt(M)
    MG=M-(M[1]-M[0])/2
    MR=M+(M[1]-M[0])/2
    dres=np.transpose(Mycode)[0]*2*M*(MR-MG)
    dres,MG,MR=surmerge(dres,MG,MR,n1,n2,nm)
    res=dres/(MR-MG)
    # resExt=np.zeros((2,len(s)))
    # tmp=np.zeros((7,len(s)))
    # kid=0
    # for muF in range(-1,2):
    #     for muR in range(-1,2):
    #         if(muR*muF>=0):
    #             Mycode=np.loadtxt("MyCodeResum/Results/"+order+"Z_muF"+str(muF)+"_muR"+str(muR)+".txt",delimiter=",")
    #             tmp[kid]=np.transpose(Mycode)[0]
    #             kid+=1
                   
    # for k in range(len(s)):
    #     resExt[0][k]=min(tmp[p][k] for p in range(len(tmp)))
    #     resExt[1][k]=max(tmp[p][k] for p in range(len(tmp)))
    return(M,res,resExt)
        

def loadataMG(order,res,resExt):
    MadG=np.loadtxt('MGfO/data/Minv'+order+'_ttx1_dXsec.txt',delimiter=" ")
    MG=np.transpose(MadG)[0]
    MR=np.transpose(MadG)[1]
    indix=np.array([10,11,12,13,14,16,18])
    resExt=np.zeros((2,len(MG)))
    for k in range(len(MG)):
        resExt[0][k]=min(np.transpose(MadG)[p][k] for p in indix)
        resExt[1][k]=max(np.transpose(MadG)[p][k] for p in indix)

    res,MG,MR=surmerge(np.transpose(MadG)[2],np.transpose(MadG)[0],np.transpose(MadG)[1],n1,n2,nm)

    resExtf=np.zeros((2,len(MG)))
    resExtf[0],MG,MR=surmerge(resExt[0],np.transpose(MadG)[0],np.transpose(MadG)[1],n1,n2,nm)
    resExtf[1],MG,MR=surmerge(resExt[1],np.transpose(MadG)[0],np.transpose(MadG)[1],n1,n2,nm)

    res=res/(MR-MG)
    resExtf[0]=resExtf[0]/(MR-MG)
    resExtf[1]=resExtf[1]/(MR-MG)

    return(MG,MR,res,resExtf)
        

def BuildPlot():
    res=1
    resExt=1
    M=1
    MG=1
    MR=1
    resS=1
    resExtS=1
    i=0

    figure, axis=plt.subplots(2,gridspec_kw={'height_ratios': [1.67, 1]})

    M, resS, resExtS = loadataMycode("Born",resS,resExtS)
       
    resLOS=np.copy(resS)
        
    #axis[0].step(M,resS, c=color[i], linewidth=0.5,where='post')
    #axis[0].fill_between(M,resExtS[0],resExtS[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],label=r'${\rm \mathbf{'+labels[i]+'}}$ ',step='post')
    #axis[0].step(M,resExtS[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    #axis[0].step(M,resExtS[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')

        
    # for a in range(len(resExtS)):
    #     resExtS[a]/=resLOS
    # resS/=resLOS

    Dummy=np.ones(len(resS))
    #axis[1].step(M,Dummy,  label=r''+labels[i]+' ', c=color[i], linewidth=0.5,where='post')
    #axis[1].step(M,resExtS[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    #axis[1].step(M,resExtS[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    #axis[1].fill_between(M,resExtS[0],resExtS[1],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='post')
       
    i+=1
        
    MG, MR, res, resExt= loadataMG("LO",res,resExt)
    #resLOMG=np.copy(res)
    axis[0].step(MG,resS, c=color[i-1], linewidth=0.5,where='post')
    axis[1].step(MG,Dummy,  label=r''+labels[i-1]+' ', c=color[i-1], linewidth=0.5,where='post')

    axis[0].step(MR[-2:],resS[-2:], c=color[i-1], linewidth=0.5,where='pre')
    axis[1].step(MR[-2:],Dummy[-2:], c=color[i-1], linewidth=0.5,where='pre')

    axis[0].fill_between(MG,resExt[0],resExt[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],label=r'${\rm \mathbf{'+labels[i]+'}}$ ',step='post')
    axis[0].step(MG,res, c=color[i], linewidth=0.5,where='post')
    axis[0].step(MG,resExt[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    axis[0].step(MG,resExt[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')

    axis[0].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],step='pre')
    axis[0].step(MR[-2:],res[-2:], c=color[i], linewidth=0.5,where='pre')
    axis[0].step(MR[-2:],resExt[0][-2:], c=color[i][-2:], linewidth=0.5, alpha=0.3,where='pre')
    axis[0].step(MR[-2:],resExt[1][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')

    
    for a in range(len(resExt)):
        resExt[a]/=resLOS#MG
    res/=resLOS#MG

    axis[1].step(MG,res, c=color[i], linewidth=0.5,where='post')
    axis[1].step(MG,resExt[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    axis[1].step(MG,resExt[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    axis[1].fill_between(MG,resExt[0],resExt[1],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='post')

    axis[1].step(MR[-2:],res[-2:], c=color[i], linewidth=0.5,where='pre')
    axis[1].step(MR[-2:],resExt[0][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].step(MR[-2:],resExt[1][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='pre')
    
    i+=1
    
    MG, MR, res, resExt= loadataMG("NLO",res,resExt)
    
    axis[0].fill_between(MG,resExt[0],resExt[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],label=r'${\rm \mathbf{'+labels[i]+'}}$ ',step='post')
    axis[0].step(MG,res, c=color[i], linewidth=0.5,where='post')
    axis[0].step(MG,resExt[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    axis[0].step(MG,resExt[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')

    
    axis[0].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],step='pre')
    axis[0].step(MR[-2:],res[-2:], c=color[i], linewidth=0.5,where='pre')
    axis[0].step(MR[-2:],resExt[0][-2:], c=color[i][-2:], linewidth=0.5, alpha=0.3,where='pre')
    axis[0].step(MR[-2:],resExt[1][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')

    for a in range(len(resExt)):
        resExt[a]/=resLOS#MG
    res/=resLOS#MG

    axis[1].step(MG,res, c=color[i], linewidth=0.5,where='post')
    axis[1].step(MG,resExt[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    axis[1].step(MG,resExt[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    axis[1].fill_between(MG,resExt[0],resExt[1],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='post')
    
    axis[1].step(MR[-2:],res[-2:], c=color[i], linewidth=0.5,where='pre')
    axis[1].step(MR[-2:],resExt[0][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].step(MR[-2:],resExt[1][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='pre')


    axis[0].legend(loc='best', fontsize = 12)
        
    axis[0].set_yscale("log")
    #axis[0].set_xscale("log")
    #axis[1].set_xscale("log")
    axis[0].set_xticklabels([])
    
    axis[0].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='y',which='both',direction='in')
    axis[0].tick_params(axis='y',which='both',direction='in')
   
    
    axis[0].grid(linestyle='dashed',linewidth=0.5)
    axis[1].grid(linestyle='dashed',linewidth=0.5)
    
    axis[0].set_ylabel(r'$\dfrac{{\rm d} \sigma}{{\rm d} M}$ [pb/GeV]',fontsize=16)
    axis[1].set_ylabel(r'$\sigma/\sigma_{{\rm ref}}$',fontsize=16)
    axis[1].set_xlabel(r'${\rm M}$ [GeV]',fontsize=16)
    
    #axis[0].set_xlim(right=2e4, left=200)
    #axis[1].set_xlim(right=2e4, left=200)
    
    figure.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('XSec_ttxdM.svg')
    plt.savefig('XSec_ttxdM.png')
    plt.legend([])
    plt.close()
 
BuildPlot()

    
