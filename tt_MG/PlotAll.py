import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams.update({'font.size': 14})

#types=["Born","NLO","Diff"]
types=["NLO","Diff"]
labels=["LO","NLO","NLO+NLL","LO MG","NLO MG"]
labels=["NLO","NLO+NLL","NLO MG"]
color=['b','g','k','r','orange','m','y','c']

def loadataMycode(order,res,resExt):
    Mycode=np.loadtxt("MyCodeResum/Results/"+order+"Z_muF0_muR0.txt",delimiter=",")
    s=np.transpose(Mycode)[1]
    s=np.sqrt(s)
    res=np.transpose(Mycode)[0]
    resExt=np.zeros((2,len(s)))
    tmp=np.zeros((7,len(s)))
    kid=0
    for muF in range(-1,2):
        for muR in range(-1,2):
            if(muR*muF>=0):
                Mycode=np.loadtxt("MyCodeResum/Results/"+order+"Z_muF"+str(muF)+"_muR"+str(muR)+".txt",delimiter=",")
                tmp[kid]=np.transpose(Mycode)[0]
                kid+=1
                   
    for k in range(len(s)):
        resExt[0][k]=min(tmp[p][k] for p in range(len(tmp)))
        resExt[1][k]=max(tmp[p][k] for p in range(len(tmp)))
    return(s,res,resExt)
        

def loadataMG(order,res,resExt):
    MG=np.loadtxt('MGfO/dataTot/MGXsec_Z.txt',delimiter=",")
    s=np.transpose(MG)[0]
    resExt=np.zeros((2,len(s)))
    res=np.transpose(MG)[1]
    indixNLO=np.array([1,5,6,7,8,9,10])
    indixLO=indixNLO+10
    if(order=="NLO"):
        lim=indixNLO
    else:
        lim=indixLO
    
    for k in range(len(s)):
        resExt[0][k]=min(np.transpose(MG)[p][k] for p in lim)
        resExt[1][k]=max(np.transpose(MG)[p][k] for p in lim)
    return(s,res,resExt)
        

def BuildPlot():
    res=1
    resExt=1
    s=1
    resS=1
    resExtS=1
    i=0

    figure, axis=plt.subplots(2,gridspec_kw={'height_ratios': [1.67, 1]})

    for order in types:
        s, resS, resExtS = loadataMycode(order,resS,resExtS)

        if(order=="NLO"):
            resLOS=np.copy(resS)
            
        axis[0].fill_between(s,resExtS[0],resExtS[1],alpha=0.3,edgecolor=color[i],linewidth=1,antialiased=True,facecolor=color[i],label=r''+labels[i]+' ')
        axis[0].plot(s,resS,alpha=0.5,linewidth=1,antialiased=True,color=color[i])
        
        for a in range(len(resExtS)):
            resExtS[a]/=resLOS
        resS/=resLOS

        axis[1].plot(s,resS, c=color[i], linewidth=0.5)
        axis[1].fill_between(s,resExtS[0],resExtS[1],alpha=0.3,edgecolor=color[i],facecolor=color[i],linewidth=1,antialiased=True)
        print(order, " for my code  + ",np.mean((resExtS[1]-resS)/resS)*100," - ",np.mean((resS-resExtS[0])/resS)*100, "%")       
        
        i+=1
        
    # s, res, resExt= loadataMG("LO",res,resExt)
    # #resLOMG=np.copy(res)
    
    # axis[0].fill_between(s,resExt[0],resExt[1],alpha=0.3,edgecolor=color[i],linewidth=1,antialiased=True,facecolor=color[i],label=r''+labels[i]+' ')
    # axis[0].plot(s,res,alpha=0.3,linewidth=1,antialiased=True,color=color[i])    

    # for a in range(len(resExt)):
    #     resExt[a]/=resLOS#MG
    # res/=resLOS#MG

    # axis[1].plot(s,res, c=color[i], linewidth=0.5)
    # axis[1].fill_between(s,resExt[0],resExt[1],alpha=0.3,edgecolor=color[i],facecolor=color[i],linewidth=1,antialiased=True)

    # i+=1
    
    s, res, resExt= loadataMG("NLO",res,resExt)
    
    axis[0].fill_between(s,resExt[0],resExt[1],alpha=0.3,edgecolor=color[i],linewidth=1,antialiased=True,facecolor=color[i],label=r''+labels[i]+' ')
    axis[0].plot(s,res,alpha=0.3,linewidth=1,antialiased=True,color=color[i])    

    for a in range(len(resExt)):
        resExt[a]/=resLOS#MG
    res/=resLOS#MG

    axis[1].plot(s,res, c=color[i], linewidth=0.5)
    axis[1].fill_between(s,resExt[0],resExt[1],alpha=0.3,edgecolor=color[i],facecolor=color[i],linewidth=1,antialiased=True)


    axis[0].legend(loc='best', fontsize = 12)
        
    axis[0].set_yscale("log")
    axis[0].set_xscale("log")
    axis[1].set_xscale("log")
    axis[0].set_xticklabels([])
    
    axis[0].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='y',which='both',direction='in')
    axis[0].tick_params(axis='y',which='both',direction='in')
   
    
    axis[0].grid(linestyle='dashed',linewidth=0.5)
    axis[1].grid(linestyle='dashed',linewidth=0.5)
    
    axis[0].set_ylabel(r'$\sigma$ [pb]',fontsize=16)
    axis[1].set_ylabel(r'$\sigma/\sigma_{NLO}$',fontsize=16)
    axis[1].set_xlabel(r'$\sqrt{s}$ [GeV]',fontsize=16)
    
    axis[0].set_xlim(right=2e4, left=200)
    axis[1].set_xlim(right=2e4, left=200)
    
    figure.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('XSec_Zprod.svg')
    plt.savefig('XSec_Zprod.png')
    plt.legend([])
    plt.close()
 
BuildPlot()

    
