import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams.update({'font.size': 14})

labels=["LO"]+["Diff_2","Diff_1"]
types=["Born","Diff","Expand"]

color=['b','g','k','r','orange','m','y','c']


def loadgg(order,res):
    Mycode=np.loadtxt("SubProcesses/P2_gg_ttx/Results/"+order+"_alpha.txt",delimiter=",")
    alpha=np.transpose(Mycode)[-1]
    res=[np.transpose(Mycode)[p] for p in range(len(np.transpose(Mycode))-1)]
    
    return(alpha,res)

def loadqqb(order,res):
    Mycode=np.loadtxt("SubProcesses/P0_uux_ttx/Results/"+order+"_alpha.txt",delimiter=",")
    alpha=np.transpose(Mycode)[-1]
    res=[np.transpose(Mycode)[p] for p in range(len(np.transpose(Mycode))-1)]

    return(alpha,res)


def BuildPlot():
    figure, axis=plt.subplots(2,gridspec_kw={'height_ratios': [1.67, 1]})    
    for i in range(len(types)):
        resgg=1
        resqqb=1
        resTot=1
        resExt=1
        alpha=1
        
        
        alpha, resgg = loadgg(types[i],resgg)
        #alpha, resqqb = loadqqb(types[i],resqqb)

        #resTot=resqqb+resgg
        resTot=resgg
        #resTot=resqqb

        print(resTot[0])
        if(types[i]=="Born"):
            resBorn=np.copy(resTot[0])
        
            
        resExt=np.zeros((2,len(resTot[0])))
        for k in range(len(resTot[0])):
            resExt[0][k]=min(abs(resTot[p][k]) for p in range(len(resTot)))
            resExt[1][k]=max(abs(resTot[p][k]) for p in range(len(resTot)))
        

        if(types[i]!="Born"):
            axis[0].plot(alpha,abs(resTot[0]),label=r''+labels[i]+' ', c=color[i], linewidth=0.5)
        
        for k in range(len(resExt)):
            resExt[k]/=resBorn
        for p in range(len(resTot)):
            resTot[p]/=resBorn

        if(types[i]!="Born"):
            axis[1].plot(alpha,abs(resTot[0]),label=r''+labels[i]+' ', c=color[i], linewidth=0.5)
            #axis[1].fill_between(alpha,resExt[0],resExt[1],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True)

    
    axis[0].legend(loc='best', fontsize = 12)
        
    axis[0].set_yscale("log")
    axis[0].set_xscale("log")
    axis[1].set_xscale("log")
    axis[1].set_yscale("log")
    axis[0].set_xticklabels([])
    
    axis[0].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='y',which='both',direction='in')
    axis[0].tick_params(axis='y',which='both',direction='in')
   
    
    axis[0].grid(linestyle='dashed',linewidth=0.5)
    axis[1].grid(linestyle='dashed',linewidth=0.5)
    
    axis[0].set_ylabel(r'$\dfrac{{\rm d} \sigma}{{\rm d} M}$ [pb/GeV]',fontsize=16)
    axis[1].set_ylabel(r'$\sigma/\sigma_{{\rm ref}}$',fontsize=16)
    axis[1].set_xlabel(r'$\alpha_s$',fontsize=16)
        
    figure.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('alphaSTo0.svg')
    plt.savefig('alphaSTo0.png')
    plt.legend([])
    plt.close()
 
BuildPlot()

    
