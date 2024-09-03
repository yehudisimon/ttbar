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

labels=["LO"]+["NLL"]+["LO \, MG","NLO"]+["NLO \oplus NLL"]
types=["Born","Diff"]
#types=["Born","Resum"]
color=['r','k','b','b','magenta','y','b','orange']

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
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/BkpLO.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/BkpLOqq.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/LOqq2.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/BkpLOgg.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/LOgg2.txt',delimiter=" ")
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
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/BkpNLO.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/BkpNLOqq.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/NLOqq2.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/BkpNLOgg.txt',delimiter=" ")
        #MadG=np.loadtxt('../tt_VaryMG/MGfO/data/NLOgg2.txt',delimiter=" ")
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
    figure, axis=plt.subplots(2,gridspec_kw={'height_ratios': [1.67, 1]})    
    ntot=1
    resMG=resMGExt=resMGNLO=resMGExtNLO=1
    resMu=1
    a=b=aNLO=bNLO=1
    
    a, b, resMG, resMGExt = loadMGLO(resMG,resMGExt,ntot)
    aNLO, bNLO, resMGNLO, resMGExtNLO, resMu = loadMGNLO(resMGNLO,resMGExtNLO,ntot) 

    resBornMG=np.copy(resMG)
    resNLO=np.copy(resMGNLO)

    DeltaNLO=np.mean((resMGExtNLO[1]-resMGExtNLO[0])/resMGNLO)
    
    # axis[0].step(a,resMG, c=color[len(types)], linewidth=0.5,where='post')
    # axis[0].fill_between(a,resMGExt[0],resMGExt[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)],label=r'${\rm \mathbf{'+labels[len(types)]+'}}$ ',step='post')
    # axis[0].step(a,resMGExt[0], c=color[len(types)], linewidth=0.5, alpha=0.3,where='post')
    # axis[0].step(a,resMGExt[1], c=color[len(types)], linewidth=0.5, alpha=0.3,where='post')

    # axis[0].step(b[-2:],resMG[-2:], c=color[len(types)], linewidth=0.5,where='pre')
    # axis[0].fill_between(b[-2:],resMGExt[0][-2:],resMGExt[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)],step='pre')
    # axis[0].step(b[-2:],resMGExt[0][-2:], c=color[len(types)], linewidth=0.5, alpha=0.3,where='pre')
    # axis[0].step(b[-2:],resMGExt[1][-2:], c=color[len(types)], linewidth=0.5, alpha=0.3,where='pre')

    axis[0].step(aNLO,resMGNLO, c=color[len(types)+1], linewidth=0.5,where='post')
    axis[0].fill_between(aNLO,resMGExtNLO[0],resMGExtNLO[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+1],label=r'${\rm \mathbf{'+labels[len(types)+1]+'}}$ ',step='post')
    axis[0].step(aNLO,resMGExtNLO[0], c=color[len(types)+1], linewidth=0.5, alpha=0.3,where='post')
    axis[0].step(aNLO,resMGExtNLO[1], c=color[len(types)+1], linewidth=0.5, alpha=0.3,where='post')

    axis[0].step(bNLO[-2:],resMGNLO[-2:], c=color[len(types)+1], linewidth=0.5,where='pre')
    axis[0].fill_between(bNLO[-2:],resMGExtNLO[0][-2:],resMGExtNLO[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[len(types)+1],step='pre')
    axis[0].step(bNLO[-2:],resMGExtNLO[0][-2:], c=color[len(types)+1], linewidth=0.5, alpha=0.3,where='pre')
    axis[0].step(bNLO[-2:],resMGExtNLO[1][-2:], c=color[len(types)+1], linewidth=0.5, alpha=0.3,where='pre')

   
    resgg=resqqb=resTot=resExt=1
    ML=MR=1
    i=0    
    ML, MR, resgg = loadgg("Born",resgg)
    ML, MR, resqqb = loadqqb("Born",resqqb)
    
    resTot=resgg+resqqb
    #resTot=resgg
    #resTot=resqqb
    
    #print(ML,a)
    
    #resBorn=np.copy(resTot[1])
    resBorn=np.copy(resTot[3])           
    
    resExt=np.zeros((2,len(resTot[0])))
    for k in range(len(resTot[0])):
        resExt[0][k]=min(resTot[p][k] for p in indixMe)
        resExt[1][k]=max(resTot[p][k] for p in indixMe)
        

    #axis[0].step(ML,resTot[1], c=color[i], linewidth=0.5,where='post')
    axis[0].step(ML,resTot[3], c=color[i], linewidth=0.5,where='post')
    axis[0].fill_between(ML,resExt[0],resExt[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],label=r'${\rm \mathbf{'+labels[i]+'}}$ ',step='post')
    axis[0].step(ML,resExt[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')
    axis[0].step(ML,resExt[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')

    #axis[0].step(MR[-2:],resTot[1][-2:], c=color[i], linewidth=0.5,where='pre')
    axis[0].step(MR[-2:],resTot[3][-2:], c=color[i], linewidth=0.5,where='pre')
    axis[0].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[i],step='pre')
    axis[0].step(MR[-2:],resExt[0][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
    axis[0].step(MR[-2:],resExt[1][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')
    
    for k in range(len(resExt)):
        resExt[k]/=resBorn
        resMGExt[k]/=resBorn
        resMGExtNLO[k]/=resBorn
        
    for p in range(len(resTot)):
        resTot[p]/=resBorn

    resMG/=resBornMG
    resMGNLO/=resBornMG

    print(resMGNLO)
    
    print(np.mean(resBornMG/resBorn))
    #axis[1].step(ML,resTot[1],label=r''+labels[i]+' ', c=color[i], linewidth=0.5,where='post')
    axis[1].step(ML,resTot[3],label=r''+labels[i]+' ', c=color[i], linewidth=0.5,where='post')
    axis[1].step(ML,resExt[0], c=color[i], linewidth=0.5, alpha=0.3,where='post')             
    axis[1].step(ML,resExt[1], c=color[i], linewidth=0.5, alpha=0.3,where='post')             
    axis[1].fill_between(ML,resExt[0],resExt[1],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='post')                                
    
    #axis[1].step(MR[-2:],resTot[1][-2:], c=color[i], linewidth=0.5,where='pre')
    axis[1].step(MR[-2:],resTot[3][-2:], c=color[i], linewidth=0.5,where='pre')               
    axis[1].step(MR[-2:],resExt[0][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')    
    axis[1].step(MR[-2:],resExt[1][-2:], c=color[i], linewidth=0.5, alpha=0.3,where='pre')    
    axis[1].fill_between(MR[-2:],resExt[0][-2:],resExt[1][-2:],alpha=0.5,facecolor=color[i],linewidth=1,antialiased=True,step='pre')
    
    
    # axis[1].step(a,resBornMG/resBorn, c=color[i+len(types)], linewidth=0.5, where='post')
    # axis[1].step(a,resMGExt[0], c=color[i+len(types)], linewidth=0.5, alpha=0.3,where='post')
    # axis[1].step(a,resMGExt[1], c=color[i+len(types)], linewidth=0.5, alpha=0.3,where='post')
    # axis[1].fill_between(a,resMGExt[0],resMGExt[1],alpha=0.5,facecolor=color[i+len(types)],linewidth=1,antialiased=True,step='post')

    # axis[1].step(b[-2:],resBornMG[-2:]/resBorn[-2:], c=color[i+len(types)], linewidth=0.5, where='pre')
    # axis[1].step(b[-2:],resMGExt[0][-2:], c=color[i+len(types)], linewidth=0.5, alpha=0.3,where='pre')
    # axis[1].step(b[-2:],resMGExt[1][-2:], c=color[i+len(types)], linewidth=0.5, alpha=0.3,where='pre')
    # axis[1].fill_between(b[-2:],resMGExt[0][-2:],resMGExt[1][-2:],alpha=0.5,facecolor=color[i+len(types)],linewidth=1,antialiased=True,step='pre')

    axis[1].step(aNLO,resMGNLO/resBorn*resBornMG, c=color[i+len(types)+1], linewidth=0.5, where='post')
    axis[1].step(aNLO,resMGExtNLO[0], c=color[i+len(types)+1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].step(aNLO,resMGExtNLO[1], c=color[i+len(types)+1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].fill_between(aNLO,resMGExtNLO[0],resMGExtNLO[1],alpha=0.5,facecolor=color[i+len(types)+1],linewidth=1,antialiased=True,step='post')

    axis[1].step(bNLO[-2:],resMGNLO[-2:]/resBorn[-2:]*resBornMG[-2:], c=color[i+len(types)+1], linewidth=0.5, where='pre')
    axis[1].step(bNLO[-2:],resMGExtNLO[0][-2:], c=color[i+len(types)+1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].step(bNLO[-2:],resMGExtNLO[1][-2:], c=color[i+len(types)+1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].fill_between(bNLO[-2:],resMGExtNLO[0][-2:],resMGExtNLO[1][-2:],alpha=0.5,facecolor=color[i+len(types)+1],linewidth=1,antialiased=True,step='pre')
    

    ML, MR, resgg = loadgg("Diff",resgg)
    ML, MR, resqqb = loadqqb("Diff",resqqb)

    resTot=resgg+resqqb
    #resTot=resgg
    #resTot=resqqb
    
    #print([resTot[indixMe[p]] for p in range(7)],resMu[0])
    
    combine=[resMu[p]+resTot[indixMe[p]] for p in range(len(indixMe))]

    
    combineExt=np.zeros((2,len(combine[0])))
    for k in range(len(combine[0])):
        combineExt[0][k]=min(combine[p][k] for p in range(len(indixMe)))
        combineExt[1][k]=max(combine[p][k] for p in range(len(indixMe)))

    DeltaComb=np.mean((combineExt[1]-combineExt[0])/combine[0])
    print("Delta Combination = ",DeltaComb,"Delta NLO = ",DeltaNLO)
    print("Enhancement = ",(combine[0]/resNLO-1)*100," %, Mean = ",np.mean((combine[0]/resNLO-1)*100))

    
    axis[0].step(ML,combine[0], c=color[-1], linewidth=0.5,where='post')
    axis[0].fill_between(ML,combineExt[0],combineExt[1],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[-1],label=r"${\rm \mathbf{NLO \oplus NLL'}}$ ",step='post')
    axis[0].step(ML,combineExt[0], c=color[-1], linewidth=0.5, alpha=0.3,where='post')
    axis[0].step(ML,combineExt[1], c=color[-1], linewidth=0.5, alpha=0.3,where='post')

    axis[0].step(MR[-2:],combine[0][-2:], c=color[-1], linewidth=0.5,where='pre')
    axis[0].fill_between(MR[-2:],combineExt[0][-2:],combineExt[1][-2:],alpha=0.5,linewidth=1,antialiased=True,facecolor=color[-1],step='pre')
    axis[0].step(MR[-2:],combineExt[0][-2:], c=color[-1], linewidth=0.5, alpha=0.3,where='pre')
    axis[0].step(MR[-2:],combineExt[1][-2:], c=color[-1], linewidth=0.5, alpha=0.3,where='pre')


    for k in range(len(combineExt)):
        combineExt[k]/=resBorn
       
    for p in range(len(combine)):
        combine[p]/=resBorn

    print(np.mean(combine))

    axis[1].step(ML,combine[0],label=r''+labels[-1]+' ', c=color[-1], linewidth=0.5,where='post')
    axis[1].step(ML,combineExt[0], c=color[-1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].step(ML,combineExt[1], c=color[-1], linewidth=0.5, alpha=0.3,where='post')
    axis[1].fill_between(ML,combineExt[0],combineExt[1],alpha=0.5,facecolor=color[-1],linewidth=1,antialiased=True,step='post')

    axis[1].step(MR[-2:],combine[0][-2:], c=color[-1], linewidth=0.5,where='pre')
    axis[1].step(MR[-2:],combineExt[0][-2:], c=color[-1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].step(MR[-2:],combineExt[1][-2:], c=color[-1], linewidth=0.5, alpha=0.3,where='pre')
    axis[1].fill_between(MR[-2:],combineExt[0][-2:],combineExt[1][-2:],alpha=0.5,facecolor=color[-1],linewidth=1,antialiased=True,step='pre')

        
    axis[0].legend(loc='best', fontsize = 12)
    #axis[0].set_yscale("log")
    axis[0].set_xticklabels([])

    #axis[0].set_ylim([0,0.9])
    axis[0].set_xlim([2*173.2,500])
    axis[1].set_xlim([2*173.2,500])
    
    axis[0].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='x',which='both',direction='in',bottom=True)
    axis[1].tick_params(axis='y',which='both',direction='in')
    axis[0].tick_params(axis='y',which='both',direction='in')
   
    axis[1].xaxis.set_minor_locator(MultipleLocator(5))
    axis[0].xaxis.set_minor_locator(MultipleLocator(5))
    axis[1].yaxis.set_minor_locator(MultipleLocator(0.1))
    axis[0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    axis[0].grid(linestyle='dashed',linewidth=0.5)
    axis[1].grid(linestyle='dashed',linewidth=0.5)
    
    
    axis[0].set_ylabel(r'$\dfrac{{\rm d} \sigma}{{\rm d} M}$ [pb/GeV]',fontsize=16)
    axis[1].set_ylabel(r'$\sigma/\sigma_{{\rm LO}}$',fontsize=16)
    axis[1].set_xlabel(r'${\rm M_{t\bar{t}}}$ [GeV]',fontsize=16)

    axis[0].set_title(r'$t \bar{t} \; {\rm production \; cross \; section} \sqrt{s}=7 \; {\rm TeV}$', y=1.1, pad=0, x=0.5)

    figure.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    plt.savefig('ComparisonSlices.svg')
    plt.savefig('ComparisonSlices.png',dpi=800)

    plt.savefig('/home/yehudi/Documents/Thesis/Manuscript/chapters/images/chapter05/ComparisonSlices.png',dpi=800)
    plt.legend([])
    plt.close()
 
BuildPlot()

    
