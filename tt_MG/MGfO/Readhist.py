import os
import sys
import numpy as np

def ReadHisTot(path,k):
    name=path[:3]
    wordC="+-"
    wordS="ebeam1"
    track=2
    a = []
    cont=""
    cp=""
    
    if(k<10):
        fileS = open("/home/yehudi/MG5_aMC_v3_3_0/"+path+"/Events/run_0"+str(int(k))+"/run_0"+str(int(k))+"_tag_1_banner.txt", "r")
    else:
        fileS = open("/home/yehudi/MG5_aMC_v3_3_0/"+path+"/Events/run_"+str(int(k))+"/run_"+str(int(k))+"_tag_1_banner.txt", "r")
    
    for line in fileS:
        if(line.find(wordS)>0):
            for word in line.split():
                try:
                    a.append(float(word)*2)
                except ValueError:
                    pass            
                
    fileS.close()
    
    file = open("dataTot/dXsec_"+path+str(k)+".txt", "r")
    tr=["",""]
    histlabel=["MinvNLO","MinvLO"]
    limit="histogram>"
    ml="total rate "
    lm=False
    ct=np.array([7,9])
    for line in file:
        if(line.find(limit)>0):
            ct-=1
        for kit in range(len(tr)):
            if(ct[kit]==0):
                tr[kit]+=line[3:]
        if(lm):
            for word in line.split():
                try:
                    a.append(float(word))
                except ValueError:
                    pass
            lm=False
        if(line.find(ml)>0):
            lm=True

    file.close()
    
    for kut in range(len(tr)):
        tr[kut] = tr[kut].split('\n', 1)[1]                                                       
        tr[kut]=tr[kut].replace("  +"," ")                                                        
        tr[kut]=tr[kut].replace("  "," ")                                                         
        tr[kut]=tr[kut].replace("  "," ")                                                         
        
        write_file = open("data/"+histlabel[kut]+"_"+name+str(k)+"_dXsec.txt", "w")               
        write_file.write(tr[kut])                                                                 
        write_file.close()             


    #indexes of interest: 0 (s) 3 (NLO central) 4 (dy NLO) 9 (PDF_min) 10 (PDF_max) 12-15 (mus) 17 (mu) 19 (mu) ((NON entre 20-79 for PDF)) 82 (LO) 83 (dy LO) 90 (PDF_min) 91 (PDF_max) 93-96 98 100
    l=[0,3,4,9,10,12,13,14,15,17,19] #NLO
    Nlo=78
    l+=[l[p]+Nlo for p in range(1,len(l))] #LO
    b=[a[p] for p in l]
    a=np.copy(b)
    
    cont=""
    for l in range(len(a)-1):
        cont+=str(a[l])+","
    cont+=str(a[len(a)-1])+"\n"
    
    namef="dataTot/MGXsec_"+name+".txt"
    if(os.path.isfile(namef)):
        readf=open(namef,"r")
        for line in readf:
            cp+=line
        cp+=cont
        write_file = open(namef, "w")
        write_file.write(cp)
        write_file.close()

    else:
        write_file = open(namef, "w")
        write_file.write(cont)
        write_file.close()


k=int(sys.argv[2])
path=sys.argv[1]
ReadHisTot(path,k)
