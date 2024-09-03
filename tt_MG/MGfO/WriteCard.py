import sys
import copy

def EnergyCard(path,sqrt_s,j):
    path="/home/yehudi/MG5_aMC_v3_3_0/"+path+"/Cards/"
    file = open(path+"run_card.dat", "r")
    replaced_content = ""
    i=0
    ji=str(j)
    #ij=str((j+1)*1000)
    word=" ebeam"
    pdf="pdlabel"
    pdfw="reweight_PDF"
    pdfID="lhaid"
    acc="req_acc_FO"
    Nev=" nevents "
    cuts="mxx_min"
    #cutmax="mxx_max"
    dynscale=" dynamical_scale_choice "
    #PS="parton_shower"
    scale="muF_over_ref"
    
    for line in file:
        if(line.find(word)>0):
            i+=1
            new_line = "     " + str(sqrt_s/2) + " = ebeam" + str(i) +" "
            replaced_content += new_line + "\n"
        elif(line.find(pdf)>0 and line.find(pdfID)<0):
            replaced_content += line.replace("nn23nlo","lhapdf")
        elif(line.find(pdfID)>0):
            #replaced_content += line.replace("244600","14400")
            replaced_content += line.replace("244600","21200")
        elif(line.find(dynscale)>0):
            replaced_content += line.replace("-1","10")
        elif(line.find(acc)>0):
            replaced_content += line.replace("0.01","0.0001")
        elif(line.find(Nev)>0):
            replaced_content += line.replace(" 10000 "," 100000 ")
        elif(line.find(pdfw)>0):
            replaced_content += line.replace("False","True")
        elif(line.find(scale)>0):
            replaced_content += line.replace("1.0","1.0")
        elif(line.find(cuts)>0):
            replaced_content += line.replace("{}","{6: "+ji+"}")
            #replaced_content += "  {6: "+ji+"} = mxx_min_pdg"+"\n"
        #     if(j==1):
        #         replaced_content += "  {6: "+ij+"} = mxx_max_pdg"+"\n"
        # elif(line.find(cutmax)>0):
        #     replaced_content += "  {6: "+ij+"} = mxx_max_pdg"+"\n"
        else:
            replaced_content += line 
            
    file.close()
    
    write_file = open(path+"run_card.dat", "w")
    write_file.write(replaced_content)
    write_file.close()
    

ic=int(sys.argv[2])
ntot=int(sys.argv[3])
path=sys.argv[1]
#s0=13000
s0=7000
#s0=8000
s=s0

# xmi=1000
# xms=1500
#xmi=2*172.7
xmi=2*173.2
xms=500
#xms=4000+2*173.2

EnergyCard(path,s,xmi)


def FOCard(path):
    path="/home/yehudi/MG5_aMC_v3_3_0/"+path+"/Cards/"

    fin = open(path+"FO_analyse_card.dat", "rt")
    fout = open(path+"FO_analyse_card2.dat", "wt")

    for line in fin:
        fout.write(line.replace('template', 'pp_ttx_v2'))

    fin.close()
    fout.close()

FOCard(path)


def FOHwU(path,xmi,xms,nbins):
    path="/home/yehudi/MG5_aMC_v3_3_0/"+path+"/FixedOrderAnalysis/"
    word="xms="
    wordi="xmi="
    setbins=",xmi,"
    ct=""
    fin = open(path+"analysis_HwU_pp_ttx_v2.f", "rt")

    for line in fin:
        if(line.find(word)>0):
            ct += "      xms="+str(xms)+"d0"+"\n"
        elif(line.find(wordi)>0):
            ct += "      xmi="+str(xmi)+"d0"+"\n"
            #ct += "      xmi=1000d0"+"\n"
        elif(line.find(setbins)>0):
            ct +="          call HwU_book(l+ 5,'m inv        '//cc(i),"+str(nbins)+",xmi,xms)"+"\n"
        else:
            ct+=line
    fin.close()
    write_file= open(path+"analysis_HwU_pp_ttx_v2.f","w")
 
    write_file.write(ct)
    write_file.close()

#FOHwU(path,s)
n=int(sys.argv[4])
FOHwU(path,xmi,xms,n)


def scales(path,n,mi,mf):
    path="/home/yehudi/MG5_aMC_v3_3_0/"+path+"/SubProcesses/"
    word="nbin="
    wordi="mi="
    wordf="mf="
    ct=""
    fin = open(path+"setscales.f", "rt")

    for line in fin:
        if(line.find(wordi)>0):
            ct += "         mi="+str(mi)+"d0"+"\n"
        elif(line.find(wordf)>0):
            ct += "         mf="+str(mf)+"d0"+"\n"
        elif(line.find(word)>0):
            ct += "         nbin="+str(n)+"\n"
        else:
            ct+=line
    fin.close()
    write_file= open(path+"setscales.f","w")
 
    write_file.write(ct)
    write_file.close()

#xmi=1000
#xms=500
scales(path,ntot*n,xmi,xms)

def cuts(path,xms):
    path="/home/yehudi/MG5_aMC_v3_3_0/"+path+"/SubProcesses/"
    wordf="Minv2.gt"
    ct=""
    fin = open(path+"cuts.f", "rt")

    for line in fin:
        if(line.find(wordf)>0):
            ct +="      if (Minv2.gt."+str(xms)+"d0**2) then"+"\n"
        else:
            ct+=line
    fin.close()
    write_file= open(path+"cuts.f","w")
 
    write_file.write(ct)
    write_file.close()

cuts(path,xms)
