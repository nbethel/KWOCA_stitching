import shutil
fin=open('logger1.log','r')
import numpy as np
lline=fin.readline()
thresh=6.0
RR=100.0
while len(lline)>0:
    trim1=lline.split()[1]
    dim1=lline.split()[2]
    rms1=float(lline.split()[0])
    if not (rms1<thresh):
        lline=fin.readline()
        continue
    fin2=open('logger2.log','r')
    lline2=fin2.readline()
    while len(lline2)>0:
        rms2=float(lline2.split()[0])
        dim2=lline2.split()[2]
        if not (rms2<thresh and trim1 in lline2): #or err>3.05:
            lline2=fin2.readline()
            continue
        fin3=open('logger3.log','r')
        lline3=fin3.readline()
        while len(lline3)>0:
            rms3=float(lline3.split()[0])
            trim1b=lline3.split()[1]
            if not (rms3<thresh and dim1 in lline3): #or err>3.05:
                lline3=fin3.readline()
                continue
            elif RR<np.sqrt((rms1*rms1+rms2*rms2+rms3*rms3)/3.0):
                lline3=fin3.readline()
                continue
            else:
                outs=[trim1,dim1,dim2,trim1b]
                RR=np.sqrt((rms1*rms1+rms2*rms2+rms3*rms3)/3.0)
                print(np.sqrt((rms1*rms1+rms2*rms2+rms3*rms3)/3.0),rms1,rms2,rms3,trim1,dim1,trim1b,dim2)
            lline3=fin3.readline()
        fin3.close()
        lline2=fin2.readline()
    fin2.close()
    lline=fin.readline()
fin.close()
shutil.copyfile(outs[0],'T1.pdb')
shutil.copyfile(outs[1],'D1.pdb')
shutil.copyfile(outs[2],'D2.pdb')
shutil.copyfile(outs[3],'T1b.pdb')
