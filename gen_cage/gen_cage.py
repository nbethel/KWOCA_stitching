from glob import glob
import numpy as np
z1=0
z2=0
zc1=0
zc2=0
box_cent=184.0
chains=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
pdb='full_asu.pdb'
fout=open('D5_cage.pdb','w')
for i in range(5):
  for j in range(2):
    xxx=np.zeros((1,3))
    cnt=0
    fin=open(pdb,'r')
    lline=fin.readline()
    rotz=i*360.0/5.0*np.pi/180.0
    roty=np.pi*j
    #rot=0
    while len(lline)>0:
        if len(lline)>4 and lline[:4]=="ATOM":
            if lline[21]=='A':
                ind=0
            elif lline[21]=='B':
                ind=1
            elif lline[21]=='C':
                ind=2
            elif lline[21]=='D':
                ind=3
            elif lline[21]=='E':
                ind=4
            elif lline[21]=='F':
                ind=5
            elif lline[21]=='G':
                ind=6
            elif lline[21]=='H':
                ind=7
            elif lline[21]=='I':
                ind=8
            xxx[0,0]=float(lline[30:38])-box_cent
            xxx[0,1]=float(lline[38:46])-box_cent
            xxx[0,2]=float(lline[46:54])-box_cent
            Rz=np.array([[np.cos(rotz),np.sin(rotz),0],[-np.sin(rotz),np.cos(rotz),0],[0.0,0.0,1.0]])
            Ry=np.array([[1,0,0],[0,np.cos(roty),np.sin(roty)],[0,-np.sin(roty),np.cos(roty)]])
            xxx=np.transpose(np.matmul(Ry,np.transpose(xxx)))
            xxx=np.transpose(np.matmul(Rz,np.transpose(xxx)))+box_cent
            fout.write(lline[:21]+chains[ind+i*2*3+j*3]+lline[22:30]+'%8.3f%8.3f%8.3f'%(xxx[0,0],xxx[0,1],xxx[0,2])+lline[54:])
            #if rotx==0:
            #    fout2.write(lline[:30]+'%8.3f%8.3f%8.3f'%(xxx[0,0],xxx[0,1],xxx[0,2])+lline[54:])
            #else:
            #    fout.write(lline[:30]+'%8.3f%8.3f%8.3f'%(xxx[0,0],xxx[0,1],xxx[0,2])+lline[54:])
        elif len(lline)>2 and lline[:3]=='TER':
            fout.write(lline)
        lline=fin.readline()
    fin.close()

fout.close()
