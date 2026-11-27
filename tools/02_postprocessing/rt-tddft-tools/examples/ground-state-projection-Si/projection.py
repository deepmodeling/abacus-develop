import numpy as np
import matplotlib.pyplot as plt
class Projection:
    def __init__(self, stepref, klist, steps, fdir='./OUT.ABACUS', wfc_dir='', s_dir=''):
        self.stepref = stepref
        self.klist = klist
        self.steps = steps
        self.wfc_dir = fdir + wfc_dir
        self.s_dir = fdir + s_dir
        wfc_ref, Ocp_ref = self.read_wfc(klist[0]+1, stepref+1, dir=self.wfc_dir)
        self.nband = len(Ocp_ref)
        self.nlocal = len(wfc_ref[0])

    def read_wfc(self, kpoint, nstep, dir='.'):
        file = dir+'/WFC_NAO_K'+str(kpoint)+'_ION'+str(nstep)+'.txt'
        read_flag=0
        with open(file,'r') as f:
            for line in f:
                if('bands' in line):
                    bands=int(line.split()[0])
                elif('orbitals' in line):
                    orbitals=int(line.split()[0])
                    wavef0=np.zeros([bands,orbitals],dtype=complex)
                    Ocp=np.zeros(bands,dtype=float)
                elif('(band)' in line):
                    read_flag=0
                    i=int(line.split()[0])-1
                    j=0
                elif('Occupations' in line):
                    Ocp[i]=float(line.split()[0])
                    read_flag=1
                    continue
                elif(read_flag==1):
                    tmp=line.split()
                    for l in range(len(tmp)//2):
                        wavef0[i][j]=complex(float(tmp[2*l]),float(tmp[2*l+1]))
                        j+=1
        return wavef0,Ocp

    def CSC(self, Cd,Sk,C):
        Cdag=np.conjugate(Cd).transpose()
        return np.matmul(Cdag,np.matmul(Sk,C))
    def CSC2(self, Cd,Sk,C):
        CSC1=self.CSC(Cd,Sk,C)
        return CSC1*CSC1.conjugate()

    def S_read(self, kpoint,nstep,dir='.'):
        file = dir+'/'+str(nstep)+'_data-'+str(kpoint)+'-S'
        count=0
        fir=1
        with open(file, 'r') as f1:
            for line in f1:
                tmp=line.split()
                if(fir==1):
                    dim=int(tmp[0])
                    i=0
                    Sk=np.zeros([dim,dim],dtype=complex)
                    fir=0
                    for j in range(dim):
                        c_s=eval(tmp[j+1])
                        Sk[i][j]=complex(c_s[0],c_s[1])
                else:
                    for j in range(dim-i):
                        c_s=eval(tmp[j])
                        Sk[i][i+j]=complex(c_s[0],c_s[1])
                i+=1
                if(i==0):
                    break 
        for i in range(dim):
            for j in range(i):
                Sk[i][j]=Sk[j][i].conjugate()
        return Sk

    def cal_Pnm(self, ik, nstep, wfc_ref):
        wfc_t,Ocp_t=self.read_wfc(ik+1, nstep+1, dir=self.wfc_dir)
        S_t=self.S_read(ik,nstep,dir=self.s_dir)
        Pnm=np.zeros([self.nband,self.nband],dtype=float)
        for a in range(self.nband):
            for b in range(self.nband):
                Pnm[a][b] = self.CSC2(wfc_ref[a],S_t,wfc_t[b]).real
        return Pnm, Ocp_t

    def cal_On_single(self, ik, nstep, wfc_ref):
        Pnm, Ocp_t = self.cal_Pnm(ik, nstep, wfc_ref)
        On = np.zeros(self.nband,dtype=float)
        On = Ocp_t @ Pnm
        return On
    
    def save_On_all(self):
        for ik in self.klist:
            wfc_ref,Ocp_ref=self.read_wfc(ik+1, self.stepref+1, dir=self.wfc_dir)
            On_tot = np.zeros([len(steps),self.nband],dtype=float)
            for n, nstep in enumerate(steps):
                On_tot[n] = self.cal_On_single(ik, nstep, wfc_ref)
            np.savetxt('On_'+str(ik)+'.dat', On_tot)


if __name__ == "__main__":
    #the kpoints you need, check kpoints file to get the index, for this example, 0 means gamma point
    klist=[0, 1]
    #the steps you need, check STRU_MD file to get the index
    start_step = 0
    end_step = 10005
    out_interval = 25
    steps = range(start_step, end_step, out_interval)
    #the ground state step
    stepref=0
    #Ground state projection
    pro=Projection(0, klist, steps)
    #save the Occupation infos to file
    pro.save_On_all()