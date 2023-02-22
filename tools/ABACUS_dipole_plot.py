from re import split
import numpy as np
import matplotlib.pyplot as plt
from ase.units import Bohr
from math import pi
from scipy.constants import c
from os import listdir
from os import mkdir
from shutil import move
from shutil import copyfile
from os import path
def Sp_plot(t_name,out,file,move=0,N_t=-1,dt=0.0034,direction=-1):
    Dipole=dipole_read(t_name,out,file,direction=direction)
    if(N_t==-1):
        N_t=Dipole.shape[1]
        N_t-=move
    shift_omega_step=1
    lambda_x=(c/np.linspace(shift_omega_step,N_t-1,N_t-shift_omega_step)*dt*N_t*1e-6)
    fig,ax=plt.subplots()
    if(direction==1 or direction==-1):
        alpha_x=np.fft.fft(Dipole[0][move:N_t])
        S_x=np.abs(alpha_x.imag)
        ax.plot(lambda_x,S_x[shift_omega_step:],label='X')
    if(direction==2 or direction==-1):
        alpha_y=np.fft.fft(Dipole[1][move:N_t])
        S_y=np.abs(alpha_y.imag)
        ax.plot(lambda_x,S_y[shift_omega_step:],label='Y')
    if(direction==3 or direction==-1):
        alpha_z=np.fft.fft(Dipole[2][move:N_t])
        S_z=np.abs(alpha_z.imag)
        ax.plot(lambda_x,S_x[shift_omega_step:],label='Z') 
    ax.set_ylabel('Strength')
    ax.set_xlabel('frequency(eV)')
    ax.set_title("Spectra")
    ax.set_xlim(80,200)
    plt.legend()
    plt.savefig('./'+out+'/'+file+'spec.png')

def Sp_field_plot():
    pass
def dipole_read(total_name,out,file,dt=0.0034,direction=-1):
    #Initial
    list_d=[[] for i in range(3)]
    #N_t=nstep-move_step                 #step used in fft
    #Read in cycle
    m=0                                 #Dimension label
    line=1
    with open(total_name,mode='r')as f1:
        while(line!=''):
            while(m<3):
                line=f1.readline()
                if(line==''):break
                tmp=line.split()
                if (m==0):
                    list_d[0].append(float(tmp[1]))
                elif (m==1):
                    list_d[1].append(float(tmp[1]))
                elif (m==2):
                    list_d[2].append(float(tmp[1]))
                m+=1
            m=0
    #Plot part
    Dipole=np.array(list_d)
    N_t=Dipole.shape[1]
    Ti=np.linspace(0,(N_t-1)*dt,N_t)
    fig,ax=plt.subplots()
    if(direction==1 or direction==-1):
        ax.plot(Ti,Dipole[0],label='X')
    if(direction==2 or direction==-1):
        ax.plot(Ti,Dipole[1],label='Y')
    if(direction==3 or direction==-1):
        ax.plot(Ti,Dipole[2],label='Z')
    ax.set_ylabel('Strength')
    ax.set_xlabel('time(fs)')
    ax.set_title('Dipole')
    plt.legend()
    plt.savefig('./'+out+'/'+file+'dipole.png')
    return Dipole
#reading para
def autoanlysis(filename,dir,out):
    list_dir=listdir("./"+dir+"/")
    count=0
    print("total file: "+str(len(list_dir)))
    #main cycle
    for file in list_dir:
        if(filename in file):
            count+=1
            t_name="./"+dir+"/"+file
            Sp_plot(t_name,out,file)

#main part
filename='dipole.dat'
dir='dipole_nomove'
out='PIC'
file='C2H2_nomove_X_dipole.dat'
if(not path.exists('./'+out)):
    mkdir(out)
#direction='X'
#autoanlysis(filename,dir,out)
#single plot
t_name="./"+dir+"/"+file
Sp_plot(t_name,out,file,direction=1)