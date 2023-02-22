from re import split
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os import mkdir
from os import path
#Function to plot spectrum
#t_name: the path of dipolefile
#out: The pics would be stored in "./out"
#move,N_t:The range of steps to do the fft
def Sp_plot(t_name,out,file,move=0,N_t=-1,dt=0.0034,direction=-1):
    Dipole=dipole_read(t_name,out,file,direction=direction)
    if(N_t==-1):
        N_t=Dipole.shape[1]
    N_t-=move
    lambda_x=4.13*np.linspace(0,N_t-1,N_t)/(dt*N_t)
    fig,ax=plt.subplots()
    if(direction==1 or direction==-1):
        alpha_x=np.fft.fft(Dipole[0][move:N_t+move])
        S_x=np.abs(alpha_x.imag)
        ax.plot(lambda_x,S_x[:],label='X')
    if(direction==2 or direction==-1):
        alpha_y=np.fft.fft(Dipole[1][move:N_t+move])
        S_y=np.abs(alpha_y.imag)
        ax.plot(lambda_x,S_y[:],label='Y')
    if(direction==3 or direction==-1):
        alpha_z=np.fft.fft(Dipole[2][move:N_t+move])
        S_z=np.abs(alpha_z.imag)
        ax.plot(lambda_x,S_z[:],label='Z')
    #write peaks
    ax.set_ylabel('Strength')
    ax.set_xlabel('frequency(eV)')
    ax.set_title("Spectra")
    ax.set_xlim(0,25)
    plt.legend()
    plt.savefig('./'+out+'/'+file+'spec.png')
#function to read in dipole_file, and return 3*nstep 2D array
def dipole_read(total_name,out,file,dt=0.0034,direction=-1):
    #Initial
    list_d=[[] for i in range(3)]
    #N_t=nstep-move_step                 #step used in fft
    #Read in cycle
    m=0                                 #Dimension label
    line=1
    with open(total_name,mode='r')as f1:
        while(line!=''):
            line=f1.readline()
            if(line==''):break
            tmp=line.split()
            list_d[0].append(float(tmp[1]))
            list_d[1].append(float(tmp[2]))
            list_d[2].append(float(tmp[3]))
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

#main part
dir='OUT.ABACUS'                        #The directory of the dipole file
out='PIC'                               #The path you want to ouput your file, equals to directory "./out"
if(not path.exists('./'+out)):
    mkdir(out)
#single plot
file='SPIN1_DIPOLE'                     #filename of dipole file
t_name="./"+dir+"/"+file
Sp_plot(t_name,out,file,dt=0.0034)
