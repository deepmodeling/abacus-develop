# SIAB Package Description


**S**ystematically

**I**mprovable

**A**tomic orbital

**B**asis generator based on spillage formula


# HOW TO USE SIAB

The optimization can choose one of the three minimization methods:

- Simulated Annealing (**SA**),
- PyTorch Gradient (**PTG**),
- PyTorch Gradient with dpsi (**PTG_dpsi**).

The executable files for the three methods are:

- ./SimulatedAnnealing/source/SIA_s.exe, 
- ./PyTorchGradient/source/main.py, 
- ../opt_orb_pytorch_dpsi/main.py, 

respectively.


##  1. Write input file

Firstly, write the input file, such as **ORBITAL_INPUT_DZP** in example-directories, for script **Generate_Orbital_AllInOne.sh**.
All three approachs work with the same bash script and use the same input file.
Please use **absolute path** for each file/directory in input file.


##  2. Set up dependence environment

Secondly, we set up the dependence environment for ABACUS and SIAB, such as:

```bash
$ module load hpcx/2.9.0/hpcx-intel-2019.update5 mkl/2019.update5 elpa/2019.05.002/hpcx-intel-2019.update5
```

Especially for SIAB with **PyTorch Gradient** approach, we need pytorch v1.1.0.


### How to install pytorch:

Take the HanHai20@USTC system for example:

```bash
$ module load gcc/7.5.0min      #optional, maybe unnecessary.
$ module load anaconda3_nompi
$ module list
Currently Loaded Modulefiles:
  1) elpa/2019.05.002/hpcx-intel-2019.update5   4) hpcx/2.9.0/hpcx-intel-2019.update5         7) libxc/4.3.4/hpcx-intel-2019.update5
  2) gcc/7.5.0min                               5) mkl/2019.update5
  3) intel/2019.update5                         6) anaconda3_nompi
$ python3 -V
Python 3.7.4

$ conda create -n pytorch110 python=3.7
$ source activate pytorch110    #or: conda activate pytorch110
$ conda install pytorch torchvision torchaudio cpuonly -c pytorch
$ source deactivate             #or: conda deactivate

$ source activate pytorch110  #or: conda activate pytorch110
$ pip3 install --user scipy numpy
$ pip3 install --user torch_optimizer
```


## 3. Run generation

Finally, `cd` into an example folder, and run command like this:

```bash
$ ../Generate_Orbital_AllInOne.sh ORBITAL_INPUT_DZP
 or
$ bsub -q idle -n 8 -oo running.log ../Generate_Orbital_AllInOne.sh ORBITAL_INPUT_DZP
```

