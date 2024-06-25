#include "local_orbital_wfc.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/read_wfc_nao.h"

Local_Orbital_wfc::Local_Orbital_wfc()
{
    allocate_flag = false;
    wfck_flag = false;
    complex_flag = false;
    nks = 0;
}

Local_Orbital_wfc::~Local_Orbital_wfc()
{

    // used for k-points.
    if (this->complex_flag)
    {
        delete[] this->wfc_k_grid2;
    }
    if (this->wfck_flag)
    {
        for (int i = 0; i < nks; i++)
        {
            delete[] this->wfc_k_grid[i];
        }
        delete[] this->wfc_k_grid;
    }
}

int Local_Orbital_wfc::globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

int Local_Orbital_wfc::localIndex(int globalindex, int nblk, int nprocs, int& myproc)
{
    myproc = int((globalindex % (nblk * nprocs)) / nblk);
    return int(globalindex / (nblk * nprocs)) * nblk + globalindex % nblk;
}
