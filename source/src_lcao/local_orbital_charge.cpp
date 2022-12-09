#include "local_orbital_charge.h"

#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"
#include "../src_pw/global.h"

// Shen Yu add 2019/5/9
extern "C"
{
    void Cblacs_gridinfo(int icontxt, int *nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
    int Cblacs_pnum(int icontxt, int prow, int pcol);
}

int Local_Orbital_Charge::out_dm = 0;
int Local_Orbital_Charge::out_dm1 = 0;

Local_Orbital_Charge::Local_Orbital_Charge()
{
    // for gamma algorithms.
    this->init_DM = false;
    this->lgd_now = 0;
    this->lgd_last = 0;

    // for k-dependent algorithms.
    this->init_DM_R = false;

    // xiaohui add 2014-06-19
    // band_local = nullptr;
    // Z_wg = nullptr;
    // Z_LOC = nullptr;
    sender_2D_index = nullptr;
    sender_size_process = nullptr;
    sender_displacement_process = nullptr;
    sender_buffer = nullptr;

    receiver_local_index = nullptr;
    receiver_size_process = nullptr;
    receiver_displacement_process = nullptr;
    receiver_buffer = nullptr;
}

Local_Orbital_Charge::~Local_Orbital_Charge()
{
    // with gamma point only
    if (this->init_DM)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            delete[] DM[is];
            delete[] DM_pool[is];
        }
        delete[] DM;
        delete[] DM_pool;
    }
    delete[] sender_2D_index;
    delete[] sender_size_process;
    delete[] sender_displacement_process;

    delete[] receiver_local_index;
    delete[] receiver_size_process;
    delete[] receiver_displacement_process;
    // with k points
    if (this->init_DM_R)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            delete[] DM_R[is];
        }
        delete[] DM_R;
    }
}

void Local_Orbital_Charge::allocate_dm_wfc(const int &lgd,
                                           elecstate::ElecState *pelec,
                                           Local_Orbital_wfc &lowf,
                                           psi::Psi<double> *psid,
                                           psi::Psi<std::complex<double>> *psi)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "allocate_dm_wfc");

    this->LOWF = &lowf;
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        // here we reset the density matrix dimension.
        this->allocate_gamma(lgd, psid, pelec);
    }
    else
    {
        lowf.allocate_k(lgd, psi, pelec);
        this->allocate_DM_k();
    }

    return;
}
