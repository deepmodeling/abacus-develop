#include "sto_wf.h"
#include<cassert>
#include "time.h"

//---------Temporary------------------------------------
#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//------------------------------------------------------

Stochastic_WF::Stochastic_WF()
{
    chiortho = nullptr;
    chi0 = nullptr;
    shchi = nullptr;
    nchip = nullptr;
}

Stochastic_WF::~Stochastic_WF()
{
    delete[] chi0;
    delete[] shchi;
    delete[] chiortho;
    delete[] nchip;
}

void Stochastic_WF::init(K_Vectors* p_kv, const int npwx_in)
{
    this->nks = p_kv->nks;
    this->ngk = p_kv->ngk.data();
    this->npwx = npwx_in;
    chi0 = new ModuleBase::ComplexMatrix[nks];
    shchi = new ModuleBase::ComplexMatrix[nks];
    chiortho = new ModuleBase::ComplexMatrix[nks];
    nchip = new int[nks];
    
    if (nks <= 0)
    {
        ModuleBase::WARNING_QUIT("Stochastic_WF", "nks <=0!");
    }
}

void Init_Sto_Orbitals
(
    Stochastic_WF& stowf, 
    const int seed_in,
    const int nks
)
{
    if (seed_in == 0 || seed_in == -1)
    {
        srand((unsigned)time(nullptr) + GlobalV::MY_RANK * 10000); // GlobalV global variables are reserved 
    }
    else
    {
        srand((unsigned)std::abs(seed_in) + GlobalV::MY_RANK * 10000);
    }

    bool firstrankmore = false;
    int igroup = 0;
    // I am not sure which is better.
    // former processor calculate more bands
    if (firstrankmore)
    {
        igroup = GlobalV::MY_STOGROUP;
    }
    // latter processor calculate more bands
    else
    {
        igroup = GlobalV::NSTOGROUP - GlobalV::MY_STOGROUP - 1;
    }
    const int nchi = INPUT.nbands_sto;
    const int ndim = stowf.npwx;
    const int ngroup = GlobalV::NSTOGROUP;
    if (ngroup <= 0)
    {
        ModuleBase::WARNING_QUIT("Init_Sto_Orbitals","ngroup <= 0!");
    }
    int tmpnchip = int(nchi / ngroup);
    if (igroup < nchi % ngroup)
    {
        ++tmpnchip;
    }
    size_t memory = 0;
    for (int ik = 0; ik < nks; ++ik)
    {
        stowf.nchip[ik] = tmpnchip;
        stowf.chi0[ik].create(tmpnchip, ndim, false);
        memory += ndim * tmpnchip;
        if (seed_in >= 0)
        {
            for (int i = 0; i < stowf.chi0[ik].size; ++i)
            {
                const double phi = 2 * ModuleBase::PI * rand() / double(RAND_MAX);
                stowf.chi0[ik].c[i] = std::complex<double>(cos(phi), sin(phi)) / sqrt(double(nchi));
            }
        }
        else
        {
            for (int i = 0; i < stowf.chi0[ik].size; ++i)
            {
                if (rand() / double(RAND_MAX) < 0.5)
                {
                    stowf.chi0[ik].c[i] = -1.0 / sqrt(double(nchi));
                }
                else
                {
                    stowf.chi0[ik].c[i] = 1.0 / sqrt(double(nchi));
                }
            }
        }
    }
    ModuleBase::Memory::record("SDFT::chi0", memory * sizeof(std::complex<double>));
    stowf.nchip_max = tmpnchip;
}

void Update_Sto_Orbitals
(
    Stochastic_WF& stowf, 
    const int seed_in,
    const int nks   // GlobalC variables are deleted and become an input parameter
)
{
    const int nchi = INPUT.nbands_sto;

    for (int ik = 0; ik < nks; ++ik)
    {
        if (seed_in >= 0)
        {
            for (int i = 0; i < stowf.chi0[ik].size; ++i)
            {
                const double phi = 2 * ModuleBase::PI * rand() / double(RAND_MAX);
                stowf.chi0[ik].c[i] = std::complex<double>(cos(phi), sin(phi)) / sqrt(double(nchi));
            }
        }
        else
        {
            for (int i = 0; i < stowf.chi0[ik].size; ++i)
            {
                if (rand() / double(RAND_MAX) < 0.5)
                {
                    stowf.chi0[ik].c[i] = -1.0 / sqrt(double(nchi));
                }
                else
                {
                    stowf.chi0[ik].c[i] = 1.0 / sqrt(double(nchi));
                }
            }
        }
    }
}

#ifdef __MPI
void Init_Com_Orbitals
(
    Stochastic_WF& stowf,
    const int ndim  // GlobalC variables are deleted and become an input parameter
)
{
    const bool firstrankmore = false;
    int igroup;
    // former processor calculate more bands
    if (firstrankmore)
    {
        igroup = GlobalV::MY_STOGROUP;
    }
    // latter processor calculate more bands
    else
    {
        igroup = GlobalV::NSTOGROUP - GlobalV::MY_STOGROUP - 1;
    }
    const int nks = stowf.nks;
    const int ngroup = GlobalV::NSTOGROUP;
    const int n_in_pool = GlobalV::NPROC_IN_POOL;
    const int i_in_group = GlobalV::RANK_IN_STOGROUP;
    const int i_in_pool = GlobalV::RANK_IN_POOL;

    int* totnpw = new int[nks];
    for (int ik = 0; ik < nks; ++ik)
    {
        int* npwip = new int[n_in_pool];
        int* rec = new int[n_in_pool];
        int* displ = new int[n_in_pool];
        const int npw = stowf.ngk[ik];
        totnpw[ik] = 0;

        for (int i_in_p = 0; i_in_p < n_in_pool; ++i_in_p)
        {
            rec[i_in_p] = 1;
            displ[i_in_p] = i_in_p;
        }
        MPI_Allgatherv(&npw, 1, MPI_INT, npwip, rec, displ, MPI_INT, POOL_WORLD);
        for (int i_in_p = 0; i_in_p < n_in_pool; ++i_in_p)
        {
            totnpw[ik] += npwip[i_in_p];
        }

        int tmpnchip = int(totnpw[ik] / ngroup);
        if (igroup < totnpw[ik] % ngroup)
        {
            ++tmpnchip;
        }
        stowf.nchip[ik] = tmpnchip;
        stowf.chi0[ik].create(tmpnchip, ndim, true);
        stowf.nchip_max = std::max(tmpnchip,stowf.nchip_max);

        const int re = totnpw[ik] % ngroup;
        int ip = 0, ig0 = 0;
        // give value to orbitals in one parallel group one by one.
        for (int ichi = 0; ichi < tmpnchip; ++ichi)
        {
            int ig;
            if (igroup < re)
            {
                // It has more nchip.
                ig = igroup * tmpnchip + ichi - ig0;
            }
            else
            {
                // It has less nchip and should add re.
                ig = igroup * tmpnchip + re + ichi - ig0;
            }
            // get which ip stores this ig.
            while (ig >= npwip[ip])
            {
                ig -= npwip[ip];
                ig0 += npwip[ip];
                ++ip;
            }
            if (i_in_pool == ip)
            {
                stowf.chi0[ik](ichi, ig) = 1;
            }
        }

        delete[] npwip;
        delete[] rec;
        delete[] displ;
    }
    delete[] totnpw;
}
#else
void Init_Com_Orbitals
(
    Stochastic_WF& stowf, 
    const int ndim // GlobalC variables are deleted and become an input parameter
)
{
    for (int ik = 0; ik < stowf.nks; ++ik)
    {
        const int npw = stowf.ngk[ik];
        stowf.nchip[ik] = ndim;
        stowf.chi0[ik].create(stowf.nchip[ik], ndim, true);
        stowf.nchip_max = ndim;
        for (int ichi = 0; ichi < npw; ++ichi)
        {
            stowf.chi0[ik](ichi, ichi) = 1;
        }
    }
}
#endif
