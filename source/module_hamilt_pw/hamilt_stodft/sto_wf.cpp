#include "sto_wf.h"

#include "module_base/memory.h"
#include "module_parameter/parameter.h"

#include <cassert>
#include <ctime>

//---------Temporary------------------------------------
#include "module_base/global_function.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//------------------------------------------------------

template <typename T, typename Device>
Stochastic_WF<T, Device>::Stochastic_WF()
{
}

template <typename T, typename Device>
Stochastic_WF<T, Device>::~Stochastic_WF()
{
    delete chi0_cpu;
    Device* ctx = {};
    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        delete chi0;
    }
    delete shchi;
    delete chiortho;
    delete[] nchip;
    delete[] chiallorder;
}

template <typename T, typename Device>
void Stochastic_WF<T, Device>::init(K_Vectors* p_kv, const int npwx_in)
{
    this->nks = p_kv->get_nks();
    this->ngk = p_kv->ngk;
    this->npwx = npwx_in;
    nchip = new int[nks];

    if (nks <= 0)
    {
        ModuleBase::WARNING_QUIT("Stochastic_WF", "nks <=0!");
    }
}

template <typename T, typename Device>
void Stochastic_WF<T, Device>::allocate_chiallorder(const int& norder)
{
    this->chiallorder = new psi::Psi<T, Device>[this->nks];
    for (int ik = 0; ik < this->nks; ++ik)
    {
        chiallorder[ik].resize(1, this->nchip[ik] * this->npwx, norder);
        setmem_complex_op()(chiallorder[ik].get_pointer(), 0, chiallorder[ik].size());
    }
}

template <typename T, typename Device>
void Stochastic_WF<T, Device>::clean_chiallorder()
{
    delete[] chiallorder;
    chiallorder = nullptr;
}

template <typename T, typename Device>
void Stochastic_WF<T, Device>::init_sto_orbitals(const int seed_in)
{
    if (seed_in == 0 || seed_in == -1)
    {
        srand((unsigned)time(nullptr) + GlobalV::MY_RANK * 10000); // GlobalV global variables are reserved
    }
    else
    {
        srand((unsigned)std::abs(seed_in) + (GlobalV::MY_STOGROUP * GlobalV::NPROC_IN_STOGROUP + GlobalV::RANK_IN_STOGROUP) * 10000);
    }

    this->allocate_chi0();
    this->update_sto_orbitals(seed_in);
}

template <typename T, typename Device>
void Stochastic_WF<T, Device>::allocate_chi0()
{
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
        igroup = PARAM.inp.bndpar - GlobalV::MY_STOGROUP - 1;
    }
    const int nchi = PARAM.inp.nbands_sto;
    const int npwx = this->npwx;
    const int nks = this->nks;
    const int ngroup = PARAM.inp.bndpar;
    if (ngroup <= 0)
    {
        ModuleBase::WARNING_QUIT("init_sto_orbitals", "ngroup <= 0!");
    }
    int tmpnchip = int(nchi / ngroup);
    if (igroup < nchi % ngroup)
    {
        ++tmpnchip;
    }

    this->nchip_max = tmpnchip;
    size_t size = this->nchip_max * npwx * nks;
    this->chi0_cpu = new psi::Psi<T>(nks, this->nchip_max, npwx, this->ngk, true);
    ModuleBase::Memory::record("SDFT::chi0_cpu", size * sizeof(T));

    for (int ik = 0; ik < nks; ++ik)
    {
        this->nchip[ik] = tmpnchip;
    }

    // allocate chi0
    Device* ctx = {};
    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        this->chi0 = new psi::Psi<T, Device>(nks, this->nchip_max, npwx, this->ngk, true);
    }
    else
    {
        this->chi0 = reinterpret_cast<psi::Psi<T, Device>*>(this->chi0_cpu);
    }
}

template <typename T, typename Device>
void Stochastic_WF<T, Device>::update_sto_orbitals(const int seed_in)
{
    const int nchi = PARAM.inp.nbands_sto;
    this->chi0_cpu->fix_k(0);
    if (seed_in >= 0)
    {
        for (int i = 0; i < this->chi0_cpu->size(); ++i)
        {
            const double phi = 2 * ModuleBase::PI * rand() / double(RAND_MAX);
            this->chi0_cpu->get_pointer()[i] = std::complex<double>(cos(phi), sin(phi)) / sqrt(double(nchi));
        }
    }
    else
    {
        for (int i = 0; i < this->chi0_cpu->size(); ++i)
        {
            if (rand() / double(RAND_MAX) < 0.5)
            {
                this->chi0_cpu->get_pointer()[i] = -1.0 / sqrt(double(nchi));
            }
            else
            {
                this->chi0_cpu->get_pointer()[i] = 1.0 / sqrt(double(nchi));
            }
        }
    }
    this->sync_chi0();
}

#ifdef __MPI
template <typename T, typename Device>
void Stochastic_WF<T, Device>::init_com_orbitals()
{
    const bool firstrankmore = false;
    const int npwx = this->npwx;
    const int nks = this->nks;
    int igroup;
    // former processor calculate more bands
    if (firstrankmore)
    {
        igroup = GlobalV::MY_STOGROUP;
    }
    // latter processor calculate more bands
    else
    {
        igroup = PARAM.inp.bndpar - GlobalV::MY_STOGROUP - 1;
    }
    const int ngroup = PARAM.inp.bndpar;
    const int n_in_pool = GlobalV::NPROC_IN_POOL;
    const int i_in_group = GlobalV::RANK_IN_STOGROUP;
    const int i_in_pool = GlobalV::RANK_IN_POOL;

    int* totnpw = new int[nks];
    for (int ik = 0; ik < nks; ++ik)
    {
        int* npwip = new int[n_in_pool];
        const int npw = this->ngk[ik];
        totnpw[ik] = 0;

        MPI_Allgather(&npw, 1, MPI_INT, npwip, 1, MPI_INT, POOL_WORLD);
        for (int i_in_p = 0; i_in_p < n_in_pool; ++i_in_p)
        {
            totnpw[ik] += npwip[i_in_p];
        }

        int tmpnchip = int(totnpw[ik] / ngroup);
        if (igroup < totnpw[ik] % ngroup)
        {
            ++tmpnchip;
        }
        this->nchip[ik] = tmpnchip;
        this->nchip_max = std::max(tmpnchip, this->nchip_max);
        delete[] npwip;
    }
    size_t size = this->nchip_max * npwx * nks;
    this->chi0_cpu = new psi::Psi<std::complex<double>>(nks, this->nchip_max, npwx, this->ngk, true);
    this->chi0_cpu->zero_out();
    ModuleBase::Memory::record("SDFT::chi0_cpu", size * sizeof(std::complex<double>));
    for (int ik = 0; ik < nks; ++ik)
    {
        int* npwip = new int[n_in_pool];
        const int npw = this->ngk[ik];
        MPI_Allgather(&npw, 1, MPI_INT, npwip, 1, MPI_INT, POOL_WORLD);
        const int re = totnpw[ik] % ngroup;
        int ip = 0, ig0 = 0;
        const int nchipk = this->nchip[ik];
        // give value to orbitals in one parallel group one by one.
        for (int ichi = 0; ichi < nchipk; ++ichi)
        {
            int ig;
            if (igroup < re)
            {
                // It has more nchip.
                ig = igroup * nchipk + ichi - ig0;
            }
            else
            {
                // It has less nchip and should add re.
                ig = igroup * nchipk + re + ichi - ig0;
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
                this->chi0_cpu->operator()(ik, ichi, ig) = 1;
            }
        }

        delete[] npwip;
    }
    delete[] totnpw;
    // allocate chi0
    Device* ctx = {};
    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        this->chi0 = new psi::Psi<T, Device>(nks, this->nchip_max, npwx, this->ngk, true);
    }
    else
    {
        this->chi0 = reinterpret_cast<psi::Psi<T, Device>*>(this->chi0_cpu);
    }
}
#else
template <typename T, typename Device>
void Stochastic_WF<T, Device>::init_com_orbitals()
{
    const int npwx = this->npwx;
    const int nks = this->nks;
    size_t size = this->nchip_max * npwx * nks;
    this->chi0_cpu = new psi::Psi<std::complex<double>>(nks, npwx, npwx, this->ngk, true);
    this->chi0_cpu->zero_out();
    ModuleBase::Memory::record("SDFT::chi0_cpu", size * sizeof(std::complex<double>));
    for (int ik = 0; ik < nks; ++ik)
    {
        const int npw = this->ngk[ik];
        this->nchip[ik] = npwx;
        this->nchip_max = npwx;
        for (int ichi = 0; ichi < npw; ++ichi)
        {
            this->chi0_cpu->operator()(ik, ichi, ichi) = 1;
        }
    }

    // allocate chi0
    Device* ctx = {};
    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        this->chi0 = new psi::Psi<T, Device>(nks, this->nchip_max, npwx, this->ngk, true);
    }
    else
    {
        this->chi0 = reinterpret_cast<psi::Psi<T, Device>*>(this->chi0_cpu);
    }
}
#endif
template <typename T, typename Device>
void Stochastic_WF<T, Device>::init_sto_orbitals_Ecut(const int seed_in,
                                                      const K_Vectors& kv,
                                                      const ModulePW::PW_Basis_K& wfcpw,
                                                      const int max_ecut)
{
    this->allocate_chi0();

    ModulePW::PW_Basis pwmax;
#ifdef __MPI
    pwmax.initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
    pwmax.initgrids(wfcpw.lat0, wfcpw.latvec, max_ecut);
    const int nx = pwmax.nx;
    const int ny = pwmax.ny;
    const int nz = pwmax.nz;
    const int nkstot = kv.get_nkstot();
    const int nks = kv.get_nks();
    const int nchitot = PARAM.inp.nbands_sto;
    bool* updown = new bool[nx * ny * nz];
    int* nrecv = new int[PARAM.inp.bndpar];
    const int nchiper = this->nchip[0];
#ifdef __MPI
    MPI_Allgather(&nchiper, 1, MPI_INT, nrecv, 1, MPI_INT, PARAPW_WORLD);
#endif
    int ichi_start = 0;
    for (int i = 0; i < GlobalV::MY_STOGROUP; ++i)
    {
        ichi_start += nrecv[i];
    }

    for (int ik = 0; ik < nks; ++ik)
    {
        const int iktot = kv.ik2iktot[ik];
        const int npw = wfcpw.npwk[ik];
        int* ig2ixyz = new int[npw];

        for (int ig = 0; ig < npw; ++ig)
        {
            ModuleBase::Vector3<double> gdirect = wfcpw.getgdirect(ik, ig);
            int ix = static_cast<int>(gdirect.x);
            int iy = static_cast<int>(gdirect.y);
            int iz = static_cast<int>(gdirect.z);
            ix = (ix + nx) % nx;
            iy = (iy + ny) % ny;
            iz = (iz + nz) % nz;
            ig2ixyz[ig] = ix * ny * nz + iy * nz + iz;
        }

        for (int ichi = 0; ichi < nchiper; ++ichi)
        {
            unsigned int seed = std::abs(seed_in) * (nkstot * nchitot) + iktot * nchitot + ichi_start + ichi;
            srand(seed);
            for (int i = 0; i < nx * ny * nz; ++i)
            {
                updown[i] = (rand() / double(RAND_MAX) < 0.5);
            }

            for (int ig = 0; ig < npw; ++ig)
            {
                if (updown[ig2ixyz[ig]])
                {
                    this->chi0_cpu->operator()(ik, ichi, ig) = -1.0 / sqrt(double(nchitot));
                }
                else
                {
                    this->chi0_cpu->operator()(ik, ichi, ig) = 1.0 / sqrt(double(nchitot));
                }
            }
        }
        delete[] ig2ixyz;
    }
    delete[] nrecv;
    delete[] updown;
}

template <typename T, typename Device>
void Stochastic_WF<T, Device>::sync_chi0()
{
    Device* ctx = {};
    if (base_device::get_device_type<Device>(ctx) == base_device::GpuDevice)
    {
        syncmem_h2d_op()(this->chi0->get_pointer(),
                         this->chi0_cpu->get_pointer(),
                         this->chi0_cpu->size());
    }
}

template class Stochastic_WF<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stochastic_WF<std::complex<double>, base_device::DEVICE_GPU>;
#endif
