// wenfei 2022-1-5
// This file contains constructor and destructor of the class LCAO_deepks,
#include "module_parameter/parameter.h"
// as well as subroutines for initializing and releasing relevant data structures

// Other than the constructor and the destructor, it contains 3 types of subroutines:
// 1. subroutines that are related to calculating descriptors:
//   - init : allocates some arrays
//   - init_index : records the index (inl)
// 2. subroutines that are related to V_delta:
//   - allocate_V_delta : allocates V_delta; if calculating force, it also allocates F_delta

#ifdef __MLALGO

#include "LCAO_deepks.h"
#include "deepks_iterate.h"
#include "source_pw/hamilt_pwdft/global.h"

// Constructor of the class
template <typename T>
LCAO_Deepks<T>::LCAO_Deepks()
{
    inl_index = new ModuleBase::IntArray[1];
    gedm = nullptr;
    this->phialpha.resize(1);
}

// Desctructor of the class
template <typename T>
LCAO_Deepks<T>::~LCAO_Deepks()
{
    delete[] inl_index;

    //=======1. to use deepks, pdm is required==========
    pdm.clear();
    //=======2. "deepks_scf" part==========
    // if (PARAM.inp.deepks_scf)
    if (gedm)
    {
        // delete gedm**
        for (int inl = 0; inl < this->inlmax; inl++)
        {
            delete[] gedm[inl];
        }
        delete[] gedm;
    }
}

template <typename T>
void LCAO_Deepks<T>::init(const LCAO_Orbitals& orb,
                          const int nat,
                          const int ntype,
                          const int nks,
                          const Parallel_Orbitals& pv_in,
                          std::vector<int> na,
                          std::ofstream& ofs)
{
    ModuleBase::TITLE("LCAO_Deepks", "init");
    ModuleBase::timer::tick("LCAO_Deepks", "init");

    ofs << " Initialize the descriptor index for DeePKS (lcao line)" << std::endl;

    const int lm = orb.get_lmax_d();
    const int nm = orb.get_nchimax_d();
    const int tot_inl_per_atom = orb.Alpha[0].getTotal_nchi();

    assert(lm >= 0);
    assert(nm >= 0);
    assert(tot_inl_per_atom >= 0);

    int tot_inl = tot_inl_per_atom * nat;

    if (PARAM.inp.deepks_equiv)
    {
        tot_inl = nat;
    }

    this->lmaxd = lm;
    this->nmaxd = nm;

    ofs << " lmax of descriptor = " << this->lmaxd << std::endl;
    ofs << " nmax of descriptor = " << nmaxd << std::endl;

    int pdm_size = 0;
    this->inlmax = tot_inl;
    this->pdm.resize(this->inlmax);

    // cal n(descriptor) per atom , related to Lmax, nchi(L) and m. (not total_nchi!)
    if (!PARAM.inp.deepks_equiv)
    {
        this->des_per_atom = 0; // mohan add 2021-04-21
        for (int l = 0; l <= this->lmaxd; l++)
        {
            this->des_per_atom += orb.Alpha[0].getNchi(l) * (2 * l + 1);
        }
        this->n_descriptor = nat * this->des_per_atom;

        this->init_index(ntype, nat, na, tot_inl, orb, ofs);
    }

    if (!PARAM.inp.deepks_equiv)
    {
        ofs << " total basis (all atoms) for descriptor = " << std::endl;

        // init pdm
        for (int inl = 0; inl < this->inlmax; inl++)
        {
            int nm = 2 * inl2l[inl] + 1;
            pdm_size += nm * nm;
            this->pdm[inl] = torch::zeros({nm, nm}, torch::kFloat64);
        }
    }
    else
    {
        for (int il = 0; il < this->lmaxd + 1; il++)
        {
            pdm_size += (2 * il + 1) * orb.Alpha[0].getNchi(il);
        }
        pdm_size = pdm_size * pdm_size;
        this->des_per_atom = pdm_size;
        ofs << " Equivariant version, size of pdm matrices : " << pdm_size << std::endl;
        for (int iat = 0; iat < nat; iat++)
        {
            this->pdm[iat] = torch::zeros({pdm_size}, torch::kFloat64);
        }
    }

    this->pv = &pv_in;

    ModuleBase::timer::tick("LCAO_Deepks", "init");
    return;
}

template <typename T>
void LCAO_Deepks<T>::init_index(const int ntype,
                                const int nat,
                                std::vector<int> na,
                                const int Total_nchi,
                                const LCAO_Orbitals& orb,
                                std::ofstream& ofs)
{
    delete[] this->inl_index;
    this->inl_index = new ModuleBase::IntArray[ntype];
    this->inl2l.resize(this->inlmax, 0);

    int inl = 0;
    int alpha = 0;
    for (int it = 0; it < ntype; it++)
    {
        this->inl_index[it].create(na[it], this->lmaxd + 1, this->nmaxd);

        ofs << " Type " << it + 1 << " number_of_atoms " << na[it] << std::endl;

        for (int ia = 0; ia < na[it]; ia++)
        {
            // alpha
            for (int l = 0; l < this->lmaxd + 1; l++)
            {
                for (int n = 0; n < orb.Alpha[0].getNchi(l); n++)
                {
                    this->inl_index[it](ia, l, n) = inl;
                    this->inl2l[inl] = l;
                    inl++;
                }
            }
        } // end ia
    }     // end it
    assert(Total_nchi == inl);
    ofs << " descriptors_per_atom " << this->des_per_atom << std::endl;
    ofs << " total_descriptors " << this->n_descriptor << std::endl;
    return;
}

template <typename T>
void LCAO_Deepks<T>::allocate_V_delta(const int nat, const int nks)
{
    ModuleBase::TITLE("LCAO_Deepks", "allocate_V_delta");
    ModuleBase::timer::tick("LCAO_Deepks", "allocate_V_delta");

    // initialize the H matrix V_delta
    V_delta.resize(nks);
    for (int ik = 0; ik < nks; ik++)
    {
        this->V_delta[ik].resize(pv->nloc);
        ModuleBase::GlobalFunc::ZEROS(this->V_delta[ik].data(), pv->nloc);
    }

    // init gedm**
    int pdm_size = 0;
    if (!PARAM.inp.deepks_equiv)
    {
        pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    }
    else
    {
        pdm_size = this->des_per_atom;
    }

    this->gedm = new double*[this->inlmax];
    for (int inl = 0; inl < this->inlmax; inl++)
    {
        this->gedm[inl] = new double[pdm_size];
        ModuleBase::GlobalFunc::ZEROS(this->gedm[inl], pdm_size);
    }

    ModuleBase::timer::tick("LCAO_Deepks", "allocate_V_delta");
    return;
}

template <typename T>
void LCAO_Deepks<T>::init_DMR(const UnitCell& ucell,
                              const LCAO_Orbitals& orb,
                              const Parallel_Orbitals& pv,
                              const Grid_Driver& GridD)
{
    this->dm_r = new hamilt::HContainer<double>(&pv);
    DeePKS_domain::iterate_ad2(
        ucell,
        GridD,
        orb,
        false, // no trace_alpha
        [&](const int iat,
            const ModuleBase::Vector3<double>& tau0,
            const int ibt1,
            const ModuleBase::Vector3<double>& tau1,
            const int start1,
            const int nw1_tot,
            ModuleBase::Vector3<int> dR1,
            const int ibt2,
            const ModuleBase::Vector3<double>& tau2,
            const int start2,
            const int nw2_tot,
            ModuleBase::Vector3<int> dR2) 
        {
            auto row_indexes = pv.get_indexes_row(ibt1);
            auto col_indexes = pv.get_indexes_col(ibt2);
            if (row_indexes.size() * col_indexes.size() == 0)
            {
                return; // to next loop
            }

            int dRx = 0;
            int dRy = 0;
            int dRz = 0;
            if (std::is_same<T, std::complex<double>>::value)
            {
                dRx = (dR1 - dR2).x;
                dRy = (dR1 - dR2).y;
                dRz = (dR1 - dR2).z;
            }
            hamilt::AtomPair<double> dm_pair(ibt1, ibt2, dRx, dRy, dRz, &pv);
            this->dm_r->insert_pair(dm_pair);
        }
    );
    this->dm_r->allocate(nullptr, true);
}

template <typename T>
void LCAO_Deepks<T>::dpks_cal_e_delta_band(const std::vector<std::vector<T>>& dm, const int nks)
{
    DeePKS_domain::cal_e_delta_band(dm, this->V_delta, nks, PARAM.inp.nspin, this->pv, this->e_delta_band);
}

template class LCAO_Deepks<double>;
template class LCAO_Deepks<std::complex<double>>;

#endif
