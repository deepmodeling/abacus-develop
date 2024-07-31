#include "realspace_proj_op_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/radial_projection.h"
#include "module_basis/module_nao/projgen.h"
#include <cassert>
#include <numeric>
namespace hamilt
{
    template<typename T, typename Device>
    void DFTU<OperatorPW<T, Device>>::read_abacus_orb(const std::string& forb,
                                                      int& lmax,
                                                      std::vector<int>& nzeta,
                                                      int& nr,
                                                      double& dr,
                                                      std::vector<std::vector<double>>& radials)
    {
        // waiting for AtomicRadials::read_abacus_orb to be static or being overloaded a
        // static version. Then this implementation will be removed.
    }

    template<typename T, typename Device>
    DFTU<OperatorPW<T, Device>>::DFTU(const std::vector<int> isk,
                                      const std::vector<int>& l_hubbard,
                                      const std::vector<double>& u,
                                      const std::vector<double>& rgrid,
                                      const std::vector<std::vector<double>>& projs,
                                      const std::vector<int>& natom,
                                      const std::vector<ModuleBase::Vector3<double>*>& tau,
                                      const double& omega,
                                      const double& tpiba,
                                      const std::vector<ModuleBase::Vector3<double>>& q,
                                      const double& dq,
                                      const double& nq)
    {
        RadialProjection::RadialProjector rp;
        const int nr = rgrid.size();
#ifdef __DEBUG
        for(auto &proj: projs)
        {
            assert(proj.size() == nr);
        }
#endif

        rp._build_sbt_tab(rgrid,     // std::vector<double>
                          projs,     // std::vector<std::vector<double>>
                          l_hubbard, // std::vector<int>
                          nq,        // int
                          dq);       // double
        
        std::vector<std::complex<double>> out;
        rp.sbtft(q,     // std::vector<ModuleBase::Vector3<double>>
                 out,   // std::vector<std::complex<double>>
                 omega, // double
                 tpiba, // double
                 'r');  // char, 'r' for <G+k|p>, 'l' for <p|G+k>
    }

    template<typename T, typename Device>
    DFTU<OperatorPW<T, Device>>::DFTU(const std::vector<int>& isk,
                                      const UnitCell* ucell_in,
                                      const ModulePW::PW_Basis_K* pw_basis)
    {
        const double omega = ucell_in->omega;
        const double tpiba = ucell_in->tpiba;
        std::vector<ModuleBase::Vector3<double>> q;
        for(int ig = 0; ig < pw_basis->npw; ++ig)
        {   // the following line is incorrect, because I am not clear about what the isk is!
            q.push_back(pw_basis->getgpluskcar(isk[0], ig)); // not isk[0], but actually loop on isk
        }
        const int ntype = ucell_in->ntype;

        std::vector<int> l_hubbard(ntype);   // read from UnitCell or somewhere else
        std::vector<double> u(ntype);        // read from UnitCell or somewhere else
        std::vector<double> onsite_r(ntype); // read from UnitCell or somewhere else
        
        std::vector<std::vector<double>> projs(ntype);
        int nr;
        double dr;
        for(int it = 0; it < ntype; ++it)
        {
            std::vector<std::vector<double>> temp_;
            int lmax;
            std::vector<int> nzeta;
            read_abacus_orb(ucell_in->orbital_fn[it], lmax, nzeta, nr, dr, temp_);
            // then get the first one of given l in l_hubbard[it]
            int idx = 0;
            // accumulate nzeta up to l_hubbard[it] (exclusive)
            std::accumulate(nzeta.begin(), nzeta.begin() + l_hubbard[it], idx);
            std::vector<double> r(nr);
            std::iota(r.begin(), r.end(), 0);
            std::for_each(r.begin(), r.end(), [dr](double& r_i) { r_i *= dr; });

            smoothgen(nr, r.data(), temp_[idx].data(), onsite_r[it], projs[it]); 
            // but... does this always hold??? forever only one proj for one type???
        }
        std::vector<int> natom(ntype);
        std::vector<ModuleBase::Vector3<double>*> tau(ucell_in->nat);
        int iat = 0;
        for(int it = 0; it < ntype; ++it)
        {
            natom[it] = ucell_in->atoms[it].na;
            for(int ia = 0; ia < natom[it]; ++ia)
            {
                tau[iat] = &ucell_in->atoms[it].tau[ia];
                ++iat;
            }
        }

        std::vector<double> rgrid(nr);
        std::iota(rgrid.begin(), rgrid.end(), 0);
        std::for_each(rgrid.begin(), rgrid.end(), [dr](double& r_i) { r_i *= dr; });

        DFTU<OperatorPW<T, Device>>(isk,        // what is this?
                                    l_hubbard,  // std::vector<int>
                                    u,          // std::vector<double>
                                    rgrid,      // std::vector<double>
                                    projs,      // std::vector<std::vector<double>>
                                    natom,      // std::vector<int>
                                    tau,        // std::vector<ModuleBase::Vector3<double>*>
                                    omega,      // double
                                    tpiba,      // double
                                    q);         // std::vector<ModuleBase::Vector3<double>>
    }
}