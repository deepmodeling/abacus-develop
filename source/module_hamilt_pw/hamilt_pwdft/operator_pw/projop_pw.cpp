#include "projop_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/radial_proj.h"
#include "module_basis/module_nao/projgen.h"
#include "module_basis/module_nao/atomic_radials.h"
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
        // the following two variables will actually be deprecated, in Python it can be "_" but not in C++
        std::string elem = "__invalid__";
        double ecut = -0.1;
        std::ifstream ifs(forb);
        AtomicRadials::read_abacus_orb(ifs, elem, ecut, nr, dr, nzeta, radials);
#ifdef __DEBUG
        assert(elem != "__invalid__");
        assert(ecut != -0.1);
#endif 
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
                                      const int& nq)
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
        
        rp.sbtft(q,             // std::vector<ModuleBase::Vector3<double>>
                 proj_q_tab_,   // std::vector<std::complex<double>>
                 'r',           // char, 'r' for <G+k|p>, 'l' for <p|G+k>
                 omega,         // double
                 tpiba          // double
                 );
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
        int nr = -1;
        double dr = -1.0;
        bool padding = false;
        for(int it = 0; it < ntype; ++it)
        {
            std::vector<std::vector<double>> temp_;
            int lmax = -1;
            int nr_ = -1;
            std::vector<int> nzeta;
            read_abacus_orb(ucell_in->orbital_fn[it], lmax, nzeta, nr_, dr, temp_);
#ifdef __DEBUG
            assert(lmax != -1);
            assert(nr_ != -1);
            assert(dr != -1.0);
#endif
            padding = padding || (nr != -1 && nr != nr_);
            nr = std::max(nr, nr_);
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
        if(padding)
        {
            std::for_each(projs.begin(), projs.end(), [nr](std::vector<double>& proj) { proj.resize(nr, 0.0); });
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

        RadialProjection::RadialProjector rp;
        rp._build_sbt_tab(rgrid,     // std::vector<double>
                          projs,     // std::vector<std::vector<double>>
                          l_hubbard, // std::vector<int>
                          10000,     // int
                          0.01);     // double
        
        rp.sbtft(q,             // std::vector<ModuleBase::Vector3<double>>
                 proj_q_tab_,   // std::vector<std::complex<double>>
                 'r',           // char, 'r' for <G+k|p>, 'l' for <p|G+k>
                 omega,         // double
                 tpiba          // double
                 );
    }
}

// I am sorry but what does becp mean?...
#include "module_base/blas_connector.h"
void cal_becp(const int ik,
              Structure_Factor& sf,
              const ModulePW::PW_Basis_K& pw_basis,
              const psi::Psi<std::complex<double>, base_device::DEVICE_CPU>& psi,
              const std::vector<std::vector<int>>& it2ia,
              const std::vector<std::vector<int>>& it2iproj,
              const std::vector<double>& rgrid,
              const std::vector<std::vector<double>>& projs,
              const std::vector<int>& iproj2l,
              std::vector<std::complex<double>>& becp,
              const double& omega,
              const double& tpiba,
              const int nq = 10000,
              const double dq = 0.01)
{
    // For given ik, can have set of q
    const int npw = pw_basis.npwk[ik];
    std::vector<ModuleBase::Vector3<double>> q(npw);
    for(int ig = 0; ig < npw; ++ig)
    {
        q[ig] = pw_basis.getgpluskcar(ik, ig);
    }

    RadialProjection::RadialProjector rp;
    std::vector<int> irow2it;
    std::vector<int> irow2iproj;
    std::vector<int> irow2m;
    std::map<std::tuple<int, int, int, int>, int> itiaiprojm2irow;
    RadialProjection::RadialProjector::_build_backward_map(it2iproj, iproj2l, irow2it, irow2iproj, irow2m);
    RadialProjection::RadialProjector::_build_forward_map(it2ia, it2iproj, iproj2l, itiaiprojm2irow);
    const int nrow = irow2it.size();
    std::vector<std::complex<double>> tab_(nrow*npw);
    rp._build_sbt_tab(rgrid, projs, iproj2l, nq, dq);
    rp.sbtft(q, tab_, 'r', omega, tpiba);
    q.clear();
    q.shrink_to_fit(); // release memory

    std::vector<int> na(it2ia.size());
    for(int it = 0; it < it2ia.size(); ++it)
    {
        na[it] = it2ia[it].size();
    }
    // make_atomic
    const int nrow_out = itiaiprojm2irow.size();
    std::vector<std::complex<double>> tab_atomic_(nrow_out*npw);
    tab_atomic_.resize(nrow_out*npw);
    for(int irow = 0; irow < nrow; ++irow)
    {
        const int it = irow2it[irow];
        const int iproj = irow2iproj[irow];
        const int m = irow2m[irow];
        for(int ia = 0; ia < na[it]; ++ia)
        {
            // why Structure_Factor needs the FULL pw_basis???
            std::complex<double>* sk = sf.get_sk(ik, it, ia, &pw_basis);
            const int irow_out = itiaiprojm2irow.at(std::make_tuple(it, ia, iproj, m));
            for(int ig = 0; ig < npw; ++ig)
            {
                tab_atomic_[irow_out*npw + ig] = sk[ig]*tab_[irow*npw + ig];
            }
        }
    }
    tab_.clear();
    tab_.shrink_to_fit(); // release memory

    // cal_becp
    const int nbands = psi.get_nbands();
    const char transa = 'N';
    const char transb = 'N';
    const int one = 1;
    const int lda = nrow_out;
    const int ldb = npw;
    const int ldc = nrow_out;
    const std::complex<double> alpha = 1.0;
    const std::complex<double> beta = 0.0;

    becp.resize(nbands*nrow_out);
    BlasConnector::gemm(transa,                 // const char transa
                        transb,                 // const char transb
                        nrow_out,               // const int m
                        nbands,                 // const int n
                        npw,                    // const int k
                        alpha,                  // const std::complex<double> alpha
                        tab_atomic_.data(),     // const std::complex<double>* a
                        lda,                    // const int lda
                        psi.get_pointer(),      // const std::complex<double>* b
                        ldb,                    // const int ldb
                        beta,                   // const std::complex<double> beta
                        becp.data(),            // std::complex<double>* c
                        ldc);                   // const int ldc
    tab_atomic_.clear();
    tab_atomic_.shrink_to_fit(); // release memory
}