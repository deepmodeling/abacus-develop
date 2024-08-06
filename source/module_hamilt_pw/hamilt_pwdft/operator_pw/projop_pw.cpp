#include <cassert>
#include <numeric>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>
#include <complex>

#include "module_hamilt_pw/hamilt_pwdft/radial_proj.h"
#include "module_basis/module_nao/projgen.h"
#include "module_basis/module_nao/atomic_radials.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/unitcell.h"
#include "module_base/blas_connector.h"

/**
 * @brief initialize the radial projector for real-space projection involving operators
 * 
 * @param orbital_dir You know what it is
 * @param ucell You know what it is
 * @param nproj # of projectors for each type defined in UnitCell, can be zero
 * @param lproj angular momentum for each projector
 * @param iproj index of zeta function that each projector generated from
 * @param onsite_r onsite-radius for all valid projectors
 * @param rgrid [out] the radial grid shared by all projectors
 * @param projs [out] projectors indexed by `iproj`
 * @param it2iproj [out] for each type, the projector index (across all types)
 */
void init_proj(const std::string& orbital_dir,
               const UnitCell& ucell,
               const std::vector<int>& nproj,           // for each type, the number of projectors
               const std::vector<int>& lproj,           // angular momentum of projectors within the type (l of zeta function)
               const std::vector<int>& iproj,           // index of projectors within the type (izeta)
               const std::vector<double>& onsite_r,     // for each projector, the "onsite_radius"
               std::vector<double>& rgrid,              // the radial grid shared by all projectors
               std::vector<std::vector<double>>& projs, // projectors indexed by `iproj`
               std::vector<std::vector<int>>& it2iproj) // for each type, the projector index (across all types)
{
    // extract the information from ucell
    const int ntype = nproj.size();
    assert(ucell.ntype == ntype); // nproj can have 0 for some types
    int nproj_tot = 0;
    std::accumulate(nproj.begin(), nproj.end(), nproj_tot);
    assert(nproj_tot == lproj.size());
    assert(nproj_tot == iproj.size());
    assert(nproj_tot == onsite_r.size());
    projs.resize(nproj_tot);

    int idx = 0;
    int nr = -1;
    double dr = -1.0;
    for(int it = 0; it < ntype; ++it)
    {
        const int nproj_it = nproj[it];
        it2iproj[it].resize(nproj_it);
        std::ifstream ifs(orbital_dir + "/" + ucell.orbital_fn[it]);
        std::string elem;
        double ecut;
        int nr_;
        double dr_;
        std::vector<int> nzeta; // number of radials for each l
        std::vector<std::vector<double>> radials; // radials arranged in serial
        AtomicRadials::read_abacus_orb(ifs, elem, ecut, nr_, dr_, nzeta, radials);
        nr = std::max(nr, nr_); // the maximal nr
        assert(dr == -1.0 || dr == dr_); // the dr should be the same for all types
        dr = (dr == -1.0) ? dr_ : dr;
        for(int ip = 0; ip < nproj_it; ++ip)
        {
            int l = lproj[idx];
            int izeta = iproj[idx];
            int irad = 0;
            std::accumulate(nzeta.begin(), nzeta.begin() + l, irad);
            irad += izeta;
            std::vector<double> temp = radials[irad];
            smoothgen(nr, rgrid.data(), temp.data(), onsite_r[idx], projs[idx]);
            it2iproj[it][ip] = idx;
            ++idx;
        }
    }
    // do zero padding
    if(nr != -1)
    {
        std::for_each(projs.begin(), projs.end(), [nr](std::vector<double>& proj) { proj.resize(nr, 0.0); });
    }
    // generate the rgrid
    rgrid.resize(nr);
    std::iota(rgrid.begin(), rgrid.end(), 0);
    std::for_each(rgrid.begin(), rgrid.end(), [dr](double& r_i) { r_i *= dr; });
}

// I am sorry but what does becp mean?...
void cal_becp(const std::vector<std::vector<int>>& it2ia,       // level0: for given type `it`, the atom indices `ia`
              const std::vector<std::vector<int>>& it2iproj,    // level0: for given type `it`, the proj indices `iproj`
              const std::vector<double>& rgrid,                 // level0: the radial grid shared by all projectors
              const std::vector<std::vector<double>>& projs,    // level0: projectors indexed by `iproj`
              const std::vector<int>& iproj2l,                  // level0: for given proj index `iproj`, the angular momentum `l`
              const int nq,                                     // level0: GlobalV::NQX
              const double& dq,                                 // level0: GlobalV::DQ
              const int ik,                                     // level1: the k-point index
              const ModulePW::PW_Basis_K& pw_basis,             // level1: the plane wave basis, need ik
              const double& omega,                              // level1: the cell volume
              const double& tpiba,                              // level1: 2*pi/lat0
              Structure_Factor& sf,                             // level2: the structure factor calculator
              const psi::Psi<std::complex<double>, base_device::DEVICE_CPU>& psi,
              std::vector<std::complex<double>>& becp
              )
{
    // STAGE 0 - making the interpolation table
    // CACHE 0 - if cache the irow2it, irow2iproj, irow2m, itiaiprojm2irow, <G+k|p> can be reused for 
    //           SCF, RELAX and CELL-RELAX calculation
    // [in] rgrid, projs, iproj2l, it2ia, it2iproj, nq, dq
    RadialProjection::RadialProjector rp;
    std::vector<int> irow2it;
    std::vector<int> irow2iproj;
    std::vector<int> irow2m;
    std::map<std::tuple<int, int, int, int>, int> itiaiprojm2irow;
    RadialProjection::RadialProjector::_build_backward_map(it2iproj, iproj2l, irow2it, irow2iproj, irow2m);
    RadialProjection::RadialProjector::_build_forward_map(it2ia, it2iproj, iproj2l, itiaiprojm2irow);
    rp._build_sbt_tab(rgrid, projs, iproj2l, nq, dq);


    // STAGE 1 - calculate the <G+k|p> for the given G+k vector
    // CACHE 1 - if cache the tab_, <G+k|p> can be reused for SCF and RELAX calculation
    // [in] pw_basis, ik, omega, tpiba, irow2it
    const int npw = pw_basis.npwk[ik];
    std::vector<ModuleBase::Vector3<double>> q(npw);
    for(int ig = 0; ig < npw; ++ig)
    {
        q[ig] = pw_basis.getgpluskcar(ik, ig); // get the G+k vector, G+k will change during CELL-RELAX
    }
    const int nrow = irow2it.size();
    std::vector<std::complex<double>> tab_(nrow*npw);
    rp.sbtft(q, tab_, 'l', omega, tpiba); // l: <p|G+k>, r: <G+k|p>
    q.clear();
    q.shrink_to_fit(); // release memory


    // STAGE 2 - make_atomic: multiply e^iqtau and extend the <G+k|p> to <G+k|pi> for each atom
    // CACHE 2 - if cache the tab_atomic_, <G+k|p> can be reused for SCF calculation
    // [in] it2ia, itiaiprojm2irow, tab_, npw, sf
    std::vector<int> na(it2ia.size());
    for(int it = 0; it < it2ia.size(); ++it)
    {
        na[it] = it2ia[it].size();
    }
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


    // STAGE 3 - cal_becp
    // CACHE 3 - it is no use to cache becp, it will change in each SCF iteration
    // [in] psi, tab_atomic_, npw, becp, ik
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
    psi.fix_k(ik);
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
