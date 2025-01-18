#include <cassert>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <map>
#include <tuple>
#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"

#include "module_base/projgen.h"
#include "module_base/blas_connector.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#ifdef __MPI
#include "module_base/parallel_reduce.h"
#include "module_base/parallel_common.h"
#endif
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_base/formatter.h"



/**
 * ===============================================================================================
 * 
 *                                          README
 * 
 * ===============================================================================================
 * 
 * This is a code demo for illustrating how to use unified radial projection in implementation of
 * Operators involving local radial projectors on PW-expanded wavefunctions.
 * 
 * Example usage:
 * ```c++
 * // select the range of atoms that impose the operator in std::vector<std::vector<int>> it2ia like
 * // it2ia[it] = {ia1, ia2, ...} for each type
 * // if all atoms in present kind is "selected", just set it2ia[it].resize(na) and call 
 * // std::iota(it2ia[it].begin(), it2ia[it].end(), 0)
 * 
 * std::vector<std::vector<int>> it2ia; // as if we have given its value...
 * 
 * // you should have the `orbital_dir` as the directory containing the orbital files, then those
 * // will be read by a static function `AtomicRadials::read_abacus_orb` to get the radial orbitals
 * 
 * // call `init_proj` to initialize the radial projector, this function only needs to be called
 * // once during the runtime.
 * // its input... 
 * // the `nproj`, is for specifying number of projectors of each atom type, can be zero,
 * // but cannot be the value larger than the number of zeta functions for the given angular momentum.
 * // the `lproj` is the angular momentum of the projectors, and `iproj` is the index of zeta function
 * // that each projector generated from.
 * // the `lproj` along with `iproj` can enable radial projectors in any number developer wants.
 * 
 * // the `onsite_r` is the onsite-radius for all valid projectors, it is used to generate the new
 * // radial function that more localized than the original one, which is expected to have enhanced
 * // projection efficiency.
 * 
 * std::vector<double> rgrid;
 * std::vector<std::vector<double>> projs;
 * std::vector<std::vector<int>> it2iproj;
 * init_proj(orbital_dir, ucell, nproj, lproj, iproj, onsite_r, rgrid, projs, it2iproj);
 * 
 * // then call the function `cal_becp` to calculate the becp. HOWEVER, there are quantities that
 * // can be calculated in advance and reused in the following calculations. Please see the function
 * // implementation, especially the comments about CACHE 0, CACHE 1, CACHE 2..., etc.
 * 
 * // the input param of `cal_becp`...
 * // the `it2ia` has been explained above
 * // the `it2iproj` is the output of function `init_proj`, so you do not need to worry about it
 * // the `rgrid` and `projs` are also the output of function `init_proj`
 * // the `lproj` is the angular momentum for each projector, actually you have used it in `init_proj`, it
 * // is the same as `lproj`
 * // the `nq` is the number of G+k vectors, typically it is always GlobalV::NQX
 * // the `dq` is the step size of G+k vectors, typically it is always GlobalV::DQ
 * // the `ik` is the k-point index
 * // the `pw_basis` is the plane wave basis, need ik
 * // the `omega` is the cell volume
 * // the `tpiba` is 2*pi/lat0
 * // the `sf` is the structure factor calculator
 * // the `psi` is the wavefunction
 * // the `becp` is the output of the function, it is the becp
 * cal_becp(it2ia, it2iproj, rgrid, projs, lproj, nq, dq, ik, pw_basis, omega, tpiba, sf, psi, becp);
 * 
 * // About parallelization, presently, the function `AtomicRadials::read_abacus_orb` is actually parallelized
 * // by MPI, so after the reading of orbital, actually all processors have the same data. Therefore it is not
 * // needed to call functions like `Parallel_Reduce` or `Parallel_Bcast` to synchronize the data.
 * // However, what is strikingly memory-consuming is the table `tab_atomic_`. Performance optimization will
 * // be needed if the memory is not enough.
 */

template<typename T, typename Device>
projectors::OnsiteProjector<T, Device>* projectors::OnsiteProjector<T, Device>::get_instance()
{
    static projectors::OnsiteProjector<T, Device> instance;
    return &instance;
}

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::init(const std::string& orbital_dir,
                                                  const UnitCell* ucell_in,
                                                  const psi::Psi<std::complex<T>, Device>& psi,
                                                  const K_Vectors& kv,
                                                  const ModulePW::PW_Basis_K& pw_basis,             // level1: the plane wave basis, need ik
                                                  Structure_Factor& sf,                              // level2: the structure factor calculator
                                                  const double onsite_radius,
                                                  const int nq,
                                                  const double dq,
                                                  const ModuleBase::matrix& wg,
                                                  const ModuleBase::matrix& ekb)
{
    this->device = base_device::get_device_type<Device>(this->ctx);
    if(!this->initialed)
    {
        this->ucell = ucell_in;
        this->ntype = ucell_in->ntype;

        this->pw_basis_ = &pw_basis;
        this->sf_ = &sf;

        std::vector<std::string> orb_files(ntype);
        std::vector<int> nproj(ntype);
        int sum_nproj = 0;
        for(int it=0;it<ntype;++it)
        {
            orb_files[it] = ucell->orbital_fn[it];
            nproj[it] = ucell->atoms[it].nwl;
            sum_nproj += nproj[it];
        }
        this->lproj.resize(sum_nproj);
        int index = 0;
        for(int it=0;it<ntype;++it)
        {
            for(int il=0;il<nproj[it];++il)
            {
                this->lproj[index++] = il;
            }
        }
        std::vector<int> iproj(sum_nproj, 0);
        std::vector<double> onsite_r(sum_nproj, onsite_radius);

        this->it2ia.resize(this->ntype);
        this->iat_nh.resize(this->ucell->nat);
        int iat = 0;
        for(int it = 0; it < it2ia.size(); it++)
        {
            it2ia[it].resize(this->ucell->atoms[it].na);
            std::iota(it2ia[it].begin(), it2ia[it].end(), 0);
            for(int ia = 0; ia < it2ia[it].size(); ia++)
            {
                iat_nh[iat++] = nproj[it] * nproj[it];
            }
        }

        this->init_proj(PARAM.inp.orbital_dir, 
                        orb_files, 
                        nproj, 
                        lproj, 
                        iproj, 
                        onsite_r);

        ModuleBase::timer::tick("OnsiteProj", "cubspl_tabulate");
        // STAGE 0 - making the interpolation table
        // CACHE 0 - if cache the irow2it, irow2iproj, irow2m, itiaiprojm2irow, <G+k|p> can be reused for 
        //           SCF, RELAX and CELL-RELAX calculation
        // [in] rgrid, projs, lproj, it2ia, it2iproj, nq, dq
        RadialProjection::RadialProjector::_build_backward_map(it2iproj, lproj, irow2it_, irow2iproj_, irow2m_);
        RadialProjection::RadialProjector::_build_forward_map(it2ia, it2iproj, lproj, itiaiprojm2irow_);
        //rp_._build_sbt_tab(rgrid, projs, lproj, nq, dq);
        rp_._build_sbt_tab(nproj, rgrid, projs, lproj, nq, dq, ucell_in->omega, psi.npol, tab, nhtol);
        // For being compatible with present cal_force and cal_stress framework  
        // uncomment the following code block if you want to use the Onsite_Proj_tools
        if(this->tab_atomic_ == nullptr)
        {
            this->tot_nproj = itiaiprojm2irow_.size();
            this->npwx_ = this->pw_basis_->npwk_max;
            this->size_vproj = this->tot_nproj * this->npwx_;
            resmem_complex_op()(this->tab_atomic_, this->size_vproj, "OnsiteP::tab_atomic_");
        }

        delete this->fs_tools; // it is okay to delete nullptr
        this->fs_tools = new hamilt::Onsite_Proj_tools<T, Device>(
            nproj, lproj, tab, nhtol, this->tab_atomic_, ucell_in, &psi, &kv, &pw_basis, &sf, wg, ekb);      
        
        ModuleBase::timer::tick("OnsiteProj", "cubspl_tabulate");

        this->initialed = true;
    }
}

template<typename T, typename Device>
projectors::OnsiteProjector<T, Device>::~OnsiteProjector()
{
    //delete[] becp;
    delete fs_tools;
    delmem_complex_op()(this->tab_atomic_);
    if(this->device == base_device::GpuDevice)
    {
        delmem_complex_h_op()(this->h_becp);
    }
    delmem_complex_op()(this->becp);

}


template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::init_proj(const std::string& orbital_dir,
                     const std::vector<std::string>& orb_files,
                     const std::vector<int>& nproj,           // for each type, the number of projectors
                     const std::vector<int>& lproj,           // angular momentum of projectors within the type (l of zeta function)
                     const std::vector<int>& iproj,           // index of projectors within the type (izeta)
                     const std::vector<double>& onsite_r) 
{
    // extract the information from ucell
    const int ntype = nproj.size();
    assert(ntype == orb_files.size());
    this->it2iproj.resize(ntype);

    int nproj_tot = 0;
    nproj_tot = std::accumulate(nproj.begin(), nproj.end(), nproj_tot, std::plus<int>());
    assert(nproj_tot == lproj.size());
    assert(nproj_tot == iproj.size());
    assert(nproj_tot == onsite_r.size());
    this->projs.resize(nproj_tot);

    int idx = 0;
    int nr = -1;
    double dr = -1.0;
    for(int it = 0; it < ntype; ++it)
    {
        const int nproj_it = nproj[it];
        this->it2iproj[it].resize(nproj_it);
        if(nproj_it == 0)
        {
            std::cout << "BECP_PW >> No projectors defined for type " << it << std::endl;
            continue;
        }
        std::ifstream ifs(orbital_dir + orb_files[it]);
        std::string elem = "";
        double ecut = -1.0;
        int nr_ = -1;
        double dr_ = -1.0;
        std::vector<int> nzeta; // number of radials for each l
        std::vector<std::vector<double>> radials; // radials arranged in serial
        this->read_abacus_orb(ifs, elem, ecut, nr_, dr_, nzeta, radials);
#ifdef __DEBUG
        assert(elem != "");
        assert(ecut != -1.0);
        assert(nr_ != -1);
        assert(dr_ != -1.0);
#endif
        nr = std::max(nr, nr_); // the maximal nr
        assert(dr == -1.0 || dr == dr_); // the dr should be the same for all types
        dr = (dr == -1.0) ? dr_ : dr;
        for(int ip = 0; ip < nproj_it; ++ip)
        {
            int l = lproj[idx];
            int izeta = iproj[idx];
            int irad = 0;
            irad = std::accumulate(nzeta.begin(), nzeta.begin() + l, irad);
            irad += izeta;
            std::vector<double> temp = radials[irad];
            rgrid.resize(nr);
            std::iota(rgrid.begin(), rgrid.end(), 0);
            std::for_each(rgrid.begin(), rgrid.end(), [dr](double& r_i) { r_i *= dr; });
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
    this->rgrid.resize(nr);
    std::iota(rgrid.begin(), rgrid.end(), 0);
    std::for_each(rgrid.begin(), rgrid.end(), [dr](double& r_i) { r_i *= dr; });
}

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::tabulate_atomic(const int ik, const char grad)
{
    ModuleBase::timer::tick("OnsiteProj", "tabulate_atomic");
    // assert(grad == 'n' || grad == 'x' || grad == 'y' || grad == 'z');
    // grad = 'n' means no gradient, grad = 'x' means gradient along x, etc.

    // STAGE 1 - calculate the <G+k|p> for the given G+k vector
    // CACHE 1 - if cache the tab_, <G+k|p> can be reused for SCF and RELAX calculation
    // [in] pw_basis, ik, omega, tpiba, irow2it
    this->ik_ = ik;
    this->npw_ = pw_basis_->npwk[ik];
    this->npwx_ = pw_basis_->npwk_max;
    // std::vector<ModuleBase::Vector3<double>> q(this->npw_);
    // for(int ig = 0; ig < this->npw_; ++ig)
    // {
    //     q[ig] = pw_basis_->getgpluskcar(ik, ig); // get the G+k vector, G+k will change during CELL-RELAX
    // }
    // const int nrow = irow2it_.size();
    // std::vector<std::complex<double>> tab_(nrow*this->npw_);
    // // convention used here: 'l': <p|G+k>, 'r': <G+k|p>
    // // denote q=G+k, <r|q> = exp(iqr), the routine Fourier Transform written as F(q) = <q|f>
    // rp_.sbtft(q, tab_, 'l', this->ucell->omega, this->ucell->tpiba);
    // // what is calculated is <p|q> here

    // STAGE 2 - make_atomic: multiply e^iqtau and extend the <G+k|p> to <G+k|pi> for each atom
    // CACHE 2 - if cache the tab_atomic_, <G+k|p> can be reused for SCF calculation
    // [in] it2ia, itiaiprojm2irow, tab_, npw, sf
    // for(int irow = 0; irow < nrow; ++irow)
    // {
    //     const int it = irow2it_[irow];
    //     const int iproj = irow2iproj_[irow];
    //     const int m = irow2m_[irow];
    //     for(int ia = 0; ia < na[it]; ++ia)
    //     {
    //         // why Structure_Factor needs the FULL pw_basis???
    //         std::complex<double>* sk = this->sf_->get_sk(ik, it, ia, pw_basis_); // exp(-iqtau)
    //         // Note: idea on extending the param list of get_sk
    //         // the get_sk should have an extra param 'grad' to calculate the gradient of S(q), which
    //         // is actually very simple to be
    //         // d(S(q))/dq = -i S(q) * tau, for one direction it is just -i S(q) * tau_x (if x is the direction)
    //         const int irow_out = itiaiprojm2irow_.at(std::make_tuple(it, ia, iproj, m));
    //         for(int ig = 0; ig < this->npw_; ++ig)
    //         {
    //             std::complex<double> deriv = (grad == 'n')? 1.0: ModuleBase::NEG_IMAG_UNIT; // because sk is exp(-iqtau)
    //             deriv = (grad == 'n')? 1.0: (grad == 'x')? deriv * q[ig].x: (grad == 'y')? deriv * q[ig].y: deriv * q[ig].z;
    //             // there must be something twisted in ABACUS
    //             // because the tab_ is <p|G+k>, but the sk is exp(-iqtau). How can it get the 
    //             // correct result?
    //             this->tab_atomic_[irow_out*this->npw_ + ig] = sk[ig] * tab_[irow*this->npw_ + ig] * deriv;
    //         }
    //         delete[] sk;
    //     }
    // }
    // q.clear();
    // q.shrink_to_fit();    // release memory
    // tab_.clear();
    // tab_.shrink_to_fit(); // release memory
    ModuleBase::timer::tick("OnsiteProj", "tabulate_atomic");
}

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::overlap_proj_psi( 
                    const int npm,
                    const std::complex<double>* ppsi
                    )
{
    ModuleBase::timer::tick("OnsiteProj", "overlap");
    // STAGE 3 - cal_becp
    // CACHE 3 - it is no use to cache becp, it will change in each SCF iteration
    // [in] psi, tab_atomic_, npw, becp, ik
//     const char transa = 'C';
//     const char transb = 'N';
//     const int ldb = this->npwx_;
//     const int ldc = this->tot_nproj;
//     const std::complex<double> alpha = 1.0;
//     const std::complex<double> beta = 0.0;
//     if(this->becp == nullptr || this->size_becp < npm*ldc)
//     {
//         delete[] this->becp;
//         this->becp = new std::complex<double>[npm*ldc];
//         this->size_becp = npm*ldc;
//     }
//     setmem_complex_op()(ctx, this->becp, 0.0, this->size_becp);
//     gemm_op()(
//         this->ctx,
//         transa,                 // const char transa
//         transb,                 // const char transb
//         ldc,                    // const int m
//         npm,                    // const int n
//         this->npw_,             // const int k
//         &alpha,                 // const std::complex<double> alpha
//         this->tab_atomic_,      // const std::complex<double>* a
//         this->npw_,             // const int lda
//         ppsi,                   // const std::complex<double>* b
//         ldb,                    // const int ldb
//         &beta,                  // const std::complex<double> beta
//         becp,                   // std::complex<double>* c
//         ldc);                   // const int ldc
// #ifdef __MPI
//     Parallel_Reduce::reduce_pool(becp, size_becp);
// #endif

    // notes on refactor for DCU calculation
    // the npm here is nbands(occ) * npol, for calling cal_becp, the npol should be divided.
    // std::cout << "npm: " << npm << std::endl;
    // std::cout << "at " << __FILE__ << ": " << __LINE__ << " output tot_nproj: " << this->tot_nproj << std::endl;
    // std::cout << "at " << __FILE__ << ": " << __LINE__ << " output npm: " << npm << std::endl;
    // std::cout << "at " << __FILE__ << ": " << __LINE__ << " ik_: " << ik_ << std::endl;
    int npol = this->ucell->get_npol();
    if(this->becp == nullptr || this->size_becp < npm*this->tot_nproj)
    {
        this->size_becp = npm*this->tot_nproj;
        resmem_complex_op()(this->becp, this->size_becp);
        if(this->device == base_device::GpuDevice )
        {
            resmem_complex_h_op()(this->h_becp, this->size_becp);
        }
        else
        {
            this->h_becp = this->becp;
        }
    }
    this->fs_tools->cal_becp(ik_, npm/npol, this->becp, ppsi); // in cal_becp, npm should be the one not multiplied by npol
    if(this->device == base_device::GpuDevice)
    {
        syncmem_complex_d2h_op()(h_becp, this->becp, this->size_becp);
    }
    ModuleBase::timer::tick("OnsiteProj", "overlap");
}

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::read_abacus_orb(std::ifstream& ifs,
                           std::string& elem,
                           double& ecut,
                           int& nr,
                           double& dr,
                           std::vector<int>& nzeta,
                           std::vector<std::vector<double>>& radials,
                           const int rank)
{
    nr = 0; // number of grid points
    dr = 0; // grid spacing
    int lmax = 0, nchi = 0; // number of radial functions
    std::vector<std::vector<int>> radial_map_; // build a map from [l][izeta] to 1-d array index
    std::string tmp;
    // first read the header
    if (rank == 0)
    {
        if (!ifs.is_open())
        {
            ModuleBase::WARNING_QUIT("AtomicRadials::read_abacus_orb", "Couldn't open orbital file.");
        }
        while (ifs >> tmp)
        {
            if (tmp == "Element")
            {
                ifs >> elem;
            }
            else if (tmp == "Cutoff(Ry)")
            {
                ifs >> ecut;
            }
            else if (tmp == "Lmax")
            {
                ifs >> lmax;
                nzeta.resize(lmax + 1);
                for (int l = 0; l <= lmax; ++l)
                {
                    ifs >> tmp >> tmp >> tmp >> nzeta[l];
                }
            }
            else if (tmp == "Mesh")
            {
                ifs >> nr;
                continue;
            }
            else if (tmp == "dr")
            {
                ifs >> dr;
                break;
            }
        }
        radial_map_.resize(lmax + 1);
        for (int l = 0; l <= lmax; ++l)
        {
            radial_map_[l].resize(nzeta[l]);
        }
        int ichi = 0;
        for (int l = 0; l <= lmax; ++l)
        {
            for (int iz = 0; iz < nzeta[l]; ++iz)
            {
                radial_map_[l][iz] = ichi++; // return the value of ichi, then increment
            }
        }
        nchi = ichi; // total number of radial functions
        radials.resize(nchi);
        std::for_each(radials.begin(), radials.end(), [nr](std::vector<double>& v) { v.resize(nr); });
    }

    // broadcast the header information
#ifdef __MPI
    Parallel_Common::bcast_string(elem);
    Parallel_Common::bcast_double(ecut);
    Parallel_Common::bcast_int(lmax);
    Parallel_Common::bcast_int(nchi);
    Parallel_Common::bcast_int(nr);
    Parallel_Common::bcast_double(dr);
#endif

    // then adjust the size of the vectors
    if (rank != 0)
    {
        nzeta.resize(lmax + 1);
        radials.resize(nchi);
        std::for_each(radials.begin(), radials.end(), [nr](std::vector<double>& v) { v.resize(nr); });
    }
    // broadcast the number of zeta functions for each angular momentum
#ifdef __MPI
    Parallel_Common::bcast_int(nzeta.data(), lmax + 1);
#endif

    // read the radial functions by rank0
    int ichi = 0;
    for (int i = 0; i != nchi; ++i)
    {
        if (rank == 0)
        {
            int l, izeta;
            ifs >> tmp >> tmp >> tmp;
            ifs >> tmp >> l >> izeta;
            ichi = radial_map_[l][izeta];
            for (int ir = 0; ir != nr; ++ir)
            {
                ifs >> radials[ichi][ir];
            }
        }
    // broadcast the radial functions
#ifdef __MPI
        Parallel_Common::bcast_int(ichi); // let other ranks know where to store the radial function
        Parallel_Common::bcast_double(radials[ichi].data(), nr);
#endif
    }
} // end of read_abacus_orb

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::cal_occupations(const psi::Psi<std::complex<T>, Device>* psi_in, const ModuleBase::matrix& wg_in)
{
    ModuleBase::timer::tick("OnsiteProj", "cal_occupation");
    this->tabulate_atomic(0);
    std::vector<std::complex<double>> occs(this->tot_nproj * 4, 0.0);

    // loop over k-points to calculate Mi of \sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
    const int nbands = psi_in->get_nbands();
    for(int ik = 0; ik < psi_in->get_nk(); ik++)
    {
        psi_in->fix_k(ik);
        if(ik != 0)
        {
            this->tabulate_atomic(ik);
        }
        // std::cout << __FILE__ << ":" << __LINE__ << " nbands = " << nbands << std::endl;
        this->overlap_proj_psi(
                        nbands * psi_in->npol,
                        psi_in->get_pointer());
        const std::complex<double>* becp_p = this->get_h_becp();
        // becp(nbands*npol , nkb)
        // mag = wg * \sum_{nh}becp * becp
        int nkb = this->tot_nproj;
        //nkb = 18;
        //std::cout << "at " << __FILE__ << ": " << __LINE__ << " output nbands: " << nbands << std::endl;
        //std::cout << "at " << __FILE__ << ": " << __LINE__ << " output nkb: " << nkb << std::endl;
        for(int ib = 0;ib<nbands;ib++)
        {
            const double weight = wg_in(ik, ib);
            int begin_ih = 0;
            for(int iat = 0; iat < this->iat_nh.size(); iat++)
            {
                const int nh = this->get_nh(iat);
                for(int ih = 0; ih < nh; ih++)
                {
                    const int occ_index = (begin_ih + ih) * 4;
                    const int index = ib*2*nkb + begin_ih + ih;
                    occs[occ_index] += weight * conj(becp_p[index]) * becp_p[index];
                    occs[occ_index + 1] += weight * conj(becp_p[index]) * becp_p[index + nkb];
                    occs[occ_index + 2] += weight * conj(becp_p[index + nkb]) * becp_p[index];
                    occs[occ_index + 3] += weight * conj(becp_p[index + nkb]) * becp_p[index + nkb];
                }
                begin_ih += nh;
            }
        }
    }
    // reduce mag from all k-pools
    Parallel_Reduce::reduce_double_allpool(GlobalV::KPAR, GlobalV::NPROC_IN_POOL, (double*)(&(occs[0])), occs.size()*2);
    // occ has been reduced and calculate mag
    // parameters for orbital charge output
    FmtCore fmt_of_chg("%15.4f");
    FmtCore fmt_of_label("%-15s");
    GlobalV::ofs_running << std::endl;
    GlobalV::ofs_running << "-------------------------------------------------------------------------------------------" << std::endl;
    GlobalV::ofs_running << "Orbital Charge Analysis      Charge         Mag(x)         Mag(y)         Mag(z)" << std::endl;
    GlobalV::ofs_running << "-------------------------------------------------------------------------------------------" << std::endl;
    // parameters for orbital charge output
    // parameters for mag output
    std::vector<double> mag_x(this->ucell->nat, 0.0);
    std::vector<double> mag_y(this->ucell->nat, 0.0);
    std::vector<double> mag_z(this->ucell->nat,0.0);
    auto atomLabels = this->ucell->get_atomLabels();
    const std::vector<std::string> title = {"Total Magnetism (uB)", "", "", ""};
    const std::vector<std::string> fmts = {"%-26s", "%20.10f", "%20.10f", "%20.10f"};
    const std::vector<std::string> orb_names = {"s", "p", "d", "f", "g"};
    FmtTable table(title, this->ucell->nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
    // parameters for mag output
    int occ_index = 0;
    for(int iat=0;iat<this->ucell->nat;iat++)
    {
        const int it = this->ucell->iat2it[iat];
        std::string atom_label = atomLabels[it];
        int ia = this->ucell->iat2ia[iat];
        GlobalV::ofs_running << FmtCore::format("%-20s", atom_label+std::to_string(ia+1)) << std::endl;
        std::vector<double> sum(4, 0.0);
        int current_l = 1;
        std::vector<double> charge_mag(4, 0.0);
        for(int ih=0;ih<this->iat_nh[iat];ih++)
        {
            charge_mag[3] += (occs[occ_index] - occs[occ_index + 3]).real();
            charge_mag[1] += (occs[occ_index + 1] + occs[occ_index + 2]).real();
            charge_mag[2] += (occs[occ_index + 1] - occs[occ_index + 2]).imag();
            charge_mag[0] += (occs[occ_index] + occs[occ_index + 3]).real();
            if(ih == current_l * current_l - 1)
            {
                sum[0] += charge_mag[0];
                sum[1] += charge_mag[1];
                sum[2] += charge_mag[2];
                sum[3] += charge_mag[3];
                GlobalV::ofs_running << FmtCore::format("%20s", orb_names[current_l-1])
                    << fmt_of_chg.format(charge_mag[0]) << fmt_of_chg.format(charge_mag[1])
                    << fmt_of_chg.format(charge_mag[2]) << fmt_of_chg.format(charge_mag[3]) << std::endl;
                current_l++;
                charge_mag.assign(4, 0.0);
            }
            occ_index += 4;
        }
        mag_x[iat] = sum[1];
        mag_y[iat] = sum[2];
        mag_z[iat] = sum[3];
        GlobalV::ofs_running << FmtCore::format("%20s", std::string("Sum")) << ""
                    << fmt_of_chg.format(sum[0]) << fmt_of_chg.format(sum[1])
                    << fmt_of_chg.format(sum[2]) << fmt_of_chg.format(sum[3]) << std::endl;
    }
    GlobalV::ofs_running << "-------------------------------------------------------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::endl;
    table << atomLabels << mag_x << mag_y << mag_z;
    GlobalV::ofs_running << table.str() << std::endl;
    
    // print charge
    ModuleBase::timer::tick("OnsiteProj", "cal_occupation");
}

template class projectors::OnsiteProjector<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class projectors::OnsiteProjector<double, base_device::DEVICE_GPU>;
#endif