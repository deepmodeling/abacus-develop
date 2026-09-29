#include <cassert>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <map>
#include <tuple>

#include "source_base/module_out/orb_io.h"
#include "source_pw/module_proj/onsite_proj.h"
#include "source_pw/module_proj/onsite_proj_print.h"
#include "source_pw/module_proj/radial_proj.h"
#include "source_base/projgen.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/tool_quit.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

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
void projectors::OnsiteProjector<T, Device>::init(const std::string& orbital_dir,
        const UnitCell* ucell_in,
        const psi::Psi<std::complex<T>, Device>& psi,
        const K_Vectors& kv,
        const ModulePW::PW_Basis_K& pw_basis, // level1: the plane wave basis, need ik
        Structure_Factor& sf,                 // level2: the structure factor calculator
        const double onsite_radius,
        const int nq,
        const double dq,
        const ModuleBase::matrix& wg,
        const ModuleBase::matrix& ekb)
{
    this->device = base_device::get_device_type(this->ctx);

    if(!this->initialed)
    {
        this->ucell = ucell_in;
        this->ntype = ucell_in->ntype;
        this->isk_ = kv.isk.data();

        this->pw_basis_ = &pw_basis;
        this->sf_ = &sf;

        std::vector<std::string> orb_files(ntype);
        std::vector<int> nproj(ntype);
        int sum_nproj = 0;
        for(int it=0; it<ntype; ++it)
        {
            orb_files[it] = ucell->orbital_fn[it];
            nproj[it] = ucell->atoms[it].nwl;
            sum_nproj += nproj[it];
        }
        this->lproj.resize(sum_nproj);
        int index = 0;
        for(int it=0; it<ntype; ++it)
        {
            for(int il=0; il<nproj[it]; ++il)
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

        ModuleBase::timer::start("OnsiteProj", "cubspl_tabulate");
        // STAGE 0 - making the interpolation table
        // CACHE 0 - if cache the irow2it, irow2iproj, irow2m, itiaiprojm2irow, <G+k|p> can be reused for
        //           SCF, RELAX and CELL-RELAX calculation
        // [in] rgrid, projs, lproj, it2ia, it2iproj, nq, dq
        RadialProjection::build_backward_map(it2iproj, lproj, irow2it_, irow2iproj_, irow2m_);
        RadialProjection::build_forward_map(it2ia, it2iproj, lproj, itiaiprojm2irow_);
        RadialProjection::build_sbt_tab(nproj, rgrid, projs, lproj, nq, dq, ucell_in->omega, psi.get_npol(), tab, nhtol);
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

        ModuleBase::timer::end("OnsiteProj", "cubspl_tabulate");

        this->initialed = true;
    }
}

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::init_proj(const std::string& orbital_dir,
        const std::vector<std::string>& orb_files,
        const std::vector<int>& nproj,  // for each type, the number of projectors
        const std::vector<int>& lproj,  // angular momentum of projectors within the type (l of zeta function)
        const std::vector<int>& iproj,  // index of projectors within the type (izeta)
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
        print::print_proj_status(it, nproj_it);
        if(nproj_it == 0)
        {
            continue;
        }
        std::ifstream ifs(orbital_dir + orb_files[it]);
        std::string elem = "";
        double ecut = -1.0;
        int nr_ = -1;
        double dr_ = -1.0;
        std::vector<int> nzeta; // number of radials for each l
        std::vector<std::vector<double>> radials; // radials arranged in serial
        ModuleIO::read_abacus_orb(ifs, elem, ecut, nr_, dr_, nzeta, radials);
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

// explicit method instantiation
template
void projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::init(
    const std::string&,
    const UnitCell*,
    const psi::Psi<std::complex<double>, base_device::DEVICE_CPU>&,
    const K_Vectors&,
    const ModulePW::PW_Basis_K&,
    Structure_Factor&,
    const double,
    const int,
    const double,
    const ModuleBase::matrix&,
    const ModuleBase::matrix&);

template
void projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::init_proj(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<int>&,
    const std::vector<int>&,
    const std::vector<int>&,
    const std::vector<double>&);

#if ((defined __CUDA) || (defined __ROCM))
template
void projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::init(
    const std::string&,
    const UnitCell*,
    const psi::Psi<std::complex<double>, base_device::DEVICE_GPU>&,
    const K_Vectors&,
    const ModulePW::PW_Basis_K&,
    Structure_Factor&,
    const double,
    const int,
    const double,
    const ModuleBase::matrix&,
    const ModuleBase::matrix&);

template
void projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::init_proj(
    const std::string&,
    const std::vector<std::string>&,
    const std::vector<int>&,
    const std::vector<int>&,
    const std::vector<int>&,
    const std::vector<double>&);
#endif
