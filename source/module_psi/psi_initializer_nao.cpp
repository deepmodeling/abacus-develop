#include "psi_initializer_nao.h"

#include <fstream>
// numerical algorithm support
#include "module_base/math_integral.h" // for numerical integration
// numerical algorithm support
#include "module_base/math_polyint.h" // for polynomial interpolation
#include "module_base/math_ylmreal.h" // for real spherical harmonics
// basic functions support
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
// three global variables definition
#include "module_base/global_variable.h"
// parallel communication
#ifdef __MPI
#include "module_base/parallel_common.h"
#include "module_base/parallel_reduce.h"
#endif
#include "module_io/orb_io.h"
#include "module_parameter/parameter.h"
// GlobalV::NQX and GlobalV::DQ are here
#include "module_parameter/parameter.h"

#include <algorithm>
#include <numeric>

/*
I don't know why some variables are distributed while others not... for example the orbital_files...
We need not only read and import, but also distribute here
*/

// free function, not needed to be a member of psi_initializer_nao
void normalize(const std::vector<double>& r, std::vector<double>& flz)
{
    std::vector<double> flz2r2(r.size());
    std::transform(r.begin(), r.end(), flz.begin(), flz2r2.begin(), [](double r, double flz) {
        return flz * flz * r * r;
    });
    double dr = r[1] - r[0];
    double norm = ModuleBase::Integral::simpson(r.size(), flz2r2.data(), dr);
    norm = sqrt(norm);
    std::transform(flz.begin(), flz.end(), flz.begin(), [norm](double flz) { return flz / norm; });
}

template <typename T>
void psi_initializer_nao<T>::read_external_orbs(const std::string* orbital_files, const int& rank)
{
    ModuleBase::timer::tick("psi_initializer_nao", "read_external_orbs");

    this->orbital_files_.resize(this->p_ucell_->ntype);
    this->nr_.resize(this->p_ucell_->ntype);
    this->rgrid_.resize(this->p_ucell_->ntype);
    this->chi_.resize(this->p_ucell_->ntype);

#ifdef __MPI
    if (rank == 0)
    {
#endif
        std::copy(orbital_files, orbital_files + this->p_ucell_->ntype, this->orbital_files_.begin());
#ifdef __MPI
    }
    Parallel_Common::bcast_string(this->orbital_files_.data(), this->p_ucell_->ntype);
#endif
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        std::ifstream ifs_it;
        bool is_open = false;
        if (rank == 0)
        {
            ifs_it.open(PARAM.inp.orbital_dir + this->orbital_files_[it]);
            is_open = ifs_it.is_open();
        }
#ifdef __MPI
        Parallel_Common::bcast_bool(is_open);
#endif
        if (!is_open)
        {
            GlobalV::ofs_warning << "psi_initializer_nao<T>::read_orbital_files: cannot open orbital file: "
                                 << this->orbital_files_[it] << std::endl;
            ModuleBase::WARNING_QUIT("psi_initializer_nao<T>::read_orbital_files", "cannot open orbital file.");
        }
        else
        {
            GlobalV::ofs_running << "psi_initializer_nao<T>::read_orbital_files: reading orbital file: "
                                 << this->orbital_files_[it] << std::endl;
        }
        std::string elem; // garbage value, will discard
        double ecut;      // garbage value, will discard
        int nr;
        double dr;
        std::vector<int> nzeta;
        std::vector<std::vector<double>> radials;
        ModuleIO::read_abacus_orb(ifs_it, elem, ecut, nr, dr, nzeta, radials, rank);

        if (rank == 0)
        {
            ifs_it.close();
        }

        const int nchi = std::accumulate(nzeta.begin(), nzeta.end(), 0);
        // nr_
        this->nr_[it].resize(nchi);
        std::for_each(this->nr_[it].begin(), this->nr_[it].end(), [nr](int& numr) { numr = nr; });
        // rgrid_
        this->rgrid_[it].resize(nchi);
        std::for_each(this->rgrid_[it].begin(), this->rgrid_[it].end(), [nr, dr](std::vector<double>& rgrid) {
            rgrid.resize(nr);
            std::iota(rgrid.begin(), rgrid.end(), 0);
            std::for_each(rgrid.begin(), rgrid.end(), [dr](double& r) { r = r * dr; });
        });
        // chi_
        this->chi_[it].resize(nchi);
        std::for_each(this->chi_[it].begin(), this->chi_[it].end(), [nr](std::vector<double>& chi) { chi.resize(nr); });
        for (int ichi = 0; ichi < nchi; ichi++)
        {
            std::copy(radials[ichi].begin(), radials[ichi].end(), this->chi_[it][ichi].begin());
        }
    }
    ModuleBase::timer::tick("psi_initializer_nao", "read_external_orbs");
}

template <typename T>
void psi_initializer_nao<T>::allocate_ao_table()
{
    // find correct dimension for ovlp_flzjlq
    int ntype = this->p_ucell_->ntype;
    int lmaxmax = 0;   // lmaxmax
    int nzeta_max = 0; // dim3 should be the maximum number of zeta for each atomtype
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        int nzeta = 0;
        int lmax = this->p_ucell_->atoms[it].nwl;
        lmaxmax = (lmaxmax > lmax) ? lmaxmax : lmax;
        for (int l = 0; l < lmax + 1; l++)
        {
            nzeta += this->p_ucell_->atoms[it].l_nchi[l];
        }
        nzeta_max = (nzeta > nzeta_max) ? nzeta : nzeta_max;
    }
    if (nzeta_max == 0)
    {
        ModuleBase::WARNING_QUIT("psi_initializer_nao<T>::psi_initializer_nao",
                                 "there is not ANY numerical atomic orbital read in present system, quit.");
    }
    // allocate a map (it, l, izeta) -> i, should allocate memory of ntype * lmax * nzeta_max
    this->projmap_.create(ntype, lmaxmax + 1, nzeta_max);
}

template <typename T>
void psi_initializer_nao<T>::initialize(const Structure_Factor* sf,
                                        const ModulePW::PW_Basis_K* pw_wfc,
                                        const UnitCell* p_ucell,
                                        const K_Vectors* p_kv_in,
                                        const int& random_seed,
                                        const pseudopot_cell_vnl* p_pspot_nl,
                                        const int& rank)
{
    ModuleBase::timer::tick("psi_initializer_nao", "initialize_mpi");

    // import
    psi_initializer<T>::initialize(sf, pw_wfc, p_ucell, p_kv_in, random_seed, p_pspot_nl, rank);

    // allocate
    this->allocate_ao_table();
    this->read_external_orbs(this->p_ucell_->orbital_fn.data(), rank);

    // then for generate random number to fill in the wavefunction
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());

    int nbands_local = 0;
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        for (int l = 0; l < this->p_ucell_->atoms[it].nwl + 1; l++)
        {
            /* EVERY ZETA FOR (2l+1) ORBS */
            const int nchi = this->p_ucell_->atoms[it].l_nchi[l];
            const int degen_l = (l == 0) ? 1 : 2 * l + 1;
            nbands_local += nchi * degen_l * PARAM.globalv.npol * this->p_ucell_->atoms[it].na;
            /*
                non-rotate basis, nbands_local*=2 for PARAM.globalv.npol = 2 is enough
            */
            // nbands_local += this->p_ucell_->atoms[it].l_nchi[l]*(2*l+1) * PARAM.globalv.npol;
            /*
                rotate basis, nbands_local*=4 for p, d, f,... orbitals, and nbands_local*=2 for s orbitals
                risky when NSPIN = 4, problematic psi value, needed to be checked
            */
        }
    }
    this->nbands_start_ = std::max(nbands_local, PARAM.inp.nbands);
    this->nbands_complem_ = this->nbands_start_ - nbands_local;

    ModuleBase::timer::tick("psi_initializer_nao", "initialize_mpi");
}

template <typename T>
void psi_initializer_nao<T>::tabulate()
{
    ModuleBase::timer::tick("psi_initializer_nao", "tabulate");

    // a uniformed qgrid
    std::vector<double> qgrid(PARAM.globalv.nqx);
    std::iota(qgrid.begin(), qgrid.end(), 0);
    std::for_each(qgrid.begin(), qgrid.end(), [this](double& q) { q = q * PARAM.globalv.dq; });

    // only when needed, allocate memory for cubspl_
    if (this->cubspl_.get())
    {
        this->cubspl_.reset();
    }
    this->cubspl_ = std::unique_ptr<ModuleBase::CubicSpline>(new ModuleBase::CubicSpline(qgrid.size(), qgrid.data()));

    // calculate the total number of radials and call reserve to allocate memory
    int nchi = 0;
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        for (int l = 0; l < this->p_ucell_->atoms[it].nwl + 1; l++)
        {
            nchi += this->p_ucell_->atoms[it].l_nchi[l];
        }
    }
    this->cubspl_->reserve(nchi);
    ModuleBase::SphericalBesselTransformer sbt_(true); // bool: enable cache

    // tabulate the spherical bessel transform of numerical orbital function
    std::vector<double> Jlfq(PARAM.globalv.nqx, 0.0);
    int i = 0;
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        int ic = 0;
        for (int l = 0; l < this->p_ucell_->atoms[it].nwl + 1; l++)
        {
            for (int izeta = 0; izeta < this->p_ucell_->atoms[it].l_nchi[l]; izeta++)
            {
                sbt_.direct(l,
                            this->nr_[it][ic],
                            this->rgrid_[it][ic].data(),
                            this->chi_[it][ic].data(),
                            PARAM.globalv.nqx,
                            qgrid.data(),
                            Jlfq.data());
                this->cubspl_->add(Jlfq.data());
                this->projmap_(it, l, izeta) = i++; // index it
                ++ic;
            }
        }
    }
    ModuleBase::timer::tick("psi_initializer_nao", "tabulate");
}

template <typename T>
void psi_initializer_nao<T>::init_psig(T* psig, const int& ik)
{
    ModuleBase::timer::tick("psi_initializer_nao", "initialize");
    assert(ik >= 0);
    const int npw = this->pw_wfc_->npwk[ik];
    const int npwk_max = this->pw_wfc_->npwk_max;
    const int total_lm = (this->p_ucell_->lmax + 1) * (this->p_ucell_->lmax + 1);
    ModuleBase::matrix ylm(total_lm, npw);
    ModuleBase::GlobalFunc::ZEROS(psig, PARAM.globalv.npol * this->nbands_start_ * npwk_max);

    std::vector<std::complex<double>> aux(npw);
    std::vector<double> qnorm(npw);
    std::vector<ModuleBase::Vector3<double>> q(npw);

#pragma omp parallel for schedule(static, 4096 / sizeof(double))
    for (int ig = 0; ig < npw; ig++)
    {
        q[ig] = this->pw_wfc_->getgpluskcar(ik, ig);
        qnorm[ig] = q[ig].norm() * this->p_ucell_->tpiba;
    }

    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, q.data(), ylm);
    // int index = 0;
    std::vector<double> Jlfq(npw, 0.0);
    int ibasis = 0;
    for (int it = 0; it < this->p_ucell_->ntype; it++)
    {
        /* HERE LOOP OVER ALL TYPES */
        for (int ia = 0; ia < this->p_ucell_->atoms[it].na; ia++)
        {
            /* HERE LOOP OVER ALL ATOMS */
            std::complex<double>* sk = this->sf_->get_sk(ik, it, ia, this->pw_wfc_);
            int ic = 0; // ic is a flatten index of chi, therefore it is defined here.
            for (int L = 0; L < this->p_ucell_->atoms[it].nwl + 1; L++)
            {
                std::complex<double> lphase = pow(ModuleBase::NEG_IMAG_UNIT, L); // mohan 2010-04-19
                for (int N = 0; N < this->p_ucell_->atoms[it].l_nchi[L]; N++)
                {
                    /* HERE LOOP OVER ALL NAOS */
                    /*
                        for already using flattened 1d index of chi, which folds l and n, the spherical bessel
                        transformation of numerical orbital function, is indiced by it and ic, is needed to
                        interpolate everytime when ic updates, therefore everytime when present orbital is done
                    */

                    // use cublic spline instead of previous polynomial interpolation
                    this->cubspl_->eval(npw, qnorm.data(), Jlfq.data(), nullptr, nullptr, this->projmap_(it, L, N));

                    /* FOR EVERY NAO IN EACH ATOM */
                    if (PARAM.inp.nspin == 4)
                    {
                        /* FOR EACH SPIN CHANNEL */
                        for (int is_N = 0; is_N < 2; is_N++) // rotate base
                        // for(int is_N = 0; is_N < 1; is_N++)
                        {
                            if (L == 0 && is_N == 1)
                            {
                                continue;
                            }
                            else
                            {
                                const double j = fabs(double(L + is_N) - 0.5);
                                double alpha, gamma;
                                std::complex<double> fup, fdown;
                                if (fabs(j - L + 0.5) < 1e-4)
                                {
                                    continue;
                                }
                                alpha = this->p_ucell_->atoms[it].angle1[ia];
                                gamma = -1 * this->p_ucell_->atoms[it].angle2[ia] + 0.5 * ModuleBase::PI;
                                for (int m = 0; m < 2 * L + 1; m++)
                                {
                                    const int lm = L * L + m;
#pragma omp parallel for
                                    for (int ig = 0; ig < npw; ig++)
                                    {
                                        aux[ig] = sk[ig] * ylm(lm, ig) * Jlfq[ig];
                                    }

#pragma omp parallel for
                                    for (int ig = 0; ig < npw; ig++)
                                    {
                                        fup = cos(0.5 * alpha) * aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5 * alpha) * aux[ig];
                                        // build the orthogonal wfc
                                        // first rotation with angle (alpha + ModuleBase::PI) around (OX)
                                        psig[ibasis * 2 * npwk_max + ig] = this->template cast_to_T<T>(
                                            (cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fup);
                                        psig[(ibasis * 2 + 1) * npwk_max + ig] = this->template cast_to_T<T>(
                                            (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown);
                                        // second rotation with angle gamma around(OZ)
                                        fup = cos(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                        fdown = ModuleBase::IMAG_UNIT * sin(0.5 * (alpha + ModuleBase::PI)) * aux[ig];
                                        psig[(ibasis + 2 * L + 1) * 2 * npwk_max + ig] = this->template cast_to_T<T>(
                                            (cos(0.5 * gamma) + ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fup);
                                        psig[((ibasis + 2 * L + 1) * 2 + 1) * npwk_max + ig] = this->template cast_to_T<T>(
                                            (cos(0.5 * gamma) - ModuleBase::IMAG_UNIT * sin(0.5 * gamma)) * fdown);
                                    }
                                    ibasis++;
                                }
                                ibasis += 2 * L + 1;
                            }
                        } // end for is_N
                    }     // end if PARAM.inp.noncolin
                    else
                    { // LSDA and nomagnet case
                        /* DOES NOT DISTINGUISH m QUANTUM NUMBER FOR CHI */
                        for (int m = 0; m < 2 * L + 1; m++)
                        {
                            const int lm = L * L + m;
#pragma omp parallel for
                            for (int ig = 0; ig < npw; ig++)
                            {
                                psig[ibasis * npwk_max + ig] = this->template cast_to_T<T>(lphase * sk[ig] * ylm(lm, ig) * Jlfq[ig]);
                            }
                            ++ibasis;
                        }
                    }
                    ++ic;
                } // end for N
            }     // end for L
            delete[] sk;
        } // end for ia
    }     // end for it
    /* complement the rest of bands if there are */
    if (this->nbands_complem() > 0)
    {
        this->random_t(psig, ibasis, this->nbands_start_, ik);
    }
    ModuleBase::timer::tick("psi_initializer_nao", "initialize");
}

template class psi_initializer_nao<std::complex<double>>;
template class psi_initializer_nao<std::complex<float>>;
// gamma point calculation
template class psi_initializer_nao<double>;
template class psi_initializer_nao<float>;
