#include "vnl_pw.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/math_ylmreal.h"
#include "source_base/timer.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_device/device.h"
#include "source_base/parallel_reduce.h"
#include "source_pw/module_pwdft/kernels/vnl_op.h"

#include <cmath>

/**
 * @file vnl_pw_deeq.cpp
 * @brief Update the effective D coefficients (deeq / deeq_nc) for USPP
 *        self-consistent iterations.
 *
 * This file contains:
 * - cal_effective_D(): main entry, dispatch between KB and USPP paths
 * - newq(): compute integral of V_eff * Q via FFT + gemm
 * - newd_so(): rotate deeq into deeq_nc for SOC case
 * - newd_nc(): rotate deeq into deeq_nc for non-collinear case
 */

/**
 * @brief Recalculate effective coefficient matrix for non-local pseudopotential.
 *
 * For norm-conserving (KB) pseudopotentials, copy dvan into deeq directly.
 * For Vanderbilt USPP, call newq() to add the augmentation contribution,
 * then rotate into deeq_nc if SOC/non-collinear is enabled.
 *
 * @param veff     effective potential on real-space grid
 * @param rho_basis plane-wave basis for charge density FFT
 * @param cell     unit cell (read-only, but passed by reference for compatibility)
 */
void pseudopot_cell_vnl::cal_effective_D(const ModuleBase::matrix& veff,
                                         const ModulePW::PW_Basis* rho_basis,
                                         UnitCell& cell)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "cal_effective_D");

    /*
    recalculate effective coefficient matrix for non-local pseudo-potential
    1. assign to each atom from element;
    2. extend to each spin when nspin larger than 1
    3. rotate to effective matrix when spin-orbital coupling is used
    */

    if (!PARAM.globalv.use_uspp)
    {
        for (int iat = 0; iat < cell.nat; iat++)
        {
            const int it = cell.iat2it[iat];
            const int nht = cell.atoms[it].ncpp.nh;
            // nht: number of beta functions per atom type
            for (int is = 0; is < PARAM.inp.nspin; is++)
            {
                for (int ih = 0; ih < nht; ih++)
                {
                    for (int jh = ih; jh < nht; jh++)
                    {
                        if (PARAM.inp.lspinorb)
                        {
                            this->deeq_nc(is, iat, ih, jh) = this->dvan_so(is, it, ih, jh);
                            this->deeq_nc(is, iat, jh, ih) = this->dvan_so(is, it, jh, ih);
                        }
                        else if (PARAM.inp.nspin == 4)
                        {
                            if (is == 0)
                            {
                                this->deeq_nc(is, iat, ih, jh) = this->dvan(it, ih, jh);
                                this->deeq_nc(is, iat, jh, ih) = this->dvan(it, ih, jh);
                            }
                            else if (is == 1)
                            {
                                this->deeq_nc(is, iat, ih, jh) = std::complex<double>(0.0, 0.0);
                                this->deeq_nc(is, iat, jh, ih) = std::complex<double>(0.0, 0.0);
                            }
                            else if (is == 2)
                            {
                                this->deeq_nc(is, iat, ih, jh) = std::complex<double>(0.0, 0.0);
                                this->deeq_nc(is, iat, jh, ih) = std::complex<double>(0.0, 0.0);
                            }
                            else if (is == 3)
                            {
                                this->deeq_nc(is, iat, ih, jh) = this->dvan(it, ih, jh);
                                this->deeq_nc(is, iat, jh, ih) = this->dvan(it, ih, jh);
                            }
                        }
                        else
                        {
                            this->deeq(is, iat, ih, jh) = this->dvan(it, ih, jh);
                            this->deeq(is, iat, jh, ih) = this->dvan(it, ih, jh);
                            // in most of pseudopotential files, number of projections of one orbital is only one,
                            // which lead to diagonal matrix of dion
                            // when number larger than 1, non-diagonal dion should be calculated.
                            if (ih != jh && std::fabs(this->deeq(is, iat, ih, jh)) > 0.0)
                            {
                                this->multi_proj = true;
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        newq(veff, rho_basis, cell);

        for (int iat = 0; iat < cell.nat; iat++)
        {
            int it = cell.iat2it[iat];
            if (PARAM.inp.noncolin)
            {
                if (cell.atoms[it].ncpp.has_so)
                {
                    this->newd_so(iat, cell);
                }
                else
                {
                    this->newd_nc(iat, cell);
                }
            }
            else
            {
                for (int is = 0; is < PARAM.inp.nspin; is++)
                {
                    for (int ih = 0; ih < cell.atoms[it].ncpp.nh; ih++)
                    {
                        for (int jh = ih; jh < cell.atoms[it].ncpp.nh; jh++)
                        {
                            deeq(is, iat, ih, jh) += this->dvan(it, ih, jh);
                            deeq(is, iat, jh, ih) = deeq(is, iat, ih, jh);
                        }
                    }
                }
            }
        }
    }
    if (this->use_gpu_)
    {
        if (PARAM.globalv.has_float_data)
        {
            castmem_d2s_h2d_op()(this->s_deeq,
                                 this->deeq.ptr,
                                 PARAM.inp.nspin * cell.nat * this->nhm * this->nhm);
            castmem_z2c_h2d_op()(this->c_deeq_nc,
                                 this->deeq_nc.ptr,
                                 PARAM.inp.nspin * cell.nat * this->nhm * this->nhm);
        }
        if (PARAM.globalv.has_double_data)
        {
            syncmem_z2z_h2d_op()(this->z_deeq_nc,
                                 this->deeq_nc.ptr,
                                 PARAM.inp.nspin * cell.nat * this->nhm * this->nhm);
        }
        syncmem_d2d_h2d_op()(this->d_deeq,
                             this->deeq.ptr,
                             PARAM.inp.nspin * cell.nat * this->nhm * this->nhm);
    }
    else
    {
        if (PARAM.globalv.has_float_data)
        {
            castmem_d2s_h2h_op()(this->s_deeq,
                                 this->deeq.ptr,
                                 PARAM.inp.nspin * cell.nat * this->nhm * this->nhm);
            castmem_z2c_h2h_op()(this->c_deeq_nc,
                                 this->deeq_nc.ptr,
                                 PARAM.inp.nspin * cell.nat * this->nhm * this->nhm);
        }
        // There's no need to synchronize double precision pointers while in a CPU environment.
    }
}

/**
 * @brief Compute the augmentation charge contribution to deeq.
 *
 * For each USPP atom type, evaluate Q(G) on the FFT grid via radial_fft_q(),
 * then compute the integral integral[ V_eff(r) * Q_{ij}(r) dr ] using gemm.
 * The result is accumulated into deeq.
 *
 * @param veff     effective potential on real-space grid
 * @param rho_basis plane-wave basis for charge density FFT
 * @param cell     unit cell
 */
void pseudopot_cell_vnl::newq(const ModuleBase::matrix& veff, const ModulePW::PW_Basis* rho_basis, UnitCell& cell)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "newq");

    const std::complex<double> ci_tpi = ModuleBase::IMAG_UNIT * ModuleBase::TWO_PI;
    double fact = 1.0;
    if (rho_basis->gamma_only)
    {
        fact = 2.0;
    }

    deeq.zero_out();

    const int npw = rho_basis->npw;
    ModuleBase::matrix ylmk0(lmaxq * lmaxq, npw);
    ModuleBase::YlmReal::Ylm_Real(lmaxq * lmaxq, npw, rho_basis->gcar, ylmk0);

    std::vector<double> qnorm(npw);
    for (int ig = 0; ig < npw; ig++)
    {
        qnorm[ig] = rho_basis->gcar[ig].norm() * cell.tpiba;
    }

    // fourier transform of the total effective potential
    ModuleBase::ComplexMatrix vaux(PARAM.inp.nspin, npw);
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        rho_basis->real2recip(&veff.c[is * veff.nc], &vaux(is, 0));
    }

    for (int it = 0; it < cell.ntype; it++)
    {
        Atom_pseudo* upf = &cell.atoms[it].ncpp;
        if (upf->tvanp)
        {
            // nij = max number of (ih,jh) pairs per atom type
            int nij = upf->nh * (upf->nh + 1) / 2;
            ModuleBase::ComplexMatrix qg(nij, npw);

            // Compute and store Q(G) for this atomic species
            // (without structure factor)
            int ijh = 0;
            for (int ih = 0; ih < upf->nh; ih++)
            {
                for (int jh = ih; jh < upf->nh; jh++)
                {
                    radial_fft_q(npw, ih, jh, it, qnorm.data(), ylmk0, &qg(ijh, 0));
                    ijh++;
                }
            }

            // Compute and store V(G) times the structure factor e^(-iG*tau)
            const int natom = cell.atoms[it].na;
            ModuleBase::ComplexMatrix aux(natom, npw);
            ModuleBase::matrix deeaux(natom, nij);
            for (int is = 0; is < PARAM.inp.nspin; is++)
            {
                for (int ia = 0; ia < natom; ia++)
                {
                    const ModuleBase::Vector3<double> tau = cell.atoms[it].tau[ia];
                    for (int ig = 0; ig < npw; ig++)
                    {
                        const ModuleBase::Vector3<double> g = rho_basis->gcar[ig];
                        const std::complex<double> phase = ci_tpi * (g * tau);
                        aux(ia, ig) = vaux(is, ig) * exp(phase);
                    }
                }
                // here we compute the integral Q*V for all atoms of this kind
                const char transa = 'C';
                const char transb = 'N';
                const double zero = 0.0;
                const int complex_npw = 2 * npw;
                double* qg_ptr = reinterpret_cast<double*>(qg.c);
                double* aux_ptr = reinterpret_cast<double*>(aux.c);

                BlasConnector::gemm(transb,
                       transa,
                       natom,
                       nij,
                       complex_npw,
                       fact,
                       aux_ptr,
                       complex_npw,
                       qg_ptr,
                       complex_npw,
                       zero,
                       deeaux.c,
                       nij);
                // I'm not sure if this is correct for gamma_only
                if (rho_basis->gamma_only && rho_basis->ig_gge0 >= 0)
                {
                    const double neg = -1.0;
                    dger_(&nij, &natom, &neg, qg_ptr, &complex_npw, aux_ptr, &complex_npw, deeaux.c, &nij);
                }

                for (int ia = 0; ia < natom; ia++)
                {
                    int ijh = 0;
                    const int iat = cell.itia2iat(it, ia);
                    for (int ih = 0; ih < upf->nh; ih++)
                    {
                        for (int jh = ih; jh < upf->nh; jh++)
                        {
                            deeq(is, iat, ih, jh) = cell.omega * deeaux(ia, ijh);
                            if (jh > ih)
                            {
                                deeq(is, iat, jh, ih) = deeq(is, iat, ih, jh);
                            }
                            ijh++;
                        }
                    }
                }
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_pool(deeq.ptr,deeq.getSize());
#endif
}

/**
 * @brief Rotate deeq into deeq_nc for spin-orbit coupling case.
 *
 * Uses the SOC rotation coefficients soc.fcoef to transform
 * the spin-diagonal deeq into the full 2x2 spinor deeq_nc.
 *
 * @param iat  atom index
 * @param cell unit cell
 */
void pseudopot_cell_vnl::newd_so(const int& iat, UnitCell& cell)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "newd_so");

    const int it = cell.iat2it[iat];
    Atom_pseudo* upf = &cell.atoms[it].ncpp;
    int ijs = 0;
    for (int is1 = 0; is1 < 2; is1++)
    {
        for (int is2 = 0; is2 < 2; is2++)
        {
            for (int ih = 0; ih < upf->nh; ih++)
            {
                for (int jh = 0; jh < upf->nh; jh++)
                {
                    deeq_nc(ijs, iat, ih, jh) = dvan_so(ijs, it, ih, jh);

                    for (int kh = 0; kh < upf->nh; kh++)
                    {
                        for (int lh = 0; lh < upf->nh; lh++)
                        {
                            if (PARAM.globalv.domag)
                            {
                                deeq_nc(ijs, iat, ih, jh)
                                    += deeq(0, iat, kh, lh)
                                           * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 0, is2, lh, jh)
                                              + soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 1, is2, lh, jh))
                                       + deeq(1, iat, kh, lh)
                                             * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 1, is2, lh, jh)
                                                + soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 0, is2, lh, jh))
                                       + ModuleBase::NEG_IMAG_UNIT * deeq(2, iat, kh, lh)
                                             * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 1, is2, lh, jh)
                                                - soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 0, is2, lh, jh))
                                       + deeq(3, iat, kh, lh)
                                             * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 0, is2, lh, jh)
                                                - soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 1, is2, lh, jh));
                            }
                            else
                            {
                                deeq_nc(ijs, iat, ih, jh)
                                    += deeq(0, iat, kh, lh)
                                       * (soc.fcoef(it, is1, 0, ih, kh) * soc.fcoef(it, 0, is2, lh, jh)
                                          + soc.fcoef(it, is1, 1, ih, kh) * soc.fcoef(it, 1, is2, lh, jh));
                            }
                        }
                    }
                }
            }
            ijs++;
        }
    }
}

/**
 * @brief Rotate deeq into deeq_nc for non-collinear magnetic case.
 *
 * Combines the scalar deeq components into the Pauli-matrix representation
 * of deeq_nc for non-collinear magnetism without SOC.
 *
 * @param iat  atom index
 * @param cell unit cell
 */
void pseudopot_cell_vnl::newd_nc(const int& iat, UnitCell& cell)
{
    ModuleBase::TITLE("pseudopot_cell_vnl", "newd_nc");

    const int it = cell.iat2it[iat];
    Atom_pseudo* upf = &cell.atoms[it].ncpp;

    for (int ih = 0; ih < upf->nh; ih++)
    {
        for (int jh = 0; jh < upf->nh; jh++)
        {
            if (PARAM.inp.lspinorb)
            {
                deeq_nc(0, iat, ih, jh) = dvan_so(0, it, ih, jh) + deeq(0, iat, ih, jh) + deeq(3, iat, ih, jh);
                deeq_nc(3, iat, ih, jh) = dvan_so(3, it, ih, jh) + deeq(0, iat, ih, jh) - deeq(3, iat, ih, jh);
            }
            else
            {
                deeq_nc(0, iat, ih, jh) = dvan(it, ih, jh) + deeq(0, iat, ih, jh) + deeq(3, iat, ih, jh);
                deeq_nc(3, iat, ih, jh) = dvan(it, ih, jh) + deeq(0, iat, ih, jh) - deeq(3, iat, ih, jh);
            }
            deeq_nc(1, iat, ih, jh) = deeq(1, iat, ih, jh) + ModuleBase::NEG_IMAG_UNIT * deeq(2, iat, ih, jh);
            deeq_nc(2, iat, ih, jh) = deeq(1, iat, ih, jh) + ModuleBase::IMAG_UNIT * deeq(2, iat, ih, jh);
        }
    }
}
