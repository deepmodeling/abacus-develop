#include "module_hamilt_pw/hamilt_pwdft/kernels/force_op.h"
#include "module_base/libm/libm.h"
#include "module_base/tool_threading.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hamilt
{

template <typename FPTYPE>
struct cal_vkb1_nl_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nkb,
                    const int& npwx,
                    const int& vkb_nc,
                    const int& nbasis,
                    const int& ipol,
                    const std::complex<FPTYPE>& NEG_IMAG_UNIT,
                    const std::complex<FPTYPE>* vkb,
                    const FPTYPE* gcar,
                    std::complex<FPTYPE>* vkb1)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int i = 0; i < nkb; i++)
        {
            for (int ig = 0; ig < nbasis; ig++)
            {
                std::complex<FPTYPE>* pvkb1 = vkb1 + i * npwx;
                const std::complex<FPTYPE>* pvkb = vkb + i * vkb_nc;
                pvkb1[ig] = pvkb[ig] * NEG_IMAG_UNIT * gcar[ig * 3 + ipol];
            }
        }
    }
};

template <typename FPTYPE>
struct cal_force_nl_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const bool& nondiagonal,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& spin,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const FPTYPE* deeq,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force)
    {
#ifdef _OPENMP
#pragma omp parallel
        {
#endif
            int iat0 = 0;
            int sum0 = 0;
            for (int it = 0; it < ntype; it++)
            {
                const int nproj = atom_nh[it];
#ifdef _OPENMP
#pragma omp for collapse(2)
#endif
                for (int ia = 0; ia < atom_na[it]; ia++)
                {
                    for (int ib = 0; ib < nbands_occ; ib++)
                    {
                        FPTYPE local_force[3] = {0, 0, 0};
                        FPTYPE fac;
                        if(occ)
                        {
                            fac = d_wg[ib] * 2.0 * tpiba;
                        }
                        else
                        {
                            fac = d_wg[0] * 2.0 * tpiba;
                        }
                        FPTYPE ekb_now = 0.0;
                        if (d_ekb != nullptr)
                        {
                            ekb_now = d_ekb[ib];
                        }
                        int iat = iat0 + ia;
                        int sum = sum0 + ia * nproj;
                        for (int ip = 0; ip < nproj; ip++)
                        {
                            FPTYPE ps_qq = 0;
                            if(ekb_now != 0)
                            {
                                ps_qq = - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip];
                            }
                            // Effective values of the D-eS coefficients
                            FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip] + ps_qq;
                            const int inkb = sum + ip;
                            // out<<"\n ps = "<<ps;

                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                const FPTYPE dbb
                                    = (conj(dbecp[ipol * nbands * nkb + ib * nkb + inkb]) * becp[ib * nkb + inkb])
                                          .real();
                                local_force[ipol] -= ps * fac * dbb;
                                // cf[iat*3+ipol] += ps * fac * dbb;
                            }
                            if (nondiagonal)
                            {
                                for (int ip2 = 0; ip2 < nproj; ip2++)
                                {
                                    if (ip != ip2)
                                    {
                                        const int jnkb = sum + ip2;
                                        FPTYPE ps_qq = 0;
                                        if (ekb_now != 0)
                                        {
                                            ps_qq = -ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip2];
                                        }
                                        FPTYPE ps = deeq[((spin * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;
                                        for (int ipol = 0; ipol < 3; ipol++)
                                        {
                                            const FPTYPE dbb = (conj(dbecp[ipol * nbands * nkb + ib * nkb + inkb])
                                                                * becp[ib * nkb + jnkb])
                                                                   .real();
                                            local_force[ipol] -= ps * fac * dbb;
                                        }
                                    }
                                }
                            }
                        }
#ifdef _OPENMP
                        if (omp_get_num_threads() > 1)
                        {
                            for (int ipol = 0; ipol < 3; ipol++)
                            {
#pragma omp atomic
                                force[iat * forcenl_nc + ipol] += local_force[ipol];
                            }
                        }
                        else
#endif
                        {
                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                force[iat * forcenl_nc + ipol] += local_force[ipol];
                            }
                        }
                    }
                } // end ia
                iat0 += atom_na[it];
                sum0 += atom_na[it] * nproj;
            } // end it
#ifdef _OPENMP
        }
#endif
    };

    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nbands_occ,
                    const int& ntype,
                    const int& deeq_2,
                    const int& deeq_3,
                    const int& deeq_4,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const bool& occ,
                    const FPTYPE* d_ekb,
                    const FPTYPE* qq_nt,
                    const std::complex<FPTYPE>* deeq_nc,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force)
    {
#ifdef _OPENMP
#pragma omp parallel
        {
#endif
            int iat0 = 0;
            int sum0 = 0;
            for (int it = 0; it < ntype; it++)
            {
                const int nprojs = atom_nh[it];
#ifdef _OPENMP
#pragma omp for collapse(2)
#endif
                for (int ia = 0; ia < atom_na[it]; ia++)
                {
                    for (int ib = 0; ib < nbands_occ; ib++)
                    {
                        const int ib2 = ib*2;
                        FPTYPE local_force[3] = {0, 0, 0};
                        FPTYPE fac;
                        if(occ)
                        {
                            fac = d_wg[ib] * 2.0 * tpiba;
                        }
                        else
                        {
                            fac = d_wg[0] * 2.0 * tpiba;
                        }
                        FPTYPE ekb_now = 0.0;
                        if (d_ekb != nullptr)
                        {
                            ekb_now = d_ekb[ib];
                        }
                        int iat = iat0 + ia;
                        int sum = sum0 + ia * nprojs;
                        for (int ip = 0; ip < nprojs; ip++)
                        {
                            const int inkb = sum + ip;
                            // out<<"\n ps = "<<ps;
                            for (int ip2 = 0; ip2 < nprojs; ip2++)
                            {
                                // Effective values of the D-eS coefficients
                                std::complex<FPTYPE> ps_qq = 0;
                                if(ekb_now)
                                {
                                    ps_qq = - ekb_now * qq_nt[it * deeq_3 * deeq_4 + ip * deeq_4 + ip2];
                                }
                                const int jnkb = sum + ip2;
                                std::complex<FPTYPE> ps0 = deeq_nc[((0 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;
                                std::complex<FPTYPE> ps1 = deeq_nc[((1 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2];
                                std::complex<FPTYPE> ps2 = deeq_nc[((2 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2];
                                std::complex<FPTYPE> ps3 = deeq_nc[((3 * deeq_2 + iat) * deeq_3 + ip) * deeq_4 + ip2] + ps_qq;
                                

                                for (int ipol = 0; ipol < 3; ipol++)
                                {
                                    const int index0 = ipol * nbands * 2 * nkb + ib2 * nkb + inkb;
                                    const int index1 = ib2 * nkb + jnkb;
                                    const std::complex<FPTYPE> dbb0 = conj(dbecp[index0]) * becp[index1];
                                    const std::complex<FPTYPE> dbb1 = conj(dbecp[index0]) * becp[index1 + nkb];
                                    const std::complex<FPTYPE> dbb2 = conj(dbecp[index0 + nkb]) * becp[index1];
                                    const std::complex<FPTYPE> dbb3 = conj(dbecp[index0 + nkb]) * becp[index1 + nkb];

                                    local_force[ipol] -= fac * (ps0 * dbb0 + ps1 * dbb1 + ps2 * dbb2 + ps3 * dbb3).real();
                                }
                            }
                        }
#ifdef _OPENMP
                        if (omp_get_num_threads() > 1)
                        {
                            for (int ipol = 0; ipol < 3; ++ipol)
                            {
#pragma omp atomic
                                force[iat * forcenl_nc + ipol] += local_force[ipol];
                            }
                        }
                        else
#endif
                        {
                            for (int ipol = 0; ipol < 3; ++ipol)
                            {
                                force[iat * forcenl_nc + ipol] += local_force[ipol];
                            }
                        }
                    }
                } // end ia
                iat0 += atom_na[it];
                sum0 += atom_na[it] * nprojs;
            } // end it
#ifdef _OPENMP
        }
#endif
    };
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const std::complex<FPTYPE>* vu,
                    const int* orbital_corr,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force)
    {
        int iat0 = 0;
        int sum0 = 0;
        for (int it = 0; it < ntype; it++)
        {
            const int orbital_l = orbital_corr[it];
            const int nproj = atom_nh[it];
            if(orbital_l == -1)
            {
                sum0 += nproj * atom_na[it];
                continue;
            }
            const int ip_begin = orbital_l * orbital_l;
            const int ip_end = (orbital_l + 1) * (orbital_l + 1);
            const int tlp1 = 2 * orbital_l + 1;
            const int tlp1_2 = tlp1 * tlp1;
            for (int ia = 0; ia < atom_na[it]; ia++)
            {
                for (int ib = 0; ib < nbands_occ; ib++)
                {
                    const int ib2 = ib*2;
                    FPTYPE local_force[3] = {0, 0, 0};
                    FPTYPE fac = d_wg[ik * wg_nc + ib] * 2.0 * tpiba;
                    int iat = iat0 + ia;
                    int sum = sum0 + ia * nproj;
                    for (int ip = ip_begin; ip < ip_end; ip++)
                    {
                        const int inkb = sum + ip;
                        const int m = ip - ip_begin;
                        // out<<"\n ps = "<<ps;
                        for (int ip2 = ip_begin; ip2 < ip_end; ip2++)
                        {
                            const int jnkb = sum + ip2;
                            const int m2 = ip2 - ip_begin;
                            std::complex<FPTYPE> ps[4];
                            for(int i = 0; i < 4; i++)
                            {
                                ps[i] = vu[(i * tlp1_2 + m * tlp1 + m2)];
                            }

                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                const int index0 = ipol * nbands * 2 * nkb + ib2 * nkb + inkb;
                                const int index1 = ib2 * nkb + jnkb;
                                const std::complex<FPTYPE> dbb0 = conj(dbecp[index0]) * becp[index1];
                                const std::complex<FPTYPE> dbb1 = conj(dbecp[index0]) * becp[index1 + nkb];
                                const std::complex<FPTYPE> dbb2 = conj(dbecp[index0 + nkb]) * becp[index1];
                                const std::complex<FPTYPE> dbb3 = conj(dbecp[index0 + nkb]) * becp[index1 + nkb];

                                local_force[ipol] -= fac * (ps[0] * dbb0 + ps[1] * dbb1 + ps[2] * dbb2 + ps[3] * dbb3).real();
                            }
                        }
                    }
                    for (int ipol = 0; ipol < 3; ++ipol)
                    {
                        force[iat * forcenl_nc + ipol] += local_force[ipol];
                    }
                }
                vu += 4 * tlp1_2;// step for vu
            } // end ia
            iat0 += atom_na[it];
            sum0 += atom_na[it] * nproj;
        } // end it
    };

    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nbands_occ,
                    const int& wg_nc,
                    const int& ntype,
                    const int& forcenl_nc,
                    const int& nbands,
                    const int& ik,
                    const int& nkb,
                    const int* atom_nh,
                    const int* atom_na,
                    const FPTYPE& tpiba,
                    const FPTYPE* d_wg,
                    const FPTYPE* lambda,
                    const std::complex<FPTYPE>* becp,
                    const std::complex<FPTYPE>* dbecp,
                    FPTYPE* force)
    {
        int iat0 = 0;
        int sum0 = 0;
        for (int it = 0; it < ntype; it++)
        {
            const int nproj = atom_nh[it];
            for (int ia = 0; ia < atom_na[it]; ia++)
            {
                int iat = iat0 + ia;
                int sum = sum0 + ia * nproj;
                const std::complex<FPTYPE> coefficients0(lambda[iat*3+2], 0.0);
                const std::complex<FPTYPE> coefficients1(lambda[iat*3] , lambda[iat*3+1]);
                const std::complex<FPTYPE> coefficients2(lambda[iat*3] , -1 * lambda[iat*3+1]);
                const std::complex<FPTYPE> coefficients3(-1 * lambda[iat*3+2], 0.0);
                for (int ib = 0; ib < nbands_occ; ib++)
                {
                    const int ib2 = ib*2;
                    FPTYPE local_force[3] = {0, 0, 0};
                    FPTYPE fac = d_wg[ik * wg_nc + ib] * 2.0 * tpiba;
                    for (int ip = 0; ip < nproj; ip++)
                    {
                        const int inkb = sum + ip;

                        for (int ipol = 0; ipol < 3; ipol++)
                        {
                            const int index0 = ipol * nbands * 2 * nkb + ib2 * nkb + inkb;
                            const int index1 = ib2 * nkb + inkb;
                            const std::complex<FPTYPE> dbb0 = conj(dbecp[index0]) * becp[index1];
                            const std::complex<FPTYPE> dbb1 = conj(dbecp[index0]) * becp[index1 + nkb];
                            const std::complex<FPTYPE> dbb2 = conj(dbecp[index0 + nkb]) * becp[index1];
                            const std::complex<FPTYPE> dbb3 = conj(dbecp[index0 + nkb]) * becp[index1 + nkb];

                            local_force[ipol] -= fac * (coefficients0 * dbb0 + coefficients1 * dbb1 + coefficients2 * dbb2 + coefficients3 * dbb3).real();
                        }
                    }//ip
                    for (int ipol = 0; ipol < 3; ++ipol)
                    {
                        force[iat * forcenl_nc + ipol] += local_force[ipol];
                    }
                } // end ib
            } // ia
            iat0 += atom_na[it];
            sum0 += atom_na[it] * nproj;
        }//it
    }
};

// CPU implementation of local force sincos operator
template <typename FPTYPE>
struct cal_force_loc_sincos_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ntype,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* vloc_factors,
                    const std::complex<FPTYPE>* aux,
                    FPTYPE* force)
    {
        const FPTYPE TWO_PI = 2.0 * M_PI;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int iat = 0; iat < nat; ++iat)
        {
            const int it = iat2it[iat];
            const FPTYPE tau_x = tau[iat * 3 + 0];
            const FPTYPE tau_y = tau[iat * 3 + 1];
            const FPTYPE tau_z = tau[iat * 3 + 2];

            FPTYPE local_force[3] = {0.0, 0.0, 0.0};

            for (int ig = 0; ig < npw; ig++)
            {
                const FPTYPE phase = TWO_PI * (gcar[ig * 3 + 0] * tau_x + 
                                               gcar[ig * 3 + 1] * tau_y + 
                                               gcar[ig * 3 + 2] * tau_z);
                FPTYPE sinp, cosp;
                ModuleBase::libm::sincos(phase, &sinp, &cosp);
                
                const FPTYPE vloc_factor = vloc_factors[it * npw + ig];
                const FPTYPE factor = vloc_factor * (cosp * aux[ig].imag() + sinp * aux[ig].real());
                
                local_force[0] += gcar[ig * 3 + 0] * factor;
                local_force[1] += gcar[ig * 3 + 1] * factor;
                local_force[2] += gcar[ig * 3 + 2] * factor;
            }

            force[iat * 3 + 0] = local_force[0];
            force[iat * 3 + 1] = local_force[1];
            force[iat * 3 + 2] = local_force[2];
        }
    }
};

// CPU implementation of Ewald parallel operator
template <typename FPTYPE>
struct cal_force_ew_parallel_op<FPTYPE, base_device::DEVICE_CPU>
{
    void operator()(const base_device::DEVICE_CPU* ctx,
                    const int& nat,
                    const int& npw,
                    const int& ntype,
                    const int& ig_gge0,
                    const FPTYPE* gcar,
                    const FPTYPE* tau,
                    const int* iat2it,
                    const FPTYPE* it_facts,
                    const int* atom_na,
                    const FPTYPE* zv_values,
                    const std::complex<FPTYPE>* aux,
                    const FPTYPE& alpha,
                    const FPTYPE& fact,
                    const FPTYPE& lat0,
                    const FPTYPE* latvec,
                    const FPTYPE* G,
                    const bool& use_paw,
                    FPTYPE* force)
    {
        const FPTYPE TWO_PI = 2.0 * M_PI;
        const FPTYPE e2 = 2.0; // ModuleBase::e2 equivalent

#ifdef _OPENMP
#pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            int thread_id = omp_get_thread_num();
#else
        int num_threads = 1;
        int thread_id = 0;
#endif

            // Task distribution for multi-thread
            int iat_beg, iat_end;
            ModuleBase::TASK_DIST_1D(num_threads, thread_id, nat, iat_beg, iat_end);
            iat_end = iat_beg + iat_end;

            // Preprocess ig_gap for skipping the ig point
            int ig_gap = (ig_gge0 >= 0 && ig_gge0 < npw) ? ig_gge0 : -1;

            // Part 1: sincos calculation (replacing original while loop with for loop)
            for (int iat = iat_beg; iat < iat_end; ++iat)
            {
                const int it = iat2it[iat];
                
                // Check if this atom type has atoms
                if (atom_na[it] == 0) continue;

                const FPTYPE tau_x = tau[iat * 3 + 0];
                const FPTYPE tau_y = tau[iat * 3 + 1];
                const FPTYPE tau_z = tau[iat * 3 + 2];
                const FPTYPE it_fact = it_facts[iat];

                FPTYPE local_force[3] = {0.0, 0.0, 0.0};

                // G-vector loop with ig_gap handling (using for loops for GPU compatibility)
                for (int ig = 0; ig < npw; ig++)
                {
                    // Skip G=0 term
                    if (ig == ig_gap) continue;

                    const FPTYPE arg = TWO_PI * (gcar[ig * 3 + 0] * tau_x + 
                                                 gcar[ig * 3 + 1] * tau_y + 
                                                 gcar[ig * 3 + 2] * tau_z);
                    FPTYPE sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    
                    const FPTYPE sumnb = -cosp * aux[ig].imag() + sinp * aux[ig].real();
                    
                    local_force[0] += gcar[ig * 3 + 0] * sumnb;
                    local_force[1] += gcar[ig * 3 + 1] * sumnb;
                    local_force[2] += gcar[ig * 3 + 2] * sumnb;
                }

                // Apply it_fact scaling
                force[iat * 3 + 0] += local_force[0] * it_fact;
                force[iat * 3 + 1] += local_force[1] * it_fact;
                force[iat * 3 + 2] += local_force[2] * it_fact;
            }

            // Part 2: Real space Ewald calculation (when processor contains G=0 term)
            if (ig_gge0 >= 0)
            {
                const FPTYPE rmax = 5.0 / (sqrt(alpha) * lat0);
                const int mxr = 200;
                
                // Allocate temporary arrays for each thread
                ModuleBase::Vector3<double>* r_flat = new ModuleBase::Vector3<double>[mxr];
                double* r2 = new double[mxr];
                int* irr = new int[mxr];
                
                // Initialize arrays
                for (int i = 0; i < mxr; i++) {
                    r2[i] = 0.0;
                    irr[i] = 0;
                    r_flat[i].x = 0.0;
                    r_flat[i].y = 0.0;
                    r_flat[i].z = 0.0;
                }

                const FPTYPE sqa = sqrt(alpha);
                const FPTYPE sq8a_2pi = sqrt(8.0 * alpha / TWO_PI);

                // Iterate over atoms in assigned range
                for (int iat1 = iat_beg; iat1 < iat_end; ++iat1)
                {
                    const int T1 = iat2it[iat1];
                    if (atom_na[T1] == 0) continue;

                    const FPTYPE tau1_x = tau[iat1 * 3 + 0];
                    const FPTYPE tau1_y = tau[iat1 * 3 + 1];
                    const FPTYPE tau1_z = tau[iat1 * 3 + 2];

                    // Iterate over all other atoms
                    for (int iat2 = 0; iat2 < nat; ++iat2)
                    {
                        const int T2 = iat2it[iat2];
                        
                        if (iat1 != iat2 && atom_na[T2] != 0 && atom_na[T1] != 0)
                        {
                            const FPTYPE tau2_x = tau[iat2 * 3 + 0];
                            const FPTYPE tau2_y = tau[iat2 * 3 + 1];
                            const FPTYPE tau2_z = tau[iat2 * 3 + 2];
                            
                            // Calculate d_tau
                            const FPTYPE d_tau_x = tau1_x - tau2_x;
                            const FPTYPE d_tau_y = tau1_y - tau2_y;
                            const FPTYPE d_tau_z = tau1_z - tau2_z;

                            // Reconstruct Matrix3 objects from flat arrays for rgen call
                            ModuleBase::Matrix3 latvec_matrix(
                                static_cast<double>(latvec[0]), static_cast<double>(latvec[1]), static_cast<double>(latvec[2]),
                                static_cast<double>(latvec[3]), static_cast<double>(latvec[4]), static_cast<double>(latvec[5]),
                                static_cast<double>(latvec[6]), static_cast<double>(latvec[7]), static_cast<double>(latvec[8])
                            );
                            
                            ModuleBase::Matrix3 G_matrix(
                                static_cast<double>(G[0]), static_cast<double>(G[1]), static_cast<double>(G[2]),
                                static_cast<double>(G[3]), static_cast<double>(G[4]), static_cast<double>(G[5]),
                                static_cast<double>(G[6]), static_cast<double>(G[7]), static_cast<double>(G[8])
                            );
                            
                            ModuleBase::Vector3<double> dtau_vec(
                                static_cast<double>(d_tau_x),
                                static_cast<double>(d_tau_y),
                                static_cast<double>(d_tau_z)
                            );

                            // Call H_Ewald_pw::rgen to generate neighbor shells
                            int nrm = 0;
                            H_Ewald_pw::rgen(dtau_vec, static_cast<double>(rmax), irr, latvec_matrix, G_matrix, r_flat, r2, nrm);
                            
                            // Process all generated r vectors
                            for (int nr = 0; nr < nrm; nr++)
                            {
                                const FPTYPE rr = sqrt(r2[nr]) * lat0;
                                
                                FPTYPE factor;
                                if (use_paw)
                                {
                                    factor = zv_values[T1] * zv_values[T2] * e2 / (rr * rr)
                                             * (erfc(sqa * rr) / rr + sq8a_2pi * exp(-alpha * rr * rr))
                                             * lat0;
                                }
                                else
                                {
                                    factor = zv_values[T1] * zv_values[T2] * e2 / (rr * rr)
                                             * (erfc(sqa * rr) / rr + sq8a_2pi * exp(-alpha * rr * rr))
                                             * lat0;
                                }

                                force[iat1 * 3 + 0] -= factor * static_cast<FPTYPE>(r_flat[nr].x);
                                force[iat1 * 3 + 1] -= factor * static_cast<FPTYPE>(r_flat[nr].y);
                                force[iat1 * 3 + 2] -= factor * static_cast<FPTYPE>(r_flat[nr].z);
                            }
                        }
                    }
                }

                delete[] r_flat;
                delete[] r2;
                delete[] irr;
            }

#ifdef _OPENMP
        }
#endif
    }
};

template struct cal_vkb1_nl_op<float, base_device::DEVICE_CPU>;
template struct cal_force_nl_op<float, base_device::DEVICE_CPU>;
template struct cal_force_loc_sincos_op<float, base_device::DEVICE_CPU>;
template struct cal_force_ew_parallel_op<float, base_device::DEVICE_CPU>;

template struct cal_vkb1_nl_op<double, base_device::DEVICE_CPU>;
template struct cal_force_nl_op<double, base_device::DEVICE_CPU>;
template struct cal_force_loc_sincos_op<double, base_device::DEVICE_CPU>;
template struct cal_force_ew_parallel_op<double, base_device::DEVICE_CPU>;

} // namespace hamilt
