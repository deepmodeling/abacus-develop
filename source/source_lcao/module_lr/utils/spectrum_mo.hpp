#pragma once
#include "source_base/tool_title.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/klist.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_lr/ao_to_mo_transformer/ao_to_mo.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "source_lcao/module_rt/velocity_op.h"
#include "source_psi/psi.h"

namespace LR_Util
{
template<typename T>
std::vector<std::complex<double>> cal_velocity_mo(const UnitCell& ucell,
                                                const Grid_Driver& gd,
                                                const TwoCenterBundle& two_center_bundle,
                                                const Parallel_Orbitals& pmat,/*nbasis×nbasis*/
                                                const Parallel_2D& pc,/*nbasis×nbands*/
                                                const K_Vectors& kv,
                                                const psi::Psi<T>& psi_ks,
                                                const int nk,
                                                const int nspin_tmp,
                                                const int nbasis, 
                                                const std::vector<int> nocc, 
                                                const std::vector<int> nvirt)
{
    ModuleBase::TITLE("LR_Util", "cal_velocity_mo");
    ModuleBase::timer::start("LR_Util", "cal_velocity_mo");
    std::cout<<"Calculating velocity matrix in KS presentation..."<<std::endl;
    // get_velocity_matrix_R(ucell, gd_, pmat, two_center_bundle_);
    LCAO_Orbitals orb;
    const auto& inp = PARAM.inp;
    two_center_bundle.to_LCAO_Orbitals(orb, inp.lcao_ecut, inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax,
                                        inp.out_element_info, inp.cal_force);

    Velocity_op<std::complex<double>> vR(&ucell, &gd, &pmat, orb, two_center_bundle.overlap_orb.get());
    vR.calculate_grad_term();   // $<\mu, 0|-i∇r|\nu, R>$
    vR.calculate_vcomm_r(); // $<\mu, 0|i[Vnl, r]|\nu, R>$

    int nks = kv.get_nks(); // include spin
    assert(nks == nk * nspin_tmp);
    assert(psi_ks.get_nk() == nks);
    int KS_num = nocc[0] + nvirt[0];
    std::vector<std::complex<double>> velocity_mo(nspin_tmp * 3 * nk * KS_num * KS_num, 0.0);
    Parallel_2D pmo;
    LR_Util::setup_2d_division(pmo, pmat.get_block_size(), KS_num, KS_num
#ifdef __MPI
        , pc.blacs_ctxt
#endif
        );

    //1. psi_ks<T> to c_psi_ks<complex<double>>, ensure complex<double> for dipole calculation
    psi::Psi<std::complex<double>> c_psi_ks(nks,
                                            pc.get_col_size(), 
                                            pc.get_row_size(), 
                                            kv.ngk, 
                                            true);
    for(int iks = 0; iks < nks; ++iks)
    {
        for(int ic = 0; ic < pc.get_col_size(); ++ic)//band
        {
            for(int ir = 0; ir < pc.get_row_size(); ++ir)//basis
            {
                c_psi_ks(iks, ic, ir) = std::complex<double>(psi_ks(iks, ic, ir));
            }
        }
    }
    
    //2. calculate v_mo = c^\dagger v c
    std::vector<ct::Tensor> vk(nks, LR_Util::newTensor<std::complex<double>>({ pmat.get_col_size(), pmat.get_row_size() }));
    for (int id = 0; id < 3; ++id)
    {
        for (auto& v : vk) v.zero();

        std::vector<std::complex<double>> v_mo(nks * pmo.get_local_size(), 0.0);

        for (int is = 0; is < nspin_tmp; ++is)
        {
            assert(KS_num == nocc[is] + nvirt[is]);

            for (int ik = is*nk; ik < (is+1)*nk; ++ik)
            {            
                hamilt::folding_HR(*vR.get_current_term_pointer(id), vk[ik].data<std::complex<double>>(), kv.kvec_d[ik], pmat.get_row_size(), 1/*column-major*/);
            }
        }
#ifdef __MPI
        LR::ao_to_mo_pblas(vk, pmat, c_psi_ks, pc, nbasis,
                        nocc[0], nvirt[0], pmo, v_mo.data(),
                        false, // add_on
                        LR_Util::MO_TYPE::ALL);
#else
        LR::ao_to_mo_blas(vk, c_psi_ks, 
                        nocc[0], nvirt[0], v_mo.data(),
                        false , //add_on
                        LR_Util::MO_TYPE::ALL);
#endif
        // gather local vk to global velocity_mo
        for (int is = 0; is < nspin_tmp; ++is)
        {
            for (int ik = 0; ik < nk; ++ik)
            {
                int glb_offset = (is * 3 * nk + id * nk + ik) * KS_num * KS_num;
                int loc_offset = ik * pmo.get_local_size();
                for (int j = 0; j < pmo.get_col_size(); ++j){
                    for (int i = 0; i < pmo.get_row_size(); ++i){
                        velocity_mo[glb_offset + pmo.local2global_col(j) * KS_num + pmo.local2global_row(i)]
                            = v_mo[loc_offset + j * pmo.get_row_size() + i];
                    }
                }
            }
        }
    }//id
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, velocity_mo.data(), velocity_mo.size(), LR_Util::MPIType<T>::value(), MPI_SUM, pmo.comm());
#endif
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "Finish velocity matrix in KS presentation.");
    ModuleBase::timer::end("LR_Util", "cal_velocity_mo");
    return velocity_mo;
}

/// @brief output the velocity matrix in KS presentation in human-read friendly format
inline void output_spectrum_mo(const std::vector<std::complex<double>>& out_spectrum_mo,
                        const std::string& filename,
                        const double* const eig_ks,
                        const int nk,
                        const int nspin_tmp,
                        const int KS_num,
                        const K_Vectors& kv)
{
    assert(out_spectrum_mo.size() == nspin_tmp * 3 * nk * KS_num * KS_num);
    std::ofstream ofs(filename);
    ofs << "Data Unit (Hartree * Bohr)" << std::endl;
    ofs << "NOTICE: KS_index are restricted in nocc and nvirt" << std::endl;

    int step = nk * KS_num * KS_num;
    int offset = 0, ipair = 0;
    for (int is = 0; is < nspin_tmp; ++is)
    {
        ofs << "ispin: " << is << std::endl;
        for (int ik = 0;ik < nk;++ik)
        {
            ofs << "k-point: " << ik << " " << kv.kvec_d[ik] << std::endl;
            ofs << std::setw(4) << "KS1" << std::setw(10) << "E1(eV)" << std::setw(6) << "KS2" << std::setw(10) << "E2(eV)"
                << std::setw(16) << "x" << std::setw(23) << "|x|^2" << std::setw(18) << "y" << std::setw(23) <<"|y|^2"
                << std::setw(18) << "z" << std::setw(23) <<"|z|^2" << std::setw(13) << "average" << std::endl;
            
            offset = (is * 3 * nk + ik) * KS_num * KS_num;
            for (int i = 0; i < KS_num; ++i)
            {
                for (int j = i; j < KS_num ; ++j)
                {
                    ipair = offset + i * KS_num + j;
                    double average = (std::norm(out_spectrum_mo[ipair]) + std::norm(out_spectrum_mo[ipair + step]) + std::norm(out_spectrum_mo[ipair + 2 * step])) / 3.0;
                    ofs << std::setw(3) << i << std::setw(12) << std::setprecision(6) << eig_ks[i + ik*KS_num] * ModuleBase::Ry_to_eV
                    << std::setw(4) << j << std::setw(12) << eig_ks[j + ik*KS_num] * ModuleBase::Ry_to_eV
                    << std::setw(28) << out_spectrum_mo[ipair] << std::setw(13) << std::norm(out_spectrum_mo[ipair])
                    << std::setw(28) << out_spectrum_mo[ipair + step] << std::setw(13) << std::norm(out_spectrum_mo[ipair + step])
                    << std::setw(28) << out_spectrum_mo[ipair + 2 * step]  << std::setw(13) << std::norm(out_spectrum_mo[ipair + 2 * step] )
                    << std::setw(13) << average << std::endl; 
                }
            }
        }
    }
    ofs.close();
}

/// @brief output the velocity matrix in KS presentation in LibRPA format
inline void output_spectrum_mo_librpa(const std::vector<std::complex<double>>& out_spectrum_mo,
                        const std::string& filename,
                        const int nk, const int nspin_tmp, const int nbands, const int nbasis,
                        const int KS_num,/*= nocc+nvirt, can be different from nbands, */
                        const K_Vectors& kv)
{
    assert(out_spectrum_mo.size() == nspin_tmp * 3 * nk * KS_num * KS_num);
    std::ofstream ofs(filename);
    ofs << std::scientific << nk << std::endl;
    ofs << nspin_tmp << std::endl;
    ofs << nbands << std::endl;
    ofs << nbasis << std::endl;
    double HaBohrToEvAng = ModuleBase::Hartree_to_eV * ModuleBase::BOHR_TO_A; // 27.211396 * 0.5291770
    int offset = 0, ipair = 0;
    for (int is = 0; is < nspin_tmp; ++is)
    {
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int id = 0; id < 3; ++id)
            {
                offset = (is * 3 * nk + id * nk + ik) * KS_num * KS_num;
                ofs <<"  " << id + 1 << "  " << ik + 1 << "  " << is + 1 << std::endl;
                for (int ib1 = 0; ib1 < KS_num; ++ib1)
                {
                    for (int ib2 = 0; ib2 < KS_num; ++ib2)
                    {
                        ipair = offset + ib1 * KS_num + ib2;
                        ofs << std::setw(26) << std::fixed << std::setprecision(16) << out_spectrum_mo[ipair].real() * HaBohrToEvAng
                            << std::setw(26) << out_spectrum_mo[ipair].imag() * HaBohrToEvAng << std::endl;
                    }
                }
            }

        }
    }
}
}