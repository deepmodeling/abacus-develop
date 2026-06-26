#include "force_rt_overlap.h"
#include "td_info.h"
#include "td_folding.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_base/parallel_reduce.h"
template <>
void cal_foverlap_rt(ModuleBase::matrix& foverlap,
                     const LCAO_domain::Setup_DM<std::complex<double>>& dmat,
                     hamilt::Hamilt<std::complex<double>>* p_hamilt,
                     const K_Vectors& kv,
                     Parallel_Orbitals& pv,
                     UnitCell& ucell)
{
#ifdef __MPI
    const int nlocal = PARAM.globalv.nlocal;
    assert(nlocal >= 0);

    TD_info* td_info = TD_info::td_vel_op;
    //get dS/dR_{x,y,z}
    std::vector<hamilt::HContainer<double>*> dsxr = td_info->get_grad_overlap();
    // allocate matrix
    const long nloc = pv.nloc;
    const int nrow = pv.nrow;
    std::complex<double>* Htmp = new std::complex<double>[nloc];
    std::complex<double>* Sinv = new std::complex<double>[nloc];
    std::complex<double>* dsxk = new std::complex<double>[nloc];
    std::complex<double>* pdsxk = new std::complex<double>[nloc];
    std::complex<double>* tmp1 = new std::complex<double>[nloc];
    std::complex<double>* tmp2 = new std::complex<double>[nloc];
    std::complex<double>* tmp3 = new std::complex<double>[nloc];
    std::complex<double>* Hybridtmp = new std::complex<double>[nloc];
    std::vector<std::complex<double>*> tmp_out = {nullptr, nullptr, nullptr};
    for(int dir = 0; dir<3; dir++)
    {
        tmp_out[dir] = new std::complex<double>[nloc];
        ModuleBase::GlobalFunc::ZEROS(tmp_out[dir], nloc);
    }
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        p_hamilt->updateHk(ik);
        // get dmk
        std::complex<double>* tmp_dmk = dmat.dm->get_DMK_pointer(ik);

        ModuleBase::GlobalFunc::ZEROS(Htmp, nloc);
        ModuleBase::GlobalFunc::ZEROS(Sinv, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp1, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp2, nloc);
        ModuleBase::GlobalFunc::ZEROS(Hybridtmp, nloc);
        
        const int inc = 1;

        hamilt::MatrixBlock<std::complex<double>> h_mat;
        hamilt::MatrixBlock<std::complex<double>> s_mat;
        //get Hk Sk
        p_hamilt->matrix(h_mat, s_mat);
        BlasConnector::copy(nloc, h_mat.p, inc, Htmp, inc);
        BlasConnector::copy(nloc, s_mat.p, inc, Sinv, inc);

        vector<int> ipiv(nloc, 0);
        int info = 0;
        const int one_int = 1;
        const std::complex<double> one_complex = {1.0, 0.0};
        const std::complex<double> two_complex = {2.0, 0.0};
        const std::complex<double> mone_complex = {-1.0, 0.0};
        const std::complex<double> zero_complex = {0.0, 0.0};

        ScalapackConnector::getrf(nlocal, nlocal, Sinv, one_int, one_int, pv.desc, ipiv.data(), &info);

        int lwork = -1;
        int liwork = -1;

        // if lwork == -1, then the size of work is (at least) of length 1.
        std::vector<std::complex<double>> work(1, 0);

        // if liwork = -1, then the size of iwork is (at least) of length 1.
        std::vector<int> iwork(1, 0);

        ScalapackConnector::getri(nlocal,
                                  Sinv,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  ipiv.data(),
                                  work.data(),
                                  &lwork,
                                  iwork.data(),
                                  &liwork,
                                  &info);

        lwork = work[0].real();
        work.resize(lwork, 0);
        liwork = iwork[0];
        iwork.resize(liwork, 0);

        ScalapackConnector::getri(nlocal,
                                  Sinv,
                                  one_int,
                                  one_int,
                                  pv.desc,
                                  ipiv.data(),
                                  work.data(),
                                  &lwork,
                                  iwork.data(),
                                  &liwork,
                                  &info);

        const char N_char = 'N';
        const char T_char = 'T';
        const char C_char = 'C';
        

        ScalapackConnector::gemm(T_char,
                                 C_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 tmp_dmk,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 Htmp,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp1,
                                 one_int,
                                 one_int,
                                 pv.desc);

        ScalapackConnector::gemm(N_char,
                                 N_char,
                                 nlocal,
                                 nlocal,
                                 nlocal,
                                 one_complex,
                                 tmp1,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 Sinv,
                                 one_int,
                                 one_int,
                                 pv.desc,
                                 zero_complex,
                                 tmp2,
                                 one_int,
                                 one_int,
                                 pv.desc);
        for(int dir = 0; dir<3; dir++)
        {
            ModuleBase::GlobalFunc::ZEROS(dsxk, nloc);
            ModuleBase::GlobalFunc::ZEROS(pdsxk, nloc);
            ModuleBase::GlobalFunc::ZEROS(tmp3, nloc);
            module_rt::folding_HR_td(*dsxr[dir], dsxk, kv.kvec_d[ik], TD_info::cart_At, TD_info::td_vel_op->get_phase_hybrid(), nrow, 1);
            module_rt::folding_partial_dot(*dsxr[dir], pdsxk, kv.kvec_d[ik], nrow, 1, &ucell, TD_info::td_vel_op->get_phase_hybrid(), TD_info::cart_At, elecstate::H_TDDFT_pw::Et);
            ScalapackConnector::gemm(N_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     two_complex,
                                     tmp2,
                                     one_int,
                                     one_int,
                                     pv.desc,
                                     dsxk,
                                     one_int,
                                     one_int,
                                     pv.desc,
                                     one_complex,
                                     tmp_out[dir],
                                     one_int,
                                     one_int,
                                     pv.desc);

            ScalapackConnector::gemm(T_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     mone_complex,
                                     tmp_dmk,
                                     one_int,
                                     one_int,
                                     pv.desc,
                                     pdsxk,
                                     one_int,
                                     one_int,
                                     pv.desc,
                                     one_complex,
                                     tmp3,
                                     one_int,
                                     one_int,
                                     pv.desc);
            
            ScalapackConnector::geadd(N_char,
                                      nlocal,
                                      nlocal,
                                      one_complex,
                                      tmp3,
                                      one_int,
                                      one_int,
                                      pv.desc,
                                      one_complex,
                                      tmp_out[dir],
                                      one_int,
                                      one_int,
                                      pv.desc);
        }
    }
    delete[] Htmp;
    delete[] Sinv;
    delete[] tmp1;
    delete[] tmp2;
    delete[] dsxk;
    delete[] Hybridtmp;
    // std::string filename = "process_debug_" + std::to_string(GlobalV::MY_RANK) + ".txt";
    // std::ofstream debug_file(filename);
    // debug_file << "=== Process " << GlobalV::MY_RANK << " ===" << std::endl;
    // debug_file << "=== Matrix Partition Information ===" << std::endl;
    // auto row_indexes = pv.get_indexes_row();
    // auto col_indexes = pv.get_indexes_col();
    // for(int iat = 0; iat < ucell.nat; iat++)
    // {
    //     int row0 = pv.atom_begin_row[iat];
    //     int col0 = pv.atom_begin_col[iat];
    //     for(int mu = 0; mu < pv.get_row_size(iat); ++mu)
    //     {
    //         for(int nu = 0; nu < pv.get_col_size(iat); ++nu)
    //         {
    //             debug_file<<"mu: "<<mu<<" nu: "<<nu<<std::endl;
    //             debug_file<<"globalmu: "<<row_indexes[row0+mu]<<" globalnu: "<<col_indexes[col0+nu]<<std::endl;
    //         }
    //     }
    // }
    auto row_indexes = pv.get_indexes_row();
    auto col_indexes = pv.get_indexes_col();
    #pragma omp parallel for
    for(int iat = 0; iat < ucell.nat; iat++)
    {
        
        double* force_tmp1 = &foverlap(iat, 0);
        int row0 = pv.atom_begin_row[iat];
        int col0 = pv.atom_begin_col[iat];
        const int row_size = pv.get_row_size();
        std::vector<std::complex<double>*> p_diag = {tmp_out[0], tmp_out[1], tmp_out[2]};
        for(int mu = 0; mu < pv.get_nrow_atom(iat); ++mu)
        {
            for(int nu = 0; nu < pv.get_ncol_atom(iat); ++nu)
            {
                if(row_indexes[row0+mu]==col_indexes[col0+nu])
                {
                    const long index = (row0 + mu) + (col0 + nu)*row_size;
                    for(int dir = 0; dir<3; dir++)
                    {
                        // p_diag[dir] = tmp_out[dir] + (row0 + mu) + (col0 + nu) * row_size;
                        force_tmp1[dir] += (tmp_out[dir][index]).real();
                    }
                }
            }
            
        }
    }
    Parallel_Reduce::reduce_all(foverlap.c, foverlap.nr * foverlap.nc);
    for(int dir = 0; dir<3; dir++)
    {
        delete[] tmp_out[dir];
    }
    return;
#endif
}
template <>
void cal_foverlap_rt(ModuleBase::matrix& foverlap,
                     const LCAO_domain::Setup_DM<double>& dmat,
                     hamilt::Hamilt<double>* p_hamilt,
                     const K_Vectors& kv,
                     Parallel_Orbitals& pv,
                     UnitCell& ucell)
{
    return;
}