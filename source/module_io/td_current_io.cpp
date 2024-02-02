#include "td_current_io.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/libm/libm.h"
#include "module_base/tool_threading.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_base/parallel_reduce.h"

#ifdef __LCAO
//init DSloc_R for current calculation
void ModuleIO::Init_DS_tmp(const Parallel_Orbitals& pv,
                              LCAO_Hamilt& UHM)
{    
    ModuleBase::TITLE("ModuleIO", "Init_DS_tmp");
    ModuleBase::timer::tick("ModuleIO", "Init_DS_tmp");
    const int nnr = pv.nnr;
    UHM.LM->DSloc_Rx = new double[nnr];
    UHM.LM->DSloc_Ry = new double[nnr];
    UHM.LM->DSloc_Rz = new double[nnr];
    const auto init_DSloc_Rxyz = [&UHM, nnr](int num_threads, int thread_id) {
        int beg, len;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, nnr, 1024, beg, len);
        ModuleBase::GlobalFunc::ZEROS(UHM.LM->DSloc_Rx + beg, len);
        ModuleBase::GlobalFunc::ZEROS(UHM.LM->DSloc_Ry + beg, len);
        ModuleBase::GlobalFunc::ZEROS(UHM.LM->DSloc_Rz + beg, len);
    };
    ModuleBase::OMP_PARALLEL(init_DSloc_Rxyz);
    bool cal_deri = true;
    UHM.genH.build_ST_new('S', cal_deri, GlobalC::ucell, UHM.genH.LM->SlocR.data());

    ModuleBase::timer::tick("ModuleIO", "Init_DS_tmp");
    return;
}
//destory DSloc_R so it can be used normally in the following force calculation
void ModuleIO::destory_DS_tmp(LCAO_Hamilt& UHM)
{
    ModuleBase::TITLE("ModuleIO", "destory_DS_tmp");
    ModuleBase::timer::tick("ModuleIO", "destory_DS_tmp");
    delete[] UHM.LM->DSloc_Rx;
    delete[] UHM.LM->DSloc_Ry;
    delete[] UHM.LM->DSloc_Rz;
    ModuleBase::timer::tick("ModuleIO", "destory_DS_tmp");
    return;
}
void ModuleIO::cal_tmp_DM(elecstate::DensityMatrix<std::complex<double>, double>& DM, const int ik, const int nspin)
{
    ModuleBase::TITLE("ModuleIO", "cal_tmp_DM");
    ModuleBase::timer::tick("ModuleIO", "cal_tmp_DM");
    int ld_hk = DM.get_paraV_pointer()->nrow;
    int ld_hk2 = 2 * ld_hk;
    for (int is = 1; is <= nspin; ++is)
    {
        int ik_begin = DM.get_DMK_nks() * (is - 1); // jump this->_nks for spin_down if nspin==2
        hamilt::HContainer<double>* tmp_DMR = DM.get_DMR_vector()[is - 1];
        // set zero since this function is called in every scf step
        tmp_DMR->set_zero();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int i = 0; i < tmp_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<double>& tmp_ap = tmp_DMR->get_atom_pair(i);
            int iat1 = tmp_ap.get_atom_i();
            int iat2 = tmp_ap.get_atom_j();
            // get global indexes of whole matrix for each atom in this process
            int row_ap = DM.get_paraV_pointer()->atom_begin_row[iat1];
            int col_ap = DM.get_paraV_pointer()->atom_begin_col[iat2];
            if (row_ap == -1 || col_ap == -1)
            {
                throw std::string("Atom-pair not belong this process");
            }
            for (int ir = 0; ir < tmp_ap.get_R_size(); ++ir)
            {
                const int* r_index = tmp_ap.get_R_index(ir);
                hamilt::BaseMatrix<double>* tmp_matrix = tmp_ap.find_matrix(r_index[0], r_index[1], r_index[2]);
#ifdef __DEBUG
                if (tmp_matrix == nullptr)
                {
                    std::cout << "tmp_matrix is nullptr" << std::endl;
                    continue;
                }
#endif
                // only ik
                if(GlobalV::NSPIN !=4 )
                {
                // cal k_phase
                // if TK==std::complex<double>, kphase is e^{ikR}
                const ModuleBase::Vector3<double> dR(r_index[0], r_index[1], r_index[2]);
                const double arg = (DM.get_kv_pointer()->kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                double sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                std::complex<double> kphase = std::complex<double>(cosp, sinp);
                // set DMR element
                double* tmp_DMR_pointer = tmp_matrix->get_pointer();
                std::complex<double>* tmp_DMK_pointer = DM.get_DMK_vector()[ik + ik_begin].data();
                double* DMK_real_pointer = nullptr;
                double* DMK_imag_pointer = nullptr;
                // jump DMK to fill DMR
                // DMR is row-major, DMK is column-major
                tmp_DMK_pointer += col_ap * DM.get_paraV_pointer()->nrow + row_ap;
                for (int mu = 0; mu < DM.get_paraV_pointer()->get_row_size(iat1); ++mu)
                {
                    DMK_real_pointer = (double*)tmp_DMK_pointer;
                    DMK_imag_pointer = DMK_real_pointer + 1;
                    BlasConnector::axpy(DM.get_paraV_pointer()->get_col_size(iat2),
                                        kphase.imag(),
                                        DMK_real_pointer,
                                        ld_hk2,
                                        tmp_DMR_pointer,
                                        1);
                    // "-" since i^2 = -1
                    BlasConnector::axpy(DM.get_paraV_pointer()->get_col_size(iat2),
                                        kphase.real(),
                                        DMK_imag_pointer,
                                        ld_hk2,
                                        tmp_DMR_pointer,
                                        1);
                    tmp_DMK_pointer += 1;
                    tmp_DMR_pointer += DM.get_paraV_pointer()->get_col_size(iat2);
                }
                }
            }
        }
    }
    ModuleBase::timer::tick("ModuleIO", "cal_tmp_DM");
}

void ModuleIO::write_current(const int istep,
                                const psi::Psi<std::complex<double>>* psi,
                                const elecstate::ElecState* pelec,
                                const K_Vectors& kv,
                                const Parallel_Orbitals* pv,
                                Record_adj& ra,
                                LCAO_Hamilt& UHM)
{

    ModuleBase::TITLE("ModuleIO", "write_current");
    ModuleBase::timer::tick("ModuleIO", "write_current");
    //Init_DS_tmp
    Init_DS_tmp(*pv, UHM);
    // construct a DensityMatrix object
    elecstate::DensityMatrix<std::complex<double>, double> DM(&kv,pv,GlobalV::NSPIN);
    
    elecstate::cal_dm_psi(DM.get_paraV_pointer(), pelec->wg, psi[0], DM);
    
    //loc.cal_dm_R(edm_k, ra, edm2d, kv);

    // cal_dm_2d
    DM.init_DMR(ra,&GlobalC::ucell);
    for (int ik = 0; ik < DM.get_DMK_nks(); ++ik)
    {
        cal_tmp_DM(DM, ik, pv->nspin);
        //check later
        double current_ik[3] = {0.0, 0.0, 0.0};
        int total_irr = 0;
#ifdef _OPENMP
#pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            //ModuleBase::matrix local_soverlap(3, 3);
            int local_total_irr = 0;
#else
        //ModuleBase::matrix& local_soverlap = soverlap;
        int& local_total_irr = total_irr;
#endif

            ModuleBase::Vector3<double> tau1, dtau, tau2;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
            {
                const int T1 = GlobalC::ucell.iat2it[iat];
                Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int I1 = GlobalC::ucell.iat2ia[iat];
                // get iat1
                int iat1 = GlobalC::ucell.itia2iat(T1, I1);

                int irr = pv->nlocstart[iat];
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                for (int cb = 0; cb < ra.na_each[iat]; ++cb)
                {
                    const int T2 = ra.info[iat][cb][3];
                    const int I2 = ra.info[iat][cb][4];

                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

                    Atom* atom2 = &GlobalC::ucell.atoms[T2];

                    // get iat2
                    int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                    double Rx = ra.info[iat][cb][0];
                    double Ry = ra.info[iat][cb][1];
                    double Rz = ra.info[iat][cb][2];
                    // get BaseMatrix
                    hamilt::BaseMatrix<double>* tmp_matrix = DM.get_DMR_pointer(1)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                    if(tmp_matrix == nullptr) continue;
                    int row_ap = pv->atom_begin_row[iat1];
                    int col_ap = pv->atom_begin_col[iat2];
                    // get DMR
                    for (int mu = 0; mu < pv->get_row_size(iat1); ++mu)
                    {
                        for (int nu = 0; nu < pv->get_col_size(iat2); ++nu)
                        {
                            // here do not sum over spin due to EDM.sum_DMR_spin();
                            double edm2d1 = tmp_matrix->get_value(mu,nu);
                            double edm2d2 = 2.0 * edm2d1;
                            current_ik[0] -= edm2d2 * UHM.LM->DSloc_Rx[irr];
                            current_ik[1] -= edm2d2 * UHM.LM->DSloc_Ry[irr];
                            current_ik[2] -= edm2d2 * UHM.LM->DSloc_Rz[irr];
                            ++local_total_irr;
                            ++irr;
                        } // end kk
                    }     // end jj
                }         // end cb
            } // end iat
#ifdef _OPENMP
#pragma omp critical(cal_foverlap_k_reduce)
            {
                total_irr += local_total_irr;
            }
        }
#endif
        Parallel_Reduce::reduce_all(current_ik, 3);
        //MPI_Reduce(local_current_ik, current_ik, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(GlobalV::MY_RANK==0)
        {
            std::string filename = GlobalV::global_out_dir + "current_" + std::to_string(ik)+".dat";
            std::ofstream fout;
            fout.open(filename, std::ios::app);
            fout << std::setprecision(16);
            fout << std::scientific;
            //fout << "# current_x current_y current_z" << std::endl;
            fout << current_ik[0] << " " << current_ik[1] << " " << current_ik[2] << std::endl;
            fout.close();
        }
        //write end
        ModuleBase::timer::tick("ModuleIO", "write_current");
    }//end nks
    destory_DS_tmp(UHM);
    return;
}
#endif //__LCAO
