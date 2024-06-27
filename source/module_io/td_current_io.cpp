#include "td_current_io.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_base/vector3.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/module_tddft/td_current.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"

#ifdef __LCAO

void ModuleIO::cal_tmp_DM(elecstate::DensityMatrix<std::complex<double>, double>& DM_real,
                          elecstate::DensityMatrix<std::complex<double>, double>& DM_imag,
                          const int ik,
                          const int nspin,
                          const int is)
{
    ModuleBase::TITLE("ModuleIO", "cal_tmp_DM");
    ModuleBase::timer::tick("ModuleIO", "cal_tmp_DM");
    int ld_hk = DM_real.get_paraV_pointer()->nrow;
    int ld_hk2 = 2 * ld_hk;
    //tmp for is
    int ik_begin = DM_real.get_DMK_nks() / nspin * (is - 1); // jump this->_nks for spin_down if nspin==2

    hamilt::HContainer<double>* tmp_DMR_real = DM_real.get_DMR_vector()[is - 1];
    hamilt::HContainer<double>* tmp_DMR_imag = DM_imag.get_DMR_vector()[is - 1];
    tmp_DMR_real->set_zero();
    tmp_DMR_imag->set_zero();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < tmp_DMR_real->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<double>& tmp_ap_real = tmp_DMR_real->get_atom_pair(i);
        hamilt::AtomPair<double>& tmp_ap_imag = tmp_DMR_imag->get_atom_pair(i);
        int iat1 = tmp_ap_real.get_atom_i();
        int iat2 = tmp_ap_real.get_atom_j();
        // get global indexes of whole matrix for each atom in this process
        int row_ap = DM_real.get_paraV_pointer()->atom_begin_row[iat1];
        int col_ap = DM_real.get_paraV_pointer()->atom_begin_col[iat2];
        if (row_ap == -1 || col_ap == -1)
        {
            throw std::string("Atom-pair not belong this process");
        }
        for (int ir = 0; ir < tmp_ap_real.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> r_index = tmp_ap_real.get_R_index(ir);
            hamilt::BaseMatrix<double>* tmp_matrix_real = tmp_ap_real.find_matrix(r_index);
            hamilt::BaseMatrix<double>* tmp_matrix_imag = tmp_ap_imag.find_matrix(r_index);
#ifdef __DEBUG
            if (tmp_matrix_real == nullptr)
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
                const ModuleBase::Vector3<double> dR(r_index.x, r_index.y, r_index.z);
                const double arg = (DM_real.get_kv_pointer()->kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                double sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                std::complex<double> kphase = std::complex<double>(cosp, sinp);
                // set DMR element
                double* tmp_DMR_real_pointer = tmp_matrix_real->get_pointer();
                double* tmp_DMR_imag_pointer = tmp_matrix_imag->get_pointer();
                std::complex<double>* tmp_DMK_pointer = DM_real.get_DMK_pointer(ik + ik_begin);
                double* DMK_real_pointer = nullptr;
                double* DMK_imag_pointer = nullptr;
                // jump DMK to fill DMR
                // DMR is row-major, DMK is column-major
                tmp_DMK_pointer += col_ap * DM_real.get_paraV_pointer()->nrow + row_ap;
                for (int mu = 0; mu < DM_real.get_paraV_pointer()->get_row_size(iat1); ++mu)
                {
                    DMK_real_pointer = (double*)tmp_DMK_pointer;
                    DMK_imag_pointer = DMK_real_pointer + 1;
                    //calculate real part
                    BlasConnector::axpy(DM_real.get_paraV_pointer()->get_col_size(iat2),
                                        -kphase.imag(),
                                        DMK_imag_pointer,
                                        ld_hk2,
                                        tmp_DMR_real_pointer,
                                        1);
                    BlasConnector::axpy(DM_real.get_paraV_pointer()->get_col_size(iat2),
                                        kphase.real(),
                                        DMK_real_pointer,
                                        ld_hk2,
                                        tmp_DMR_real_pointer,
                                        1);
                    //calculate imag part
                    BlasConnector::axpy(DM_imag.get_paraV_pointer()->get_col_size(iat2),
                                        kphase.imag(),
                                        DMK_real_pointer,
                                        ld_hk2,
                                        tmp_DMR_imag_pointer,
                                        1);
                    BlasConnector::axpy(DM_imag.get_paraV_pointer()->get_col_size(iat2),
                                        kphase.real(),
                                        DMK_imag_pointer,
                                        ld_hk2,
                                        tmp_DMR_imag_pointer,
                                        1);
                    tmp_DMK_pointer += 1;
                    tmp_DMR_real_pointer += DM_real.get_paraV_pointer()->get_col_size(iat2);
                    tmp_DMR_imag_pointer += DM_imag.get_paraV_pointer()->get_col_size(iat2);
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
                             const TwoCenterIntegrator* intor,
                             const Parallel_Orbitals* pv,
                             Record_adj& ra,
                             hamilt::HContainer<double>* sR) 
{

    ModuleBase::TITLE("ModuleIO", "write_current");
    ModuleBase::timer::tick("ModuleIO", "write_current");
    TD_current* current_term =nullptr;
    ModuleBase::Vector3<double> cart_At(0,0,0);
    if (TD_Velocity::td_vel_op!=nullptr)
    {
        current_term = TD_Velocity::td_vel_op->get_current_pointer();
        cart_At=TD_Velocity::td_vel_op->cart_At;
    }
    else
    {
        current_term = new TD_current(&GlobalC::ucell, &GlobalC::GridD, pv, cart_At, intor);
        current_term->calculate_grad_term();
        current_term->calculate_vcomm_r();
    }
    // construct a DensityMatrix object
    elecstate::DensityMatrix<std::complex<double>, double> DM_real(&kv, pv, GlobalV::NSPIN);
    elecstate::DensityMatrix<std::complex<double>, double> DM_imag(&kv, pv, GlobalV::NSPIN);
    //calculate DMK
    elecstate::cal_dm_psi(DM_real.get_paraV_pointer(), pelec->wg, psi[0], DM_real);

    // init DMR
    DM_real.init_DMR(ra,&GlobalC::ucell);
    DM_imag.init_DMR(ra,&GlobalC::ucell);

    int nks = DM_real.get_DMK_nks();
    if(GlobalV::NSPIN==2)
    {
        nks/=2;
    }
    double current_total[3]={0.0,0.0,0.0};
    for (int is = 1; is <= GlobalV::NSPIN; ++is)
    {
        for (int ik = 0; ik < nks; ++ik)
        {
            cal_tmp_DM(DM_real, DM_imag, ik, GlobalV::NSPIN, is);
            //check later
            double current_ik[3]={0.0,0.0,0.0};
            int total_irr = 0;
#ifdef _OPENMP
#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                double local_current_ik[3] = {0.0, 0.0, 0.0};
                int local_total_irr = 0;
#else
            //ModuleBase::matrix& local_soverlap = soverlap;
            double* local_current_ik = current_ik;
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
                        //std::cout<< "iat1: " << iat1 << " iat2: " << iat2 << " Rx: " << Rx << " Ry: " << Ry << " Rz: " << Rz << std::endl;
                        // get BaseMatrix
                        hamilt::BaseMatrix<double>* tmp_matrix_real = DM_real.get_DMR_pointer(is)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        hamilt::BaseMatrix<double>* tmp_matrix_imag = DM_imag.get_DMR_pointer(is)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        //refactor
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvx = current_term->get_current_term_pointer(0)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvy = current_term->get_current_term_pointer(1)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvz = current_term->get_current_term_pointer(2)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        //hamilt::BaseMatrix<std::complex<double>>* tmp_nonlocal = current_term->get_nonlocal_pointer()->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        //get SR info
                        hamilt::BaseMatrix<double>* tmp_s = sR->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        if(tmp_matrix_real == nullptr) continue;
                        int row_ap = pv->atom_begin_row[iat1];
                        int col_ap = pv->atom_begin_col[iat2];
                        // get DMR
                        for (int mu = 0; mu < pv->get_row_size(iat1); ++mu)
                        {
                            for (int nu = 0; nu < pv->get_col_size(iat2); ++nu)
                            {
                                double dm2d1_real = tmp_matrix_real->get_value(mu,nu);
                                double dm2d1_imag = tmp_matrix_imag->get_value(mu,nu);

                                double s_tmp = 0.0;
                                std::complex<double> rvx = {0,0};
                                std::complex<double> rvy = {0,0};
                                std::complex<double> rvz = {0,0};

                                if(tmp_m_rvx!=nullptr)
                                {
                                    rvx = tmp_m_rvx->get_value(mu,nu);
                                    rvy = tmp_m_rvy->get_value(mu,nu);
                                    rvz = tmp_m_rvz->get_value(mu,nu);
                                }
                                if(tmp_s!=nullptr)
                                {
                                    s_tmp = tmp_s->get_value(mu,nu);
                                }

                                local_current_ik[0] -= dm2d1_real * rvx.real() - dm2d1_imag * rvx.imag() + dm2d1_real*cart_At[0]*s_tmp/2.0;
                                local_current_ik[1] -= dm2d1_real * rvy.real() - dm2d1_imag * rvy.imag() + dm2d1_real*cart_At[1]*s_tmp/2.0;
                                local_current_ik[2] -= dm2d1_real * rvz.real() - dm2d1_imag * rvz.imag() + dm2d1_real*cart_At[2]*s_tmp/2.0;
                                ++local_total_irr;
                                ++irr;
                            } // end kk
                        }     // end jj
                    }         // end cb
                } // end iat
#ifdef _OPENMP
#pragma omp critical(cal_current_k_reduce)
                {
                    total_irr += local_total_irr;
                    for (int i = 0; i < 3; ++i)
                    {
                        current_ik[i] += local_current_ik[i];
                    }
                }
            }
#endif
            Parallel_Reduce::reduce_all(current_ik, 3);
            for(int i=0;i<3;++i)
            {
                current_total[i] += current_ik[i];
            }
            //MPI_Reduce(local_current_ik, current_ik, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(GlobalV::MY_RANK==0 && TD_Velocity::out_current_k)
            {
                std::string filename = GlobalV::global_out_dir + "current_spin"+ std::to_string(is) + "_ik" + std::to_string(ik)+".dat";
                std::ofstream fout;
                fout.open(filename, std::ios::app);
                fout << std::setprecision(16);
                fout << std::scientific;
                fout << istep << " " << current_ik[0] << " " << current_ik[1] << " " << current_ik[2] << std::endl;
                fout.close();
            }
            //write end
            ModuleBase::timer::tick("ModuleIO", "write_current");
        }//end nks
    }//end is
    if(GlobalV::MY_RANK==0)
    {
        std::string filename = GlobalV::global_out_dir + "current_total.dat";
        std::ofstream fout;
        fout.open(filename, std::ios::app);
        fout << std::setprecision(16);
        fout << std::scientific;
        fout << istep << " " << current_total[0] << " " << current_total[1] << " " << current_total[2] << std::endl;
        fout.close();
    }
    delete current_term;
    return;
}
#endif //__LCAO
