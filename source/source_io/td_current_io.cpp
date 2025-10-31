#include "td_current_io.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/libm/libm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_threading.h"
#include "source_base/vector3.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_lcao/LCAO_domain.h"
#include "source_lcao/module_rt/td_info.h"
#include "source_pw/module_pwdft/global.h"
#include "source_io/module_parameter/parameter.h"

#ifdef __LCAO
void ModuleIO::cal_tmp_DM(const UnitCell& ucell,
                        elecstate::DensityMatrix<std::complex<double>, double>& DM_real,
                        elecstate::DensityMatrix<std::complex<double>, double>& DM_imag,
                        int nspin_dm)
{
    ModuleBase::TITLE("ModuleIO", "cal_tmp_DM");
    ModuleBase::timer::tick("ModuleIO", "cal_tmp_DM");
    for (int is = 1; is <= nspin_dm; ++is)
    {
        for (int ik = 0; ik < DM_real.get_DMK_nks() / nspin_dm; ++ik)
        {
            cal_tmp_DM_k(ucell, DM_real, DM_imag, ik, nspin_dm, is, false);
        }
    }
    ModuleBase::timer::tick("ModuleIO", "cal_tmp_DM");
}
template <typename TR>
void ModuleIO::write_current(const UnitCell& ucell,
                             const int istep,
                             const psi::Psi<std::complex<double>>* psi,
                             const elecstate::ElecState* pelec,
                             const K_Vectors& kv,
                             const TwoCenterIntegrator* intor,
                             const Parallel_Orbitals* pv,
                             const LCAO_Orbitals& orb,
                             const Velocity_op<TR>* cal_current,
                             Record_adj& ra)
{

    ModuleBase::TITLE("ModuleIO", "write_current");
    ModuleBase::timer::tick("ModuleIO", "write_current");
    std::vector<hamilt::HContainer<std::complex<double>>*> current_term = {nullptr, nullptr, nullptr};
    if (PARAM.inp.td_stype!=1)
    {
        for (int dir = 0; dir < 3; dir++)
        {
            current_term[dir] = cal_current->get_current_term_pointer(dir);
        }
    }
    else
    {
        if (TD_info::td_vel_op == nullptr)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_current", "velocity gauge infos is null!");
        }
        for (int dir = 0; dir < 3; dir++)
        {
            current_term[dir] = TD_info::td_vel_op->get_current_term_pointer(dir);
        }
    }
    double omega=ucell.omega;
    // construct a DensityMatrix object
    // Since the function cal_dm_psi do not suport DMR in complex type, I replace it with two DMR in double type. Should
    // be refactored in the future.
    const int nspin0 = PARAM.inp.nspin;
    const int nspin_dm = std::map<int, int>({ {1,1},{2,2},{4,1} })[nspin0];
    elecstate::DensityMatrix<std::complex<double>, double> DM_real(pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);
    elecstate::DensityMatrix<std::complex<double>, double> DM_imag(pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);
    // calculate DMK
    elecstate::cal_dm_psi(DM_real.get_paraV_pointer(), pelec->wg, psi[0], DM_real);

    // init DMR
    DM_real.init_DMR(ra, &ucell);
    DM_imag.init_DMR(ra, &ucell);
    cal_tmp_DM(ucell, DM_real, DM_imag, nspin_dm);
    //DM_real.sum_DMR_spin();
    //DM_imag.sum_DMR_spin();

    double current_total[3] = {0.0, 0.0, 0.0};
#ifdef _OPENMP
#pragma omp parallel
    {
        double local_current[3] = {0.0, 0.0, 0.0};
#else
        // ModuleBase::matrix& local_soverlap = soverlap;
        double* local_current = current_total;
#endif
        ModuleBase::Vector3<double> tau1, dtau, tau2;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            const int T1 = ucell.iat2it[iat];
            Atom* atom1 = &ucell.atoms[T1];
            const int I1 = ucell.iat2ia[iat];
            // get iat1
            int iat1 = ucell.itia2iat(T1, I1);
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            for (int cb = 0; cb < ra.na_each[iat]; ++cb)
            {
                const int T2 = ra.info[iat][cb][3];
                const int I2 = ra.info[iat][cb][4];

                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

                Atom* atom2 = &ucell.atoms[T2];

                // get iat2
                int iat2 = ucell.itia2iat(T2, I2);
                double Rx = ra.info[iat][cb][0];
                double Ry = ra.info[iat][cb][1];
                double Rz = ra.info[iat][cb][2];
                //std::cout<< "iat1: " << iat1 << " iat2: " << iat2 << " Rx: " << Rx << " Ry: " << Ry << " Rz:" << Rz << std::endl;
                //  get BaseMatrix
                hamilt::BaseMatrix<double>* tmp_matrix_real
                    = DM_real.get_DMR_pointer(1)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                hamilt::BaseMatrix<double>* tmp_matrix_imag
                    = DM_imag.get_DMR_pointer(1)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                // refactor
                hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvx
                    = current_term[0]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvy
                    = current_term[1]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvz
                    = current_term[2]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                if (tmp_matrix_real == nullptr)
                {
                    continue;
                }
                int row_ap = pv->atom_begin_row[iat1];
                int col_ap = pv->atom_begin_col[iat2];
                // get DMR
                for (int mu = 0; mu < pv->get_row_size(iat1); ++mu)
                {
                    for (int nu = 0; nu < pv->get_col_size(iat2); ++nu)
                    {
                        double dm2d1_real = tmp_matrix_real->get_value(mu, nu);
                        double dm2d1_imag = tmp_matrix_imag->get_value(mu, nu);

                        std::complex<double> rvx = {0, 0};
                        std::complex<double> rvy = {0, 0};
                        std::complex<double> rvz = {0, 0};

                        if (tmp_m_rvx != nullptr)
                        {
                            rvx = tmp_m_rvx->get_value(mu, nu);
                            rvy = tmp_m_rvy->get_value(mu, nu);
                            rvz = tmp_m_rvz->get_value(mu, nu);
                        }
                        //std::cout<<"mu: "<< mu <<" nu: "<< nu << std::endl;
                        // std::cout<<"dm2d1_real: "<< dm2d1_real << " dm2d1_imag: "<< dm2d1_imag << std::endl;
                        //std::cout<<"rvz: "<< rvz.real() << " " << rvz.imag() << std::endl;
                        local_current[0] -= dm2d1_real * rvx.real() - dm2d1_imag * rvx.imag();    
                        local_current[1] -= dm2d1_real * rvy.real() - dm2d1_imag * rvy.imag();
                        local_current[2] -= dm2d1_real * rvz.real() - dm2d1_imag * rvz.imag();
                    } // end kk
                } // end jj
            } // end cb
        } // end iat
#ifdef _OPENMP
#pragma omp critical(cal_current_k_reduce)
        {
            for (int i = 0; i < 3; ++i)
            {
                current_total[i] += local_current[i];
            }
        }
    }
#endif
    Parallel_Reduce::reduce_all(current_total, 3);
    // write end
    if (GlobalV::MY_RANK == 0)
    {
        std::string filename = PARAM.globalv.global_out_dir + "current_total.txt";
        std::ofstream fout;
        fout.open(filename, std::ios::app);
        fout << std::setprecision(16);
        fout << std::scientific;
        fout << istep << " " << current_total[0]/omega << " " << current_total[1]/omega << " " << current_total[2]/omega << std::endl;
        fout.close();
    }

    ModuleBase::timer::tick("ModuleIO", "write_current");
    return;
}
void ModuleIO::cal_tmp_DM_k(const UnitCell& ucell,
                          elecstate::DensityMatrix<std::complex<double>, double>& DM_real,
                          elecstate::DensityMatrix<std::complex<double>, double>& DM_imag,
                          const int ik,
                          const int nspin,
                          const int is,
                          const bool reset)
{
    ModuleBase::TITLE("ModuleIO", "cal_tmp_DM_k");
    ModuleBase::timer::tick("ModuleIO", "cal_tmp_DM_k");
    int ld_hk = DM_real.get_paraV_pointer()->nrow;
    int ld_hk2 = 2 * ld_hk;
    // tmp for is
    int ik_begin = DM_real.get_DMK_nks() / nspin * (is - 1); // jump nk for spin_down if nspin==2
    //sum spin up and down into up
    hamilt::HContainer<double>* tmp_DMR_real = DM_real.get_DMR_vector()[0];
    hamilt::HContainer<double>* tmp_DMR_imag = DM_imag.get_DMR_vector()[0];
    if(reset)
    {
        tmp_DMR_real->set_zero();
        tmp_DMR_imag->set_zero();
    }
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
        // SOC
        std::vector<std::complex<double>> tmp_DMR;
        if (PARAM.inp.nspin == 4)
        {
            tmp_DMR.resize(tmp_ap_real.get_size());
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
            if (PARAM.inp.nspin != 4)
            {
                double arg_td = 0.0;
                if(elecstate::H_TDDFT_pw::stype == 2)
                {
                    //cal tddft phase for hybrid gauge
                    const int iat1 = tmp_ap_real.get_atom_i();
                    const int iat2 = tmp_ap_real.get_atom_j();
                    ModuleBase::Vector3<double> dtau = ucell.cal_dtau(iat1, iat2, r_index);
                    double& tmp_lat0 = ucell.lat0;
                    arg_td = TD_info::cart_At * dtau * tmp_lat0;
                }
                // cal k_phase
                // if TK==std::complex<double>, kphase is e^{ikR}
                const ModuleBase::Vector3<double> dR(r_index.x, r_index.y, r_index.z);
                const double arg = (DM_real.get_kvec_d()[ik] * dR) * ModuleBase::TWO_PI + arg_td;
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
                    // calculate real part
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
                    // calculate imag part
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
            // treat DMR as pauli matrix when NSPIN=4
            if (PARAM.inp.nspin == 4)
            {
                tmp_DMR.assign(tmp_ap_real.get_size(), std::complex<double>(0.0, 0.0));
                {
                    // cal k_phase
                    // if TK==std::complex<double>, kphase is e^{ikR}
                    const ModuleBase::Vector3<double> dR(r_index.x, r_index.y, r_index.z);
                    double arg_td = 0.0;
                    if(elecstate::H_TDDFT_pw::stype == 2)
                    {
                        //new
                        //cal tddft phase for mixing gauge
                        const int iat1 = tmp_ap_real.get_atom_i();
                        const int iat2 = tmp_ap_real.get_atom_j();
                        ModuleBase::Vector3<double> dtau = ucell.cal_dtau(iat1, iat2, r_index);
                        double& tmp_lat0 = ucell.lat0;
                        arg_td = TD_info::cart_At * dtau * tmp_lat0;
                    }
                    const double arg = (DM_real.get_kvec_d()[ik] * dR) * ModuleBase::TWO_PI + arg_td;
                    double sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    std::complex<double> kphase = std::complex<double>(cosp, sinp);
                    // set DMR element
                    std::complex<double>* tmp_DMR_pointer = tmp_DMR.data();
                    std::complex<double>* tmp_DMK_pointer = DM_real.get_DMK_pointer(ik + ik_begin);;
                    double* DMK_real_pointer = nullptr;
                    double* DMK_imag_pointer = nullptr;
                    // jump DMK to fill DMR
                    // DMR is row-major, DMK is column-major
                    tmp_DMK_pointer += col_ap * DM_real.get_paraV_pointer()->nrow + row_ap;
                    for (int mu = 0; mu < tmp_ap_real.get_row_size(); ++mu)
                    {
                        BlasConnector::axpy(tmp_ap_real.get_col_size(),
                                            kphase,
                                            tmp_DMK_pointer,
                                            ld_hk,
                                            tmp_DMR_pointer,
                                            1);
                        tmp_DMK_pointer += 1;
                        tmp_DMR_pointer += tmp_ap_real.get_col_size();
                    }
                }
                int npol = 2;
                // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
                int step_trace[4];
                for (int is = 0; is < npol; is++)
                {
                    for (int is2 = 0; is2 < npol; is2++)
                    {
                        step_trace[is * npol + is2] = tmp_ap_real.get_col_size() * is + is2;
                    }
                }
                std::complex<double> tmp[4];
                double* target_DMR_real = tmp_matrix_real->get_pointer();
                double* target_DMR_imag = tmp_matrix_imag->get_pointer();
                std::complex<double>* tmp_DMR_pointer = tmp_DMR.data();
                for (int irow = 0; irow < tmp_ap_real.get_row_size(); irow += 2)
                {
                    for (int icol = 0; icol < tmp_ap_real.get_col_size(); icol += 2)
                    {
                        // catch the 4 spin component value of one orbital pair
                        tmp[0] = tmp_DMR_pointer[icol + step_trace[0]];
                        tmp[1] = tmp_DMR_pointer[icol + step_trace[1]];
                        tmp[2] = tmp_DMR_pointer[icol + step_trace[2]];
                        tmp[3] = tmp_DMR_pointer[icol + step_trace[3]];
                        // transfer to Pauli matrix and save the real part
                        // save them back to the tmp_matrix
                        target_DMR_real[icol + step_trace[0]] += tmp[0].real() + tmp[3].real();
                        target_DMR_real[icol + step_trace[1]] += tmp[1].real() + tmp[2].real();
                        target_DMR_real[icol + step_trace[2]]
                            += -tmp[1].imag() + tmp[2].imag(); // (i * (rho_updown - rho_downup)).real()
                        target_DMR_real[icol + step_trace[3]] += tmp[0].real() - tmp[3].real();
                        //imag part
                        target_DMR_imag[icol + step_trace[0]] += tmp[0].imag() + tmp[3].imag();
                        target_DMR_imag[icol + step_trace[1]] += tmp[1].imag() + tmp[2].imag();
                        target_DMR_imag[icol + step_trace[2]]
                            += tmp[1].real() - tmp[2].real(); // (i * (rho_updown - rho_downup)).real()
                        target_DMR_imag[icol + step_trace[3]] += tmp[0].imag() - tmp[3].imag();
                    }
                    tmp_DMR_pointer += tmp_ap_real.get_col_size() * 2;
                    target_DMR_real += tmp_ap_real.get_col_size() * 2;
                    target_DMR_imag += tmp_ap_real.get_col_size() * 2;
                }
            }
        }
    }
    ModuleBase::timer::tick("ModuleIO", "cal_tmp_DM_k");
}
template <typename TR>
void ModuleIO::write_current_eachk(const UnitCell& ucell,
                             const int istep,
                             const psi::Psi<std::complex<double>>* psi,
                             const elecstate::ElecState* pelec,
                             const K_Vectors& kv,
                             const TwoCenterIntegrator* intor,
                             const Parallel_Orbitals* pv,
                             const LCAO_Orbitals& orb,
                             const Velocity_op<TR>* cal_current,
                             Record_adj& ra)
{

    ModuleBase::TITLE("ModuleIO", "write_current");
    ModuleBase::timer::tick("ModuleIO", "write_current");
    std::vector<hamilt::HContainer<std::complex<double>>*> current_term = {nullptr, nullptr, nullptr};
    if (PARAM.inp.td_stype != 1)
    {
        for (int dir = 0; dir < 3; dir++)
        {
            current_term[dir] = cal_current->get_current_term_pointer(dir);
        }
    }
    else
    {
        if (TD_info::td_vel_op == nullptr)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_current", "velocity gauge infos is null!");
        }
        for (int dir = 0; dir < 3; dir++)
        {
            current_term[dir] = TD_info::td_vel_op->get_current_term_pointer(dir);
        }
    }
    double omega=ucell.omega;
    // construct a DensityMatrix object
    // Since the function cal_dm_psi do not suport DMR in complex type, 
    // I replace it with two DMR in double type.
    // Should be refactored in the future.

    const int nspin0 = PARAM.inp.nspin;
    const int nspin_dm = std::map<int, int>({ {1,1},{2,2},{4,1} })[nspin0];
    elecstate::DensityMatrix<std::complex<double>, double> DM_real(pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);
    elecstate::DensityMatrix<std::complex<double>, double> DM_imag(pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);
    // calculate DMK
    elecstate::cal_dm_psi(DM_real.get_paraV_pointer(), pelec->wg, psi[0], DM_real);

    // init DMR
    DM_real.init_DMR(ra, &ucell);
    DM_imag.init_DMR(ra, &ucell);

    int nks = DM_real.get_DMK_nks() / nspin_dm;
    double current_total[3] = {0.0, 0.0, 0.0};
    for (int is = 1; is <= nspin_dm; ++is)
    {
        for (int ik = 0; ik < nks; ++ik)
        {
            cal_tmp_DM_k(ucell, DM_real, DM_imag, ik, nspin_dm, is);
            // check later
            double current_ik[3] = {0.0, 0.0, 0.0};
#ifdef _OPENMP
#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                double local_current_ik[3] = {0.0, 0.0, 0.0};
#else
            // ModuleBase::matrix& local_soverlap = soverlap;
            double* local_current_ik = current_ik;
#endif

                ModuleBase::Vector3<double> tau1, dtau, tau2;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
                for (int iat = 0; iat < ucell.nat; iat++)
                {
                    const int T1 = ucell.iat2it[iat];
                    Atom* atom1 = &ucell.atoms[T1];
                    const int I1 = ucell.iat2ia[iat];
                    // get iat1
                    int iat1 = ucell.itia2iat(T1, I1);
                    const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                    for (int cb = 0; cb < ra.na_each[iat]; ++cb)
                    {
                        const int T2 = ra.info[iat][cb][3];
                        const int I2 = ra.info[iat][cb][4];

                        const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

                        Atom* atom2 = &ucell.atoms[T2];

                        // get iat2
                        int iat2 = ucell.itia2iat(T2, I2);
                        double Rx = ra.info[iat][cb][0];
                        double Ry = ra.info[iat][cb][1];
                        double Rz = ra.info[iat][cb][2];
                        //std::cout<< "iat1: " << iat1 << " iat2: " << iat2 << " Rx: " << Rx << " Ry: " << Ry << " Rz:" << Rz << std::endl;
                        //  get BaseMatrix
                        hamilt::BaseMatrix<double>* tmp_matrix_real
                            = DM_real.get_DMR_pointer(is)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        hamilt::BaseMatrix<double>* tmp_matrix_imag
                            = DM_imag.get_DMR_pointer(is)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        // refactor
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvx
                            = current_term[0]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvy
                            = current_term[1]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvz
                            = current_term[2]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        if (tmp_matrix_real == nullptr)
                        {
                            continue;
                        }
                        int row_ap = pv->atom_begin_row[iat1];
                        int col_ap = pv->atom_begin_col[iat2];
                        // get DMR
                        for (int mu = 0; mu < pv->get_row_size(iat1); ++mu)
                        {
                            for (int nu = 0; nu < pv->get_col_size(iat2); ++nu)
                            {
                                double dm2d1_real = tmp_matrix_real->get_value(mu, nu);
                                double dm2d1_imag = tmp_matrix_imag->get_value(mu, nu);

                                std::complex<double> rvx = {0, 0};
                                std::complex<double> rvy = {0, 0};
                                std::complex<double> rvz = {0, 0};

                                if (tmp_m_rvx != nullptr)
                                {
                                    rvx = tmp_m_rvx->get_value(mu, nu);
                                    rvy = tmp_m_rvy->get_value(mu, nu);
                                    rvz = tmp_m_rvz->get_value(mu, nu);
                                }
                                // std::cout<<"mu: "<< mu <<" nu: "<< nu << std::endl;
                                // std::cout<<"dm2d1_real: "<< dm2d1_real << " dm2d1_imag: "<< dm2d1_imag << std::endl;
                                // std::cout<<"rvz: "<< rvz.real() << " " << rvz.imag() << std::endl;
                                local_current_ik[0] -= dm2d1_real * rvx.real() - dm2d1_imag * rvx.imag();    
                                local_current_ik[1] -= dm2d1_real * rvy.real() - dm2d1_imag * rvy.imag();
                                local_current_ik[2] -= dm2d1_real * rvz.real() - dm2d1_imag * rvz.imag();
                            } // end kk
                        } // end jj
                    } // end cb
                } // end iat
#ifdef _OPENMP
#pragma omp critical(cal_current_k_reduce)
                {
                    for (int i = 0; i < 3; ++i)
                    {
                        current_ik[i] += local_current_ik[i];
                    }
                }
            }
#endif
            Parallel_Reduce::reduce_all(current_ik, 3);
            for (int i = 0; i < 3; ++i)
            {
                current_total[i] += current_ik[i];
            }
            // MPI_Reduce(local_current_ik, current_ik, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (GlobalV::MY_RANK == 0 && TD_info::out_current_k)
            {
                std::string filename = PARAM.globalv.global_out_dir + "current_spin" + std::to_string(is) + "_ik"
                                       + std::to_string(ik) + ".txt";
                std::ofstream fout;
                fout.open(filename, std::ios::app);
                fout << std::setprecision(16);
                fout << std::scientific;
                fout << istep << " " << current_ik[0]/omega << " " << current_ik[1]/omega << " " << current_ik[2]/omega << std::endl;
                fout.close();
            }
            // write end
        } // end nks
    } // end is
    if (GlobalV::MY_RANK == 0)
    {
        std::string filename = PARAM.globalv.global_out_dir + "current_total.txt";
        std::ofstream fout;
        fout.open(filename, std::ios::app);
        fout << std::setprecision(16);
        fout << std::scientific;
        fout << istep << " " << current_total[0]/omega << " " << current_total[1]/omega << " " << current_total[2]/omega << std::endl;
        fout.close();
    }

    ModuleBase::timer::tick("ModuleIO", "write_current");
    return;
}
template 
void ModuleIO::write_current_eachk<double>(
                        const UnitCell& ucell,
                        const int istep,
                        const psi::Psi<std::complex<double>>* psi,
                        const elecstate::ElecState* pelec,
                        const K_Vectors& kv,
                        const TwoCenterIntegrator* intor,
                        const Parallel_Orbitals* pv,
                        const LCAO_Orbitals& orb,
                        const Velocity_op<double>* cal_current,
                        Record_adj& ra);
template 
void ModuleIO::write_current_eachk<std::complex<double>>(const UnitCell& ucell,
                        const int istep,
                        const psi::Psi<std::complex<double>>* psi,
                        const elecstate::ElecState* pelec,
                        const K_Vectors& kv,
                        const TwoCenterIntegrator* intor,
                        const Parallel_Orbitals* pv,
                        const LCAO_Orbitals& orb,
                        const Velocity_op<std::complex<double>>* cal_current,
                        Record_adj& ra);
template 
void ModuleIO::write_current<double>(const UnitCell& ucell,
                const int istep,
                const psi::Psi<std::complex<double>>* psi,
                const elecstate::ElecState* pelec,
                const K_Vectors& kv,
                const TwoCenterIntegrator* intor,
                const Parallel_Orbitals* pv,
                const LCAO_Orbitals& orb,
                const Velocity_op<double>* cal_current,
                Record_adj& ra);
template 
void ModuleIO::write_current<std::complex<double>>(const UnitCell& ucell,
                const int istep,
                const psi::Psi<std::complex<double>>* psi,
                const elecstate::ElecState* pelec,
                const K_Vectors& kv,
                const TwoCenterIntegrator* intor,
                const Parallel_Orbitals* pv,
                const LCAO_Orbitals& orb,
                const Velocity_op<std::complex<double>>* cal_current,
                Record_adj& ra);
#endif //__LCAO

