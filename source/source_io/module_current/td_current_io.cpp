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
#include "source_io/module_parameter/parameter.h"

#ifdef __LCAO
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
                             TD_info* td_p,
                             Record_adj& ra)
{

    ModuleBase::TITLE("ModuleIO", "write_current");
    ModuleBase::timer::start("ModuleIO", "write_current");
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
    elecstate::DensityMatrix<std::complex<double>, std::complex<double>> tmp_dm(pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);
    // calculate DMK
    elecstate::cal_dm_psi(pv, pelec->wg, psi[0], tmp_dm);

    // init DMR
    tmp_dm.init_DMR(ra, &ucell);

    if(PARAM.inp.td_stype!=2)
    {
        tmp_dm.cal_DMR();
    }
    else
    {
        tmp_dm.cal_DMR_td(td_p->get_phase_hybrid(),TD_info::cart_At);
    }
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
                hamilt::BaseMatrix<std::complex<double>>* tmp_matrix
                    = tmp_dm.get_DMR_pointer(1)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                // refactor
                hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvx
                    = current_term[0]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvy
                    = current_term[1]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvz
                    = current_term[2]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                if (tmp_matrix == nullptr)
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
                        std::complex<double> dm2d1 = tmp_matrix->get_value(mu, nu);

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
                        local_current[0] -= dm2d1.real() * rvx.real() - dm2d1.imag() * rvx.imag();    
                        local_current[1] -= dm2d1.real() * rvy.real() - dm2d1.imag() * rvy.imag();
                        local_current[2] -= dm2d1.real() * rvz.real() - dm2d1.imag() * rvz.imag();
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
        std::string filename = PARAM.globalv.global_out_dir + "current_tot.txt";
        std::ofstream fout;
        fout.open(filename, std::ios::app);
        fout << std::setprecision(16);
        fout << std::scientific;
        fout << istep+1 << " " << current_total[0]/omega 
             << " " << current_total[1]/omega 
             << " " << current_total[2]/omega << std::endl;
        fout.close();
    }

    ModuleBase::timer::end("ModuleIO", "write_current");
    return;
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
                             TD_info* td_p,
                             Record_adj& ra)
{

    ModuleBase::TITLE("ModuleIO", "write_current");
    ModuleBase::timer::start("ModuleIO", "write_current");
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
    elecstate::DensityMatrix<std::complex<double>, std::complex<double>> tmp_dm(pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);
    //elecstate::DensityMatrix<std::complex<double>, double> DM_real(pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);
    //elecstate::DensityMatrix<std::complex<double>, double> DM_imag(pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);
    // calculate DMK
    elecstate::cal_dm_psi(pv, pelec->wg, psi[0], tmp_dm);

    // init DMR
    tmp_dm.init_DMR(ra, &ucell);

    int nks = tmp_dm.get_DMK_nks() / nspin_dm;
    double current_total[3] = {0.0, 0.0, 0.0};
    for (int is = 1; is <= nspin_dm; ++is)
    {
        for (int ik = 0; ik < nks; ++ik)
        {
            if(PARAM.inp.td_stype!=2)
            {
                tmp_dm.cal_DMR(ik);
            }
            else
            {
                tmp_dm.cal_DMR_td(td_p->get_phase_hybrid(),TD_info::cart_At,ik);
            }
            
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
                        hamilt::BaseMatrix<std::complex<double>>* tmp_matrix
                            = tmp_dm.get_DMR_pointer(is)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        // refactor
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvx
                            = current_term[0]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvy
                            = current_term[1]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        hamilt::BaseMatrix<std::complex<double>>* tmp_m_rvz
                            = current_term[2]->find_matrix(iat1, iat2, Rx, Ry, Rz);
                        if (tmp_matrix == nullptr)
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
                                std::complex<double> dm2d1 = tmp_matrix->get_value(mu, nu);

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
                                local_current_ik[0] -= dm2d1.real() * rvx.real() - dm2d1.imag() * rvx.imag();    
                                local_current_ik[1] -= dm2d1.real() * rvy.real() - dm2d1.imag() * rvy.imag();
                                local_current_ik[2] -= dm2d1.real() * rvz.real() - dm2d1.imag() * rvz.imag();
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
                std::string filename = PARAM.globalv.global_out_dir + "current_s" + std::to_string(is) + "k"
                                       + std::to_string(ik+1) + ".txt";
                std::ofstream fout;
                fout.open(filename, std::ios::app);
                fout << std::setprecision(16);
                fout << std::scientific;
                fout << istep+1 << " " << current_ik[0]/omega 
                     << " " << current_ik[1]/omega 
                     << " " << current_ik[2]/omega << std::endl;
                fout.close();
            }
            // write end
        } // end nks
    } // end is
    if (GlobalV::MY_RANK == 0)
    {
        std::string filename = PARAM.globalv.global_out_dir + "current_tot.txt";
        std::ofstream fout;
        fout.open(filename, std::ios::app);
        fout << std::setprecision(16);
        fout << std::scientific;
        fout << istep+1 << " " << current_total[0]/omega 
                        << " " << current_total[1]/omega 
                        << " " << current_total[2]/omega << std::endl;
        fout.close();
    }

    ModuleBase::timer::end("ModuleIO", "write_current");
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
                        TD_info* td_p,
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
                        TD_info* td_p,
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
                TD_info* td_p,
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
                TD_info* td_p,
                Record_adj& ra);
#endif //__LCAO

