#include "spin_constrain.h"
#include "module_base/parallel_common.h"

/// @brief  set input parameters
template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::set_input_parameters(double sc_thr_in,
                                                                                int nsc_in,
                                                                                int nsc_min_in,
                                                                                double alpha_trial_in,
                                                                                double sccut_in,
                                                                                bool decay_grad_switch_in)
{
    this->sc_thr_ = sc_thr_in;
    this->nsc_ = nsc_in;
    this->nsc_min_ = nsc_min_in;
    this->alpha_trial_ = alpha_trial_in / ModuleBase::Ry_to_eV;
    this->restrict_current_ = sccut_in / ModuleBase::Ry_to_eV;
    this->decay_grad_switch_ = decay_grad_switch_in;
}

/// @brief  set ParaV
template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::set_ParaV(Parallel_Orbitals* ParaV_in)
{
    this->ParaV = ParaV_in;
    int nloc = this->ParaV->nloc;
    if (nloc <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::init_sc", "nloc <= 0");
    }
}

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::set_solver_parameters(
    int nspin_in,
    K_Vectors kv_in,
    hsolver::HSolver<std::complex<double>, psi::DEVICE_CPU>* phsol_in,
    hamilt::Hamilt<std::complex<double>, psi::DEVICE_CPU>* p_hamilt_in,
    psi::Psi<std::complex<double>>* psi_in,
    elecstate::ElecState* pelec_in,
    std::string KS_SOLVER_in,
    LCAO_Matrix* LM_in)
{
    /// set nspin
    this->set_nspin(nspin_in);
    this->kv_ = kv_in;
    this->phsol = phsol_in;
    this->p_hamilt = p_hamilt_in;
    this->psi = psi_in;
    this->pelec = pelec_in;
    this->KS_SOLVER = KS_SOLVER_in;
    this->LM = LM_in;
}

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::bcast_ScData(const UnitCell& ucell, int NPOL, std::string sc_file)
{
    /// set ScData
    this->clear_ScData();
    this->clear_atomCounts();
    std::map<int, int> atomCounts = ucell.get_atomCounts();
    std::map<int, int> orbitalCounts = ucell.get_orbitalCounts();
    this->set_atomCounts(atomCounts);
    this->set_orbitalCounts(orbitalCounts);
    this->set_npol(NPOL);
    // std::cout << "nw = " << this->get_nw() << std::endl;
    ModuleBase::Vector3<double>* sc_lambda;
    ModuleBase::Vector3<double>* init_mag;
    ModuleBase::Vector3<double>* target_mag;
    ModuleBase::Vector3<int>* constrain;
    double* decay_grad;
    int nat = this->get_nat();
    int ntype = this->get_ntype();
    if(GlobalV::MY_RANK == 0)
    {
        this->Set_ScData_From_Json(sc_file);
        this->set_sc_lambda();
        this->set_target_mag();
        this->set_constrain();
        this->set_decay_grad();
        sc_lambda = const_cast<ModuleBase::Vector3<double>*>(this->get_sc_lambda().data());
        target_mag = const_cast<ModuleBase::Vector3<double>*>(this->get_target_mag().data());
        constrain = const_cast<ModuleBase::Vector3<int>*>(this->get_constrain().data());
        decay_grad = const_cast<double*>(this->get_decay_grad().data());
    }
    else
    {
        sc_lambda = new ModuleBase::Vector3<double>[nat];
        target_mag = new ModuleBase::Vector3<double>[nat];
        constrain = new ModuleBase::Vector3<int>[nat];
        decay_grad = new double[ntype];
        ModuleBase::GlobalFunc::ZEROS(sc_lambda, nat);
        ModuleBase::GlobalFunc::ZEROS(target_mag, nat);
        ModuleBase::GlobalFunc::ZEROS(constrain, nat);
        ModuleBase::GlobalFunc::ZEROS(decay_grad, ntype);
    }
    for (int iat = 0; iat < nat; iat++)
    {
        Parallel_Common::bcast_double(sc_lambda[iat].x);
        Parallel_Common::bcast_double(sc_lambda[iat].y);
        Parallel_Common::bcast_double(sc_lambda[iat].z);
        Parallel_Common::bcast_double(target_mag[iat].x);
        Parallel_Common::bcast_double(target_mag[iat].y);
        Parallel_Common::bcast_double(target_mag[iat].z);
        Parallel_Common::bcast_int(constrain[iat].x);
        Parallel_Common::bcast_int(constrain[iat].y);
        Parallel_Common::bcast_int(constrain[iat].z);
    }
    for (int it = 0; it < ntype; it++)
    {
        Parallel_Common::bcast_double(decay_grad[it]);
    }
    if(GlobalV::MY_RANK != 0)
    {
        this->set_sc_lambda(sc_lambda, nat);
        this->set_target_mag(target_mag, nat);
        this->set_constrain(constrain, nat);
        this->set_decay_grad(decay_grad, ntype);
        delete [] sc_lambda;
        delete [] target_mag;
        delete [] constrain;
        delete[] decay_grad;
    }
}

// init sc
template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::init_sc(const UnitCell& ucell,
                                            int NPOL,
                                            std::string sc_file,
                                            Parallel_Orbitals* ParaV_in,
                                            int nspin_in,
                                            double sc_thr_in,
                                            int nsc_in,
                                            int nsc_min_in,
                                            double alpha_trial_in,
                                            double sccut_in,
                                            bool decay_grad_switch_in,
                                            K_Vectors kv_in,
                                            std::string KS_SOLVER_in,
                                            LCAO_Matrix* LM_in,
                                            hsolver::HSolver<std::complex<double>, psi::DEVICE_CPU>* phsol_in,
                                            hamilt::Hamilt<std::complex<double>, psi::DEVICE_CPU>* p_hamilt_in,
                                            psi::Psi<std::complex<double>>* psi_in,
                                            elecstate::ElecState* pelec_in)
{
    this->set_input_parameters(sc_thr_in, nsc_in, nsc_min_in, alpha_trial_in, sccut_in, decay_grad_switch_in);
    this->set_ParaV(ParaV_in);
    this->set_solver_parameters(nspin_in, kv_in, phsol_in, p_hamilt_in, psi_in, pelec_in, KS_SOLVER_in, LM_in);
    this->bcast_ScData(ucell, NPOL, sc_file);
}