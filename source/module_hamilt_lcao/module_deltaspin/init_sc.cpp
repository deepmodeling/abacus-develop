#include "spin_constrain.h"
#include "module_base/parallel_common.h"

// init sc
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::init_sc(const UnitCell& ucell,
                                            int NPOL,
                                            std::string sc_file,
                                            Parallel_Orbitals* ParaV_in,
                                            int nspin_in,
                                            double sc_thr_in,
                                            int nsc_in,
                                            int nsc_min_in,
                                            double alpha_trial_in,
                                            bool decay_grad_switch_in,
                                            K_Vectors kv_in,
                                            std::string KS_SOLVER_in,
                                            LCAO_Matrix* LM_in,
                                            hsolver::HSolver<FPTYPE, Device>* phsol_in,
                                            hamilt::Hamilt<FPTYPE, Device>* p_hamilt_in,
                                            psi::Psi<FPTYPE>* psi_in,
                                            elecstate::ElecState* pelec_in)
{
    // input parameters for lambda loop
    this->sc_thr_ = sc_thr_in;
    this->nsc_ = nsc_in;
    this->nsc_min_ = nsc_min_in;
    this->decay_grad_switch_ = decay_grad_switch_in;
    this->alpha_trial_ = alpha_trial_in/ModuleBase::Ry_to_eV;
    // get pointer to outter pointers
    this->ParaV = ParaV_in;
    this->phsol = phsol_in;
    this->p_hamilt = p_hamilt_in;
    this->psi = psi_in;
    this->pelec = pelec_in;
    this->KS_SOLVER = KS_SOLVER_in;
    this->kv_ = kv_in;
    this->LM = LM_in;
    // get nloc
    int nloc = this->ParaV->nloc;
    if (nloc <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::init_sc", "nloc <= 0");
    }
    /// set nspin
    this->set_nspin(nspin_in);
    /// set ScData
    this->clear_ScData();
    this->clear_atomCounts();
    std::map<int, int> atomCounts = ucell.get_atomCounts();
    std::map<int, int> orbitalCounts = ucell.get_orbitalCounts();
    this->set_atomCounts(atomCounts);
    this->set_orbitalCounts(orbitalCounts);
    this->set_npol(GlobalV::NPOL);
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
        this->Set_ScData_From_Json(GlobalV::sc_file);
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

template class SpinConstrain<std::complex<double>, psi::DEVICE_CPU>;
template class SpinConstrain<double, psi::DEVICE_CPU>;
