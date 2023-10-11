#include "spin_constrain.h"
#include "module_base/parallel_common.h"

// init sc
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::init_sc(const UnitCell& ucell,
                            int NPOL,
                            std::string sc_file,
                            Parallel_Orbitals* ParaV_in,
                            int nspin_in,
                            double sc_thr_in,
                            int nsc_in,
                            int nsc_min_in,
                            K_Vectors kv_in,
                            std::string KS_SOLVER_in,
                            LCAO_Matrix* LM_in,
                            hsolver::HSolver<FPTYPE, Device>* phsol_in,
                            hamilt::Hamilt<FPTYPE, Device>* p_hamilt_in,
                            psi::Psi<FPTYPE>* psi_in,
                            elecstate::ElecState* pelec_in)
{
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
    int nat = this->get_nat();
    if(GlobalV::MY_RANK == 0)
    {
        this->Set_ScData_From_Json(GlobalV::sc_file);
        this->set_sc_lambda();
        this->set_target_mag();
        this->set_constrain();
        sc_lambda = const_cast<ModuleBase::Vector3<double>*>(this->get_sc_lambda().data());
        target_mag = const_cast<ModuleBase::Vector3<double>*>(this->get_target_mag().data());
        constrain = const_cast<ModuleBase::Vector3<int>*>(this->get_constrain().data());
    }
    else
    {
        sc_lambda = new ModuleBase::Vector3<double>[nat];
        target_mag = new ModuleBase::Vector3<double>[nat];
        constrain = new ModuleBase::Vector3<int>[nat];
        ModuleBase::GlobalFunc::ZEROS(sc_lambda, nat);
        ModuleBase::GlobalFunc::ZEROS(target_mag, nat);
        ModuleBase::GlobalFunc::ZEROS(constrain, nat);
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
    if(GlobalV::MY_RANK != 0)
    {
        this->set_sc_lambda(sc_lambda, nat);
        this->set_target_mag(target_mag, nat);
        this->set_constrain(constrain, nat);
        delete [] sc_lambda;
        delete [] target_mag;
        delete [] constrain;
    }
    // parameters for lambda loop
    this->sc_thr_ = sc_thr_in;
    this->nsc_ = nsc_in;
    this->nsc_min_ = nsc_min_in;
}

template class SpinConstrain<std::complex<double>, psi::DEVICE_CPU>;
