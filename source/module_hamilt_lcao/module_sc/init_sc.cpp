#include "spin_constrain.h"
#include "module_base/parallel_common.h"

// init sc
void SpinConstrain::init_sc(const UnitCell& ucell, int NPOL, std::string sc_file, Parallel_Orbitals* ParaV_in, int nspin_in)
{
    // get pointer to Parallel_Oribitals
    this->ParaV = ParaV_in;
    int nloc = this->ParaV->nloc;
    if (nloc <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::init_sc", "nloc <= 0");
    }
    // initialize Wi_, which is the weight function with size nloc
    this->Wi_.resize(nloc);
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
    ModuleBase::Vector3<double>* sc_mag;
    int nat = this->get_nat();
    if(GlobalV::MY_RANK == 0)
    {
        this->Set_ScData_From_Json(GlobalV::sc_file);
        this->set_sc_lambda();
        this->set_sc_mag();
        sc_lambda = const_cast<ModuleBase::Vector3<double>*>(this->get_sc_lambda().data());
        sc_mag = const_cast<ModuleBase::Vector3<double>*>(this->get_sc_mag().data());
    }
    else
    {
        sc_lambda = new ModuleBase::Vector3<double>[nat];
        sc_mag = new ModuleBase::Vector3<double>[nat];
        ModuleBase::GlobalFunc::ZEROS(sc_lambda, nat);
        ModuleBase::GlobalFunc::ZEROS(sc_mag, nat);
    }
    for (int iat = 0; iat < nat; iat++)
    {
        Parallel_Common::bcast_double(sc_lambda[iat].x);
        Parallel_Common::bcast_double(sc_lambda[iat].y);
        Parallel_Common::bcast_double(sc_lambda[iat].z);
        Parallel_Common::bcast_double(sc_mag[iat].x);
        Parallel_Common::bcast_double(sc_mag[iat].y);
        Parallel_Common::bcast_double(sc_mag[iat].z);
    }
    if(GlobalV::MY_RANK != 0)
    {
        this->set_sc_lambda(sc_lambda, nat);
        this->set_sc_mag(sc_mag, nat);
        delete [] sc_lambda;
        delete [] sc_mag;
    }
}