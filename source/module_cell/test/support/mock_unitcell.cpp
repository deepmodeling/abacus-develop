#include "module_cell/unitcell.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
/*
    README:
    This file supports idea like "I dont need any functions of UnitCell, I want
   to avoid using UnitCell functions because there is GLobalC, which will bring
   endless compile troubles like undefined behavior"
*/
void UnitCell::set_iat2iwt(const int& npol_in) {}
UnitCell::UnitCell() {
    itia2iat.create(1, 1);
}
UnitCell::~UnitCell() {
    delete[] atom_label;
    delete[] atom_mass;
    delete[] pseudo_fn;
    delete[] pseudo_type;
    delete[] orbital_fn;
    if (set_atom_flag) {
        delete[] atoms;
    }
}
void UnitCell::print_cell(std::ofstream& ofs) const {}
void UnitCell::print_cell_xyz(const std::string& fn) const {}
int UnitCell::read_atom_species(std::ifstream& ifa,
                                std::ofstream& ofs_running) {
    return 0;
}
bool UnitCell::read_atom_positions(std::ifstream& ifpos,
                                   std::ofstream& ofs_running,
                                   std::ofstream& ofs_warning) {
    return true;
}
void UnitCell::update_pos_taud(double* posd_in) {}
void UnitCell::update_pos_taud(const ModuleBase::Vector3<double>* posd_in) {}
void UnitCell::update_vel(const ModuleBase::Vector3<double>* vel_in) {}
void UnitCell::bcast_atoms_tau() {}
bool UnitCell::judge_big_cell() const { return true; }
void UnitCell::update_stress(ModuleBase::matrix& scs) {}
void UnitCell::update_force(ModuleBase::matrix& fcs) {}
#ifdef __MPI
void UnitCell::bcast_unitcell() {}
void UnitCell::bcast_unitcell2() {}
#endif
void UnitCell::set_iat2itia() {}
void UnitCell::setup_cell(const std::string& fn, std::ofstream& log) {}
void UnitCell::read_orb_file(int it,
                             std::string& orb_file,
                             std::ofstream& ofs_running,
                             Atom* atom) {}
int UnitCell::find_type(const std::string& label) { return 0; }
void UnitCell::print_tau() const {}
void UnitCell::print_stru_file(const std::string& fn,
                               const int& nspin,
                               const bool& direct,
                               const bool& vel,
                               const bool& magmom,
                               const bool& orb,
                               const bool& dpks_desc,
                               const int& iproc) const {}
void UnitCell::check_dtau() {}
void UnitCell::cal_nwfc(std::ofstream& log) {}
void UnitCell::cal_meshx() {}
void UnitCell::cal_natomwfc(std::ofstream& log) {}
bool UnitCell::check_tau() const { return true; }
bool UnitCell::if_atoms_can_move() const { return true; }
bool UnitCell::if_cell_can_change() const { return true; }
void UnitCell::setup(const std::string& latname_in,
                     const int& ntype_in,
                     const int& lmaxmax_in,
                     const bool& init_vel_in,
                     const std::string& fixed_axes_in) {}
void cal_nelec(const Atom* atoms, const int& ntype, double& nelec) {}
void UnitCell::compare_atom_labels(std::string label1, std::string label2) {}