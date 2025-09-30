#include "source_cell/unitcell.h"
#define private public
#include "source_io/module_parameter/parameter.h"
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
    if (set_atom_flag) {
        delete[] atoms;
    }
}
SepPot::SepPot(){}
SepPot::~SepPot(){}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

void UnitCell::print_cell(std::ofstream& ofs) const {}

void UnitCell::set_iat2itia() {}

void UnitCell::setup_cell(const std::string& fn, std::ofstream& log) {}

bool UnitCell::if_atoms_can_move() const { return true; }

bool UnitCell::if_cell_can_change() const { return true; }

void UnitCell::setup(const std::string& latname_in,
                     const int& ntype_in,
                     const int& lmaxmax_in,
                     const bool& init_vel_in,
                     const std::string& fixed_axes_in) {}

void cal_nelec(const Atom* atoms, const int& ntype, double& nelec) {}

void UnitCell::compare_atom_labels(const std::string &label1, const std::string &label2) {}
