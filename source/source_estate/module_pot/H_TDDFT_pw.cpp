#include "H_TDDFT_pw.h"

#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_io/module_efield/td_efield_io.h"
#include "source_io/module_parameter/parameter.h"
#include "td_field_manager.h"

namespace elecstate
{

int H_TDDFT_pw::stype = 0;
ModuleBase::Vector3<double> H_TDDFT_pw::At;
ModuleBase::Vector3<double> H_TDDFT_pw::At_laststep;
ModuleBase::Vector3<double> H_TDDFT_pw::Et;
std::vector<double> H_TDDFT_pw::global_vext_time = {0.0, 0.0, 0.0};

H_TDDFT_pw::H_TDDFT_pw(const ModulePW::PW_Basis* rho_basis_in,
                       const UnitCell* ucell_in,
                       const std::shared_ptr<TDFieldManager>& field_manager)
    : ucell_(ucell_in), field_manager_(field_manager)
{
    this->dynamic_mode = false;
    this->fixed_mode = true;
    this->rho_basis_ = rho_basis_in;

    if (!field_manager_)
    {
        ModuleBase::WARNING_QUIT("H_TDDFT_pw", "RT-TDDFT field manager is not initialized.");
    }
    sync_compatibility_state(*field_manager_);
}

void H_TDDFT_pw::sync_compatibility_state(const TDFieldManager& manager)
{
    stype = manager.gauge();
    At = manager.vector_potential();
    At_laststep = manager.vector_potential_laststep();
    Et = manager.electric_field();
    const ModuleBase::Vector3<double>& total_field = manager.total_electric_field();
    global_vext_time = {total_field[0], total_field[1], total_field[2]};
}

void H_TDDFT_pw::cal_fixed_v(double* vl_pseudo)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_fixed_v");
    if (field_manager_->gauge() != 0)
    {
        return;
    }

    // Advance exactly once per rebuilt fixed potential. The potential then
    // consumes the same per-occurrence samples exposed to field output.
    field_manager_->advance_length_gauge();
    sync_compatibility_state(*field_manager_);
    if (!field_manager_->active())
    {
        return;
    }

    ModuleBase::timer::start("H_TDDFT_pw", "cal_fixed_v");
    const std::vector<TDField>& fields = field_manager_->fields();
    const std::vector<double>& field_values = field_manager_->field_values();
    for (std::size_t field_index = 0; field_index < fields.size(); ++field_index)
    {
        std::vector<double> vext_space(this->rho_basis_->nrxx, 0.0);
        const double field_value = field_values[field_index];

        cal_v_space_length(vext_space, fields[field_index].direction() + 1);
        for (std::size_t ir = 0; ir < static_cast<std::size_t>(this->rho_basis_->nrxx); ++ir)
        {
            vl_pseudo[ir] += vext_space[ir] * field_value;
        }
    }
    if (PARAM.inp.out_efield && GlobalV::MY_RANK == 0)
    {
        ModuleIO::write_td_field_values(*field_manager_, PARAM.globalv.global_out_dir);
    }
    ModuleBase::timer::end("H_TDDFT_pw", "cal_fixed_v");
}

void H_TDDFT_pw::cal_v_space_length(std::vector<double>& vext_space, const int direction)
{
    ModuleBase::TITLE("H_TDDFT_pw", "cal_v_space_length");
    ModuleBase::timer::start("H_TDDFT_pw", "cal_v_space_length");

    for (int ir = 0; ir < this->rho_basis_->nrxx; ++ir)
    {
        const int i = ir / (this->rho_basis_->ny * this->rho_basis_->nplane);
        const int j = ir / this->rho_basis_->nplane - i * this->rho_basis_->ny;
        const int k = ir % this->rho_basis_->nplane + this->rho_basis_->startz_current;
        const double x = static_cast<double>(i) / this->rho_basis_->nx;
        const double y = static_cast<double>(j) / this->rho_basis_->ny;
        const double z = static_cast<double>(k) / this->rho_basis_->nz;

        if (direction == 1)
        {
            vext_space[ir] = cal_v_space_length_potential(x) * this->ucell_->latvec.e11
                             + cal_v_space_length_potential(y) * this->ucell_->latvec.e21
                             + cal_v_space_length_potential(z) * this->ucell_->latvec.e31;
        }
        else if (direction == 2)
        {
            vext_space[ir] = cal_v_space_length_potential(x) * this->ucell_->latvec.e12
                             + cal_v_space_length_potential(y) * this->ucell_->latvec.e22
                             + cal_v_space_length_potential(z) * this->ucell_->latvec.e32;
        }
        else
        {
            vext_space[ir] = cal_v_space_length_potential(x) * this->ucell_->latvec.e13
                             + cal_v_space_length_potential(y) * this->ucell_->latvec.e23
                             + cal_v_space_length_potential(z) * this->ucell_->latvec.e33;
        }
    }

    ModuleBase::timer::end("H_TDDFT_pw", "cal_v_space_length");
}

double H_TDDFT_pw::cal_v_space_length_potential(const double coordinate) const
{
    const double lower_cut = field_manager_->length_cut1();
    const double upper_cut = field_manager_->length_cut2();
    if (coordinate < lower_cut)
    {
        return -((coordinate - lower_cut) * (upper_cut - lower_cut) / (lower_cut + 1.0 - upper_cut) - lower_cut) * this->ucell_->lat0;
    }
    if (coordinate < upper_cut)
    {
        return coordinate * this->ucell_->lat0;
    }
    return -((coordinate - upper_cut) * (upper_cut - lower_cut) / (lower_cut + 1.0 - upper_cut) - upper_cut) * this->ucell_->lat0;
}

void H_TDDFT_pw::compute_force(const UnitCell& cell, ModuleBase::matrix& force)
{
    int atom_index = 0;
    for (int type = 0; type < cell.ntype; ++type)
    {
        for (int atom = 0; atom < cell.atoms[type].na; ++atom)
        {
            for (int direction = 0; direction < 3; ++direction)
            {
                force(atom_index, direction) = global_vext_time[direction] * cell.atoms[type].ncpp.zv;
            }
            ++atom_index;
        }
    }
}

} // namespace elecstate
