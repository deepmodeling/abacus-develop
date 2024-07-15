#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::set_globalv_bcast()
{
    add_bool_bcast(sys.two_fermi);
    add_int_bcast(sys.nrcut);
    add_double_bcast(sys.bessel_nao_rcut);
    add_bool_bcast(sys.dos_setemin);
    add_bool_bcast(sys.dos_setemax);
    add_int_bcast(sys.ncx);
    add_int_bcast(sys.ncy);
    add_int_bcast(sys.ncz);
    add_bool_bcast(sys.out_md_control);
    add_bool_bcast(sys.rpa_setorb);
    add_bool_bcast(sys.gamma_only_local);
    add_bool_bcast(sys.double_grid);
    add_double_bcast(sys.uramping);
    add_string_bcast(sys.global_calculation);
}
} // namespace ModuleIO