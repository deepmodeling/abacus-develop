#include "init_helpers.h"
#include "utils/lr_util.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_lcao/lcao_nonlocal_info.h"
#include "source_cell/unitcell.h"
#include <set>
namespace LR
{
void setup_2center_table(TwoCenterBundle& two_center_bundle, LCAO_Orbitals& orb, UnitCell& ucell,
                         const Input_para& input, const int rank, std::ofstream& running_log)
{
    if (input.vnl_in_h)
    {
        auto* lcao_nl = new LCAONonlocalInfo();
        lcao_nl->setupNonlocal(ucell.ntype, ucell.atoms, running_log, orb,
                               input.basis_type, input.out_element_info,
                               input.lspinorb, input.nspin, rank);
        ucell.infoNL.reset(lcao_nl);
        two_center_bundle.build_beta(ucell.ntype, lcao_nl->get_nonlocal().get_Beta_data());
    }
    // NOTE: tabulate() must be called AFTER build_beta(), otherwise the
    // nonlocal (beta) two-center tables are left empty.
#ifdef __FFT_TWO_CENTER
    two_center_bundle.tabulate();
#else
    two_center_bundle.tabulate(input.lcao_ecut, input.lcao_dk, input.lcao_dr, input.lcao_rmax);
#endif
}

void check_lr_parameters(const Input_para& input, const std::string& xc_kernel, const int nspin)
{
    const std::set<std::string> lr_solvers = { "dav", "lapack" , "spectrum", "dav_subspace", "cg", "elpa", "plot" };
    // "rpa" and "bse" have no xc kernel at all; everything else is either a (semi)local
    // functional or a hybrid, both of which `LR_Util` enumerates. Listing the names a third
    // time here is what used to make a newly supported hybrid fail at input parsing.
    const std::set<std::string> kernel_less = { "rpa", "bse" };
    const std::set<std::string> abs_gauge = { "velocity", "length" };
    if (lr_solvers.find(input.lr_solver) == lr_solvers.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of lr_solver");
    }
    if (!kernel_less.count(xc_kernel) && !LR_Util::has_local_xc(xc_kernel)
        && !LR_Util::hybrid_xc_list().count(xc_kernel)) {
        throw std::invalid_argument("ESolver_LR: unknown type of xc_kernel");
    }
    if (nspin != 1 && nspin != 2) {
        throw std::invalid_argument("LR-TDDFT only supports nspin = 1 or 2 now");
    }
    if (abs_gauge.find(input.abs_gauge) == abs_gauge.end()) {
        throw std::invalid_argument("ESolver_LR: unknown type of abs_gauge");
    }
}

}
