#include "mom_op_mat.h"

#include "source_base/basis_flattened_indexing.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/timer.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_basis/module_ao/orb_read.h"
#include "source_cell/unitcell.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_io/module_parameter/input_parameter.h"
#include "velocity_op.h"

#include "source_io/module_hs/hsr_writer.h"

#include <complex>
#include <string>
#include <vector>

namespace Momentum_op
{

namespace
{

/// @brief Build the CSR filename for one momentum component.
/// @param axis        'x', 'y' or 'z'
/// @param output_dir  directory prefix (e.g. "MAT/")
/// @param istep       zero-based ionic step index
/// @param out_app_flag 0 = fixed name overwritten, 1 = fixed name appended, 2 = per-step g{step} name
std::string make_filename(const char axis,
                          const std::string& output_dir,
                          const int istep,
                          const int out_app_flag)
{
    std::string base = output_dir + "pr";
    base.push_back(axis);
    if (out_app_flag == 2 && istep >= 0)
    {
        base += "g" + std::to_string(istep + 1);
    }
    base += "_nao.csr";
    return base;
}

} // anonymous namespace

void calculate(const int istep,
               const UnitCell& ucell,
               Grid_Driver& gd,
               TwoCenterBundle& two_center_bundle,
               const Parallel_Orbitals& pv,
               const LCAO_Orbitals& orb,
               const std::string& output_dir,
               const Input_para& inp)
{
    ModuleBase::TITLE("Momentum_op", "calculate");
    ModuleBase::timer::start("Momentum_op", "calculate");

    if (inp.out_mat_p.empty() || inp.out_mat_p[0] == 0)
    {
        ModuleBase::timer::end("Momentum_op", "calculate");
        return;
    }

    int rank = 0;
#ifdef __MPI
    const MPI_Comm comm = pv.para_2d().comm();
    if (comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm, &rank);
    }
#endif
    if (rank != 0)
    {
        ModuleBase::timer::end("Momentum_op", "calculate");
        return;
    }

    // Assemble on rank 0 in a serial layout. A distributed HContainer owns
    // only rank-local matrix elements and cannot be written directly as the
    // complete CSR matrix.
    Parallel_Blkcyclic serial_layout;
    serial_layout.set_serial(pv.para_2d().nrglb(), pv.para_2d().ncglb());

    // --- Convert the two-center bundle into the legacy LCAO_Orbitals layout ---
    LCAO_Orbitals lcao_orb;
    two_center_bundle.to_LCAO_Orbitals(lcao_orb,
                                       inp.lcao_ecut,
                                       inp.lcao_dk,
                                       inp.lcao_dr,
                                       inp.lcao_rmax,
                                       inp.out_element_info,
                                       inp.cal_force);

    // --- Build the flattened basis indexing used by Velocity_op ---
    const ModuleBase::Element_Basis_Index::Range range =
        ModuleBase::Element_Basis_Index::construct_range(lcao_orb);
    const ModuleBase::Element_Basis_Index::IndexLNM hierarchy =
        ModuleBase::Element_Basis_Index::construct_index(range);
    const BasisFlattenedIndexing basis(hierarchy, ucell.iat2it, ucell.get_npol());

    // --- Construct the velocity (momentum) operator kernel ---
    Velocity_op<std::complex<double>> vel_op(
        &ucell,
        &gd,
        basis,
        serial_layout,
        lcao_orb,
        two_center_bundle.overlap_orb.get());

    // --- Compute the gradient (kinetic) term and the nonlocal commutator term.
    //     Both accumulate into the same per-component HContainer inside Velocity_op. ---
    vel_op.calculate_grad_term();
    vel_op.calculate_vcomm_r();

    // --- Write each Cartesian component ---
    const char axes[3] = {'x', 'y', 'z'};
    const std::string label[3] = {
        "Momentum operator px matrix (a.u.)",
        "Momentum operator py matrix (a.u.)",
        "Momentum operator pz matrix (a.u.)"};
    const std::string note
        = "p = -i nabla + i [V_nl, r], Hermitian operator in the local atomic basis; "
          "only the scalar (npol = 1) block is written.";

    const int precision = inp.out_mat_p.size() > 1 ? inp.out_mat_p[1] : 8;
    const int out_app_flag = inp.out_app_flag;

    for (int c = 0; c < 3; ++c)
    {
        hamilt::HContainer<std::complex<double>>* comp = vel_op.get_current_term_pointer(c);
        if (comp == nullptr)
        {
            ModuleBase::WARNING_QUIT("Momentum_op::calculate",
                                     "Velocity_op returned a null HContainer for component " + std::to_string(c));
        }

        const std::string filename = make_filename(axes[c], output_dir, istep, out_app_flag);
        const bool append = (out_app_flag == 1) && (istep > 0);
        // get_s passes istep = -1; the CSR header prints one-based step indices,
        // so normalize it to the single "first" record like the other matrices.
        const int header_istep = istep >= 0 ? istep : 0;

        ModuleIO::write_hcontainer_csr(filename,
                                       &ucell,
                                       precision,
                                       comp,
                                       header_istep,
                                       append,
                                       0,       // ispin: momentum is spin-independent, single scalar block
                                       1,       // nspin: single scalar block written
                                       label[c],
                                       note,
                                       0.0,     // efermi_eV: not applicable
                                       false);  // has_efermi
    }

    ModuleBase::timer::end("Momentum_op", "calculate");
}

} // namespace Momentum_op
