#include "source_lcao/module_lr/lr_dim.h"

#include "source_base/matrix.h"
#include "source_base/tool_quit.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_io/module_wf/read_wfc_nao.h"
#include "source_lcao/module_lr/utils/lr_io.h"
#include "source_lcao/module_lr/utils/lr_util.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>

namespace LR
{

Dim resolve_dim(const DimSrc& src)
{
    Dim dim;
    const Input_para& inp = *src.inp;

    dim.nspin = inp.nspin;
    dim.nstates = inp.lr_nstates;
    dim.nbasis = src.nlocal;
    dim.ks_nbands = inp.nbands;
    if (dim.nspin == 2)
    {
        dim.nupdown = static_cast<int>(std::lround(inp.nupdown));
        if (src.wg_ks)
        {
            dim.nupdown = LR_Util::spin_population_difference(*src.wg_ks, src.kpar, src.rank_in_pool);
        }
        else if (inp.ri_hartree_benchmark != "aims" && inp.ri_hartree_benchmark != "aims-librpa")
        {
            std::vector<double> populations;
            const bool binary = inp.init_wfc_file_format == "binary";
            if (!ModuleIO::read_wfc_nao_spin_populations(
                    src.in_dir, src.nkstot, dim.nspin, src.gamma_only, binary, src.rank, populations))
            {
                ModuleBase::WARNING_QUIT("ESolver_LR", "read complete KS occupations failed");
            }
            dim.nupdown = static_cast<int>(std::lround(populations[0]) - std::lround(populations[1]));
        }
    }
    // read_pseudo/ParamUpdater has already resolved nelec and applied nelec_delta.
    dim.nocc_max = LR_Util::cal_nocc(inp.nelec, dim.nspin, dim.nupdown);
    if (inp.ri_hartree_benchmark == "aims"
        || (inp.ri_hartree_benchmark == "aims-librpa" && !inp.aims_nbasis.empty()))
    {
        resolve_aims_basis(*src.ucell, inp, dim);
    }
    // calculate the number of occupied and unoccupied states
    // which determines the basis size of the excited states
    dim.nocc_in = LR_Util::cal_nocc_window(inp.nocc, dim.nocc_max);
    dim.nvirt_in = dim.ks_nbands - dim.nocc_max; // nbands-nocc
    if (inp.nvirt > dim.nvirt_in)
    {
        if (src.log)
        {
            *src.log << "ESolver_LR: input nvirt is too large to cover by nbands, set nvirt = nbands - nocc = "
                     << dim.nvirt_in << std::endl;
        }
    }
    else if (inp.nvirt > 0)
    {
        dim.nvirt_in = inp.nvirt;
    }
    dim.nbands = dim.nocc_in + dim.nvirt_in;
    dim.nk = inp.nspin == 2 ? src.nks / 2 : src.nks;
    dim.nocc.assign(dim.nspin, dim.nocc_in);
    dim.nvirt.assign(dim.nspin, dim.nvirt_in);
    if (dim.nstates <= 0)
    {
        dim.nstates = dim.nk * dim.nocc_in * dim.nvirt_in;
        if (src.log)
        {
            *src.log << "ESolver_LR: lr_nstates <= 0, set nstates = nk * nocc * nvirt = " << dim.nstates << std::endl;
        }
    }
    dim.npairs.resize(dim.nspin);
    for (int is = 0; is < dim.nspin; ++is)
    {
        dim.npairs[is] = dim.nocc[is] * dim.nvirt[is];
    }
    if (src.log)
    {
        std::ostream& log = *src.log;
        log << "Setting LR-TDDFT parameters: " << std::endl;
        log << "full occupied bands in largest spin channel: " << dim.nocc_max << std::endl;
        const int skipped_core_bands = dim.nocc_max - dim.nocc_in;
        log << "shared core bands omitted from LR window: " << skipped_core_bands << std::endl;
        log << "number of occupied bands: " << dim.nocc_in << std::endl;
        log << "number of virtual bands: " << dim.nvirt_in << std::endl;
        log << "number of Atom orbitals (LCAO-basis size): " << dim.nbasis << std::endl;
        // eig_ks is allocated after resolve_dim; report the resolved counts here.
        log << "number of KS bands: " << dim.ks_nbands << std::endl;
        log << "number of KS bands in LR window: " << dim.nbands << std::endl;
        log << "number of excited states to be solved: " << dim.nstates << std::endl;
    }
    return dim;
}

void resolve_aims_basis(UnitCell& ucell, const Input_para& inp, Dim& dim)
{
    // Total basis size: sum the per-atom basis size times the atom count for each type.
    int total_nbasis = 0;
    for (std::size_t it = 0; it < inp.aims_nbasis.size(); ++it)
    {
        total_nbasis += inp.aims_nbasis[it] * ucell.atoms[it].na;
    }
    dim.nbasis = total_nbasis;
    std::cout << "nbasis from aims: " << dim.nbasis << std::endl;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        ucell.atoms[it].nw = inp.aims_nbasis[it];
    }
    ucell.set_iat2iwt(1); // update iat2iwt for aims_nbasis

    int nbands_file = 0;
    int nk_file = 0;
    int nspin_file = 0;
    int nocc_file = 0;
    LR_IO::parse_band_out_file(inp.rpa_outdir, nbands_file, nk_file, nspin_file, nocc_file);
    std::cout << "nocc from band_out: " << nocc_file << std::endl;
    dim.ks_nbands = nbands_file;
    dim.nocc_max = nocc_file;
}

void apply_spin_pol(Dim& dim, const Input_para& inp)
{
    if (dim.nspin != 2)
    {
        return;
    }
    if (dim.nupdown != 0)
    {
        dim.openshell = true;
        const int minority_spin = dim.nupdown > 0 ? 1 : 0;
        const int polarization = std::abs(dim.nupdown);
        dim.nocc[minority_spin] -= polarization;
        dim.nvirt[minority_spin] += polarization;
        dim.npairs = {dim.nocc[0] * dim.nvirt[0], dim.nocc[1] * dim.nvirt[1]};
        std::cout << "** Solve the spin-up and spin-down states separately for open-shell system. **" << std::endl;
    }
    for (int is : {0, 1})
    {
        if (dim.npairs[is] <= 0)
        {
            throw std::invalid_argument(
                std::string("ESolver_LR: npairs (nocc*nvirt) <= 0 for spin") + std::string(is == 0 ? "up" : "down"));
        }
    }

    if (dim.nstates > (dim.npairs[0] + dim.npairs[1]) * dim.nk)
    {
        throw std::invalid_argument("ESolver_LR: nstates > nocc*nvirt*nk");
    }
    if (inp.lr_unrestricted)
    {
        dim.openshell = true;
    }
    if (!dim.openshell)
    {
        std::cout << " ** Assuming degenerate spin-up and spin-down states  **" << std::endl;
    }
}

} // namespace LR
