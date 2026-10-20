#include "source_lcao/module_lr/kernel_pots.h"
// Resolve LR dimensions, validate solver options, and initialize kernel potentials.
#include "esolver_lr_lcao_tddft.h"
#include "source_lcao/module_lr/init_helpers.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_io.h"
#include "source_io/module_wf/read_wfc_nao.h"
#include "source_io/module_parameter/parameter.h"
#include <algorithm>
#include <cmath>
#include <set>

using namespace LR;

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::set_dimension()
{
    this->nspin = this->inp_->nspin;
    this->nstates = this->inp_->lr_nstates;
    this->nbasis = PARAM.globalv.nlocal;
    int ks_nbands = this->inp_->nbands;
    if (this->nspin == 2)
    {
        this->nupdown = static_cast<int>(std::lround(this->inp_->nupdown));
        if (this->ks_)
        {
            this->nupdown = LR_Util::spin_population_difference(this->ks_->pelec->wg, this->kv.para_k.kpar, this->kv.para_k.rank_in_pool);
        }
        else if (this->inp_->ri_hartree_benchmark != "aims"
                 && this->inp_->ri_hartree_benchmark != "aims-librpa")
        {
            std::vector<double> populations;
            const bool gamma_only = std::is_same<T, double>::value;
            const bool binary = this->inp_->init_wfc_file_format == "binary";
            if (!ModuleIO::read_wfc_nao_spin_populations(this->in_dir, this->kv.get_nkstot(), this->nspin, gamma_only, binary,
                    this->my_rank_, populations))
            {
                ModuleBase::WARNING_QUIT("ESolver_LR", "read complete KS occupations failed");
            }
            this->nupdown = static_cast<int>(std::lround(populations[0]) - std::lround(populations[1]));
        }
    }
    // read_pseudo/ParamUpdater has already resolved nelec and applied nelec_delta.
    this->nocc_max = LR_Util::cal_nocc(this->inp_->nelec, this->nspin, this->nupdown);
    if (this->inp_->ri_hartree_benchmark == "aims" || this->inp_->ri_hartree_benchmark == "aims-librpa"
        && !this->inp_->aims_nbasis.empty())
    {
        // Total basis size: sum the per-atom basis size times the atom count for each type.
        int total_nbasis = 0;
        for (std::size_t it = 0; it < this->inp_->aims_nbasis.size(); ++it)
        {
            total_nbasis += this->inp_->aims_nbasis[it] * this->ucell_->atoms[it].na;
        }
        this->nbasis = total_nbasis;
        std::cout << "nbasis from aims: " << this->nbasis << std::endl;
        for (int it = 0; it < this->ucell_->ntype; ++it)
        {
            this->ucell_->atoms[it].nw = this->inp_->aims_nbasis[it];
        }
        const_cast<UnitCell*>(this->ucell_)->set_iat2iwt(1);  // update iat2iwt for aims_nbasis 25-05-23

        int nbands_file = 0;
        int nk_file = 0;
        int nspin_file = 0;
        int nocc_file = 0;
        LR_IO::parse_band_out_file(this->inp_->rpa_outdir, nbands_file, nk_file, nspin_file, nocc_file);
        std::cout << "nocc from band_out: " << nocc_file << std::endl;
        ks_nbands = nbands_file;
        this->nocc_max = nocc_file;
    }
    // calculate the number of occupied and unoccupied states
    // which determines the basis size of the excited states
    this->nocc_in = LR_Util::cal_nocc_window(this->inp_->nocc, this->nocc_max);
    this->nvirt_in = ks_nbands - this->nocc_max;   //nbands-nocc
    if (this->inp_->nvirt > this->nvirt_in)
    {
        this->ofs_running_ << "ESolver_LR: input nvirt is too large to cover by nbands, set nvirt = nbands - nocc = "
                           << this->nvirt_in << std::endl;
    }
    else if (this->inp_->nvirt > 0)
    {
        this->nvirt_in = this->inp_->nvirt;
    }
    this->nbands = this->nocc_in + this->nvirt_in;
    this->nk = this->inp_->nspin == 2 ? this->kv.get_nks() / 2 : this->kv.get_nks();
    this->nocc.resize(nspin, nocc_in);
    this->nvirt.resize(nspin, nvirt_in);
    if (this->nstates <= 0) {
        this->nstates = nk * nocc_in * nvirt_in;
        this->ofs_running_ << "ESolver_LR: lr_nstates <= 0, set nstates = nk * nocc * nvirt = " << this->nstates << std::endl;
    }
    for (int is = 0; is < nspin; ++is)
    {
        this->npairs.push_back(nocc[is] * nvirt[is]);
    }
    this->ofs_running_ << "Setting LR-TDDFT parameters: " << std::endl;
    this->ofs_running_ << "full occupied bands in largest spin channel: " << nocc_max << std::endl;
    const int skipped_core_bands = nocc_max - nocc_in;
    this->ofs_running_ << "shared core bands omitted from LR window: " << skipped_core_bands << std::endl;
    this->ofs_running_ << "number of occupied bands: " << nocc_in << std::endl;
    this->ofs_running_ << "number of virtual bands: " << nvirt_in << std::endl;
    this->ofs_running_ << "number of Atom orbitals (LCAO-basis size): " << this->nbasis << std::endl;
    // eig_ks is allocated after set_dimension; report the resolved counts here.
    this->ofs_running_ << "number of KS bands: " << ks_nbands << std::endl;
    this->ofs_running_ << "number of KS bands in LR window: " << this->nbands << std::endl;
    this->ofs_running_ << "number of excited states to be solved: " << this->nstates << std::endl;
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::reset_dim_spin2()
{
	if (nspin != 2)
	{
		return;
	}
	if (nupdown != 0)
    {
        this->openshell = true;
        const int minority_spin = nupdown > 0 ? 1 : 0;
        const int polarization = std::abs(nupdown);
        nocc[minority_spin] -= polarization;
        nvirt[minority_spin] += polarization;
        npairs = { nocc[0] * nvirt[0], nocc[1] * nvirt[1] };
        std::cout << "** Solve the spin-up and spin-down states separately for open-shell system. **" << std::endl;
    }
	for (int is : {0, 1})
	{
		if (npairs[is] <= 0)
		{
			throw std::invalid_argument(std::string("ESolver_LR: npairs (nocc*nvirt) <= 0 for spin") + std::string(is == 0 ? "up" : "down"));
		}
	}

	if (nstates > (npairs[0] + npairs[1]) * nk)
	{
		throw std::invalid_argument("ESolver_LR: nstates > nocc*nvirt*nk");
	}
	if (this->inp_->lr_unrestricted)
	{
		this->openshell = true;
	}
    if (!this->openshell)
	{
		std::cout << " ** Assuming degenerate spin-up and spin-down states  **" << std::endl;
	}
}

template<typename T, typename TR>
void ModuleESolver::ESolver_LR<T, TR>::init_pot(const Charge& chg_gs)
{
    this->pot.assign(nspin, nullptr);
    if (this->inp_->ri_hartree_benchmark != "none")
    {
        return;
    }
    const auto& grid = pgrid();
    auto prepared = LR::prepare_lr_potentials(*this->inp_, xc_kernel, *this->pw_rho,
        *this->ucell_, chg_gs, grid, nspin, openshell);
    this->pot = std::move(prepared.pot);
    if (this->inp_->cal_force)
    {
        this->init_pot_groundstate(chg_gs);
        this->pot_hxc_gs = LR::prepare_gs_potential(*this->inp_, xc_kernel, *this->pw_rho,
            *this->ucell_, chg_gs, grid, nspin, openshell, prepared);
    }
}

template class ModuleESolver::ESolver_LR<double, double>;
template class ModuleESolver::ESolver_LR<std::complex<double>, double>;
