#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"

#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_bz_sampling(const UnitCell& ucell)
try
{
    if (this->runtime.rank != 0)
    {
        return;
    }

    ModuleBase::TITLE("DFT_RPA_interface", "out_bz_sampling");
    const int nks_tot = this->runtime.input.nspin == 2 ? static_cast<int>(p_kv->get_nks()) / 2 : p_kv->get_nks();
    const int n_coulomb_irreducible = nks_tot;
    const int nk_full = p_kv->get_nkstot_nospin();
    if (nk_full < nks_tot || nk_full != p_kv->nmp[0] * p_kv->nmp[1] * p_kv->nmp[2]
        || p_kv->kvec_c_full.size() != static_cast<std::size_t>(nk_full)
        || p_kv->ibz_index.size() != static_cast<std::size_t>(nk_full))
    {
        throw std::runtime_error("Inconsistent full-grid metadata for bz_sample.txt");
    }

    std::ofstream ofs(this->runtime.input.rpa_outdir + "bz_sample.txt", std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open bz_sample.txt");
    }
    ofs << p_kv->nmp[0] << std::setw(6) << p_kv->nmp[1] << std::setw(6) << p_kv->nmp[2] << std::endl;
    ofs << nks_tot << std::setw(8) << n_coulomb_irreducible << std::endl;
    double weight_sum = 0.0;
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        weight_sum += p_kv->wk[ik];
    }
    if (!std::isfinite(weight_sum) || weight_sum <= 0.0)
    {
        throw std::runtime_error("Cannot write bz_sample.txt with non-positive total k-point weight.");
    }
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        const ModuleBase::Vector3<double> kvec_cartesian = p_kv->kvec_c[ik] * (ModuleBase::TWO_PI / ucell.lat0);
        ofs << std::setw(8) << ik + 1 << std::setw(24) << std::scientific << std::setprecision(15)
            << (p_kv->wk[ik] / weight_sum) << std::setw(24) << std::scientific << std::setprecision(15)
            << p_kv->kvec_d[ik].x << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_d[ik].y
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_d[ik].z << std::setw(24)
            << std::scientific << std::setprecision(15) << kvec_cartesian.x << std::setw(24) << std::scientific
            << std::setprecision(15) << kvec_cartesian.y << std::setw(24) << std::scientific << std::setprecision(15)
            << kvec_cartesian.z << std::setw(8) << ik + 1 << std::setw(8) << ik + 1 << std::endl;
    }
    // SCF, KS and Coulomb files share the reduced point numbering above.
    // The full list below explicitly connects that numbering to the BvK grid.
    // Keep established full-grid files byte-compatible (no optional tail).
    if (nk_full > nks_tot)
    {
        std::vector<int> multiplicity(nks_tot, 0);
        ofs << "full_kmap " << nk_full << '\n';
        for (int ik = 0; ik < nk_full; ++ik)
        {
            const int rep = p_kv->ibz_index[ik];
            if (rep < 0 || rep >= nks_tot)
            {
                throw std::runtime_error("Full-grid representative is outside the SCF list");
            }
            ++multiplicity[rep];
            const auto kcart = p_kv->kvec_c_full[ik] * (ModuleBase::TWO_PI / ucell.lat0);
            ofs << std::setw(8) << ik + 1 << std::setw(24) << kcart.x << std::setw(24) << kcart.y
                << std::setw(24) << kcart.z << std::setw(8) << rep + 1 << '\n';
        }
        for (int ik = 0; ik < nks_tot; ++ik)
        {
            if (multiplicity[ik] == 0
                || std::abs(p_kv->wk[ik] / weight_sum - static_cast<double>(multiplicity[ik]) / nk_full) > 1e-6)
            {
                throw std::runtime_error("Folded SCF weight disagrees with full-grid multiplicity");
            }
        }
    }
    ofs.close();
    if (!ofs.good()) throw std::runtime_error("Failed to write bz_sample.txt");
}
catch (const std::exception& error)
{
    RpaLriDetail::abort_output(this->mpi_comm, error.what());
}
catch (...)
{
    RpaLriDetail::abort_output(this->mpi_comm, "Unknown exception in out_bz_sampling");
}

template void RPA_LRI<double, double>::out_bz_sampling(const UnitCell& ucell);
template void RPA_LRI<std::complex<double>, double>::out_bz_sampling(const UnitCell& ucell);
