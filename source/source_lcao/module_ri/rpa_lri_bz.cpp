#include "rpa_lri.h"
#include "rpa_lri_bz_io.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"

#include <array>
#include <complex>
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
    RpaLriDetail::BzSampleData data;
    data.mesh = {{p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]}};
    data.reduced.resize(nks_tot);
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        const ModuleBase::Vector3<double> kcart = p_kv->kvec_c[ik] * (ModuleBase::TWO_PI / ucell.lat0);
        data.reduced[ik].weight = p_kv->wk[ik];
        data.reduced[ik].direct = {{p_kv->kvec_d[ik].x, p_kv->kvec_d[ik].y, p_kv->kvec_d[ik].z}};
        data.reduced[ik].cartesian = {{kcart.x, kcart.y, kcart.z}};
    }
    data.full_cartesian.resize(p_kv->kvec_c_full.size());
    for (std::size_t ik = 0; ik < p_kv->kvec_c_full.size(); ++ik)
    {
        const auto kcart = p_kv->kvec_c_full[ik] * (ModuleBase::TWO_PI / ucell.lat0);
        data.full_cartesian[ik] = {{kcart.x, kcart.y, kcart.z}};
    }
    data.full_map = p_kv->ibz_index;
    RpaLriDetail::write_bz_sample(this->runtime.input.rpa_outdir + "bz_sample.txt", data);
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
