#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_band_io.h"
#include "rpa_lri_detail.h"
#include "librpa_output_names.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>


template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_bands(const elecstate::ElecState* pelec)
try
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_bands");
    if (this->runtime.rank != 0)
    {
        return;
    }
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    const int nspin_tmp = this->runtime.input.nspin == 2 ? 2 : 1;
    RpaLriDetail::BandOutputData data;
    data.nk = nks_tot;
    data.nspin = nspin_tmp;
    data.nbands = this->runtime.input.nbands;
    data.nbasis = this->runtime.nlocal;
    data.fermi_ry = pelec->eferm.ef;
    data.samples.reserve(static_cast<std::size_t>(nks_tot) * static_cast<std::size_t>(nspin_tmp)
                         * static_cast<std::size_t>(this->runtime.input.nbands));
    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (int is = 0; is != nspin_tmp; is++)
        {
            for (int ib = 0; ib != this->runtime.input.nbands; ib++)
            {
                RpaLriDetail::BandSample sample;
                sample.occupation = pelec->wg(ik + is * nks_tot, ib) * nks_tot;
                sample.energy_ry = pelec->ekb(ik + is * nks_tot, ib);
                data.samples.push_back(sample);
            }
        }
    }
    std::ofstream ofs(
        RpaLriDetail::librpa_output_path(this->runtime.input.rpa_outdir, "band_out.txt"),
        std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open band_out.txt");
    }
    RpaLriDetail::write_band_output(ofs, data);
    ofs.close();
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to write band_out.txt");
    }
}
catch (const std::exception& error)
{
    RpaLriDetail::abort_output(this->mpi_comm, error.what());
}
catch (...)
{
    RpaLriDetail::abort_output(this->mpi_comm, "Unknown exception in out_bands");
}

template void RPA_LRI<double, double>::out_bands(const elecstate::ElecState* pelec);
template void RPA_LRI<std::complex<double>, double>::out_bands(const elecstate::ElecState* pelec);
