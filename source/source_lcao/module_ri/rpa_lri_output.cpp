#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
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
void RPA_LRI<T, Tdata>::out_struc(const UnitCell& ucell)
{
    if (this->runtime.rank != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_struc");
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    const ModuleBase::Matrix3 lat = ucell.latvec * ucell.lat0;                     // in unit of Bohr
    const ModuleBase::Matrix3 G_RPA = ucell.G * (ModuleBase::TWO_PI / ucell.lat0); // in unit of 1/Bohr
    std::ofstream ofs;
    ofs.open(this->runtime.input.rpa_outdir + "stru_out.txt", std::ios::out);
    ofs << std::fixed << std::setprecision(9);
    ofs << lat.e11 << std::setw(15) << lat.e12 << std::setw(15) << lat.e13 << std::endl;
    ofs << lat.e21 << std::setw(15) << lat.e22 << std::setw(15) << lat.e23 << std::endl;
    ofs << lat.e31 << std::setw(15) << lat.e32 << std::setw(15) << lat.e33 << std::endl;

    ofs << G_RPA.e11 << std::setw(15) << G_RPA.e12 << std::setw(15) << G_RPA.e13 << std::endl;
    ofs << G_RPA.e21 << std::setw(15) << G_RPA.e22 << std::setw(15) << G_RPA.e23 << std::endl;
    ofs << G_RPA.e31 << std::setw(15) << G_RPA.e32 << std::setw(15) << G_RPA.e33 << std::endl;

    ofs << ucell.nat << std::endl;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const ModuleBase::Vector3<double> position = ucell.atoms[it].tau[ia] * ucell.lat0; // in unit of Bohr
            ofs << std::setw(15) << position.x << std::setw(15) << position.y << std::setw(15) << position.z
                << std::setw(15) << it + 1 << std::endl;
        }
    }

    // Reader-v1 stores k-point sampling in bz_sample.txt. Keep the legacy
    // k-point block in stru_out.txt for reader-v0 consumers.
    if (this->runtime.input.out_librpa_ver != 1)
    {
        ofs << p_kv->nmp[0] << std::setw(6) << p_kv->nmp[1] << std::setw(6) << p_kv->nmp[2] << std::setw(6)
            << std::endl;

        for (int ik = 0; ik != nks_tot; ik++)
        {
            const ModuleBase::Vector3<double> kpoint
                = p_kv->kvec_c[ik] * (ModuleBase::TWO_PI / ucell.lat0); // in unit of 1/Bohr
            ofs << std::setw(15) << kpoint.x << std::setw(15) << kpoint.y << std::setw(15) << kpoint.z << std::endl;
        }
        if (this->runtime.input.symmetry == "-1")
        {
            for (int ik = 0; ik != nks_tot; ++ik)
            {
                ofs << (ik + 1) << std::endl;
            }
        }
    }
    ofs.close();
    return;
}

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
    if (weight_sum <= 0.0)
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
    ofs.close();
}
catch (const std::exception& error)
{
    RpaLriDetail::abort_output(this->mpi_comm, error.what());
}
catch (...)
{
    RpaLriDetail::abort_output(this->mpi_comm, "Unknown exception in out_bz_sampling");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_bands(const elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_bands");
    if (this->runtime.rank != 0)
    {
        return;
    }
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    const int nspin_tmp = this->runtime.input.nspin == 2 ? 2 : 1;
    std::ofstream ofs;
    ofs.open(this->runtime.input.rpa_outdir + "band_out.txt", std::ios::out);
    ofs << std::fixed << std::setprecision(15);
    ofs << nks_tot << std::endl;
    ofs << nspin_tmp << std::endl;
    ofs << this->runtime.input.nbands << std::endl;
    ofs << this->runtime.nlocal << std::endl;
    ofs << (pelec->eferm.ef / 2.0) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (int is = 0; is != nspin_tmp; is++)
        {
            ofs << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << std::endl;
            for (int ib = 0; ib != this->runtime.input.nbands; ib++)
            {
                ofs << std::setw(5) << ib + 1 << "   " << std::setw(8) << pelec->wg(ik + is * nks_tot, ib) * nks_tot
                    << std::setw(25) << pelec->ekb(ik + is * nks_tot, ib) / 2.0 << std::setw(25)
                    << pelec->ekb(ik + is * nks_tot, ib) * ModuleBase::Ry_to_eV << std::endl;
            }
        }
    }
    ofs.close();
    return;
}

template void RPA_LRI<double, double>::out_struc(const UnitCell& ucell);
template void RPA_LRI<std::complex<double>, double>::out_struc(const UnitCell& ucell);
template void RPA_LRI<double, double>::out_bz_sampling(const UnitCell& ucell);
template void RPA_LRI<std::complex<double>, double>::out_bz_sampling(const UnitCell& ucell);
template void RPA_LRI<double, double>::out_bands(const elecstate::ElecState* pelec);
template void RPA_LRI<std::complex<double>, double>::out_bands(const elecstate::ElecState* pelec);
