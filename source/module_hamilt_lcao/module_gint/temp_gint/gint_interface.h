#pragma once
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint_type.h"


namespace ModuleGint
{

void cal_gint_vl(
    const double* vr_eff,
    HContainer<double>* hR);

void cal_gint_vl(
    std::vector<const double*> vr_eff,
    HContainer<std::complex<double>>* hR);

void cal_gint_vl_metagga(
    const double* vr_eff,
    const double* vfork,
    HContainer<double>* hR);

void cal_gint_vl_metagga(
    std::vector<const double*> vr_eff,
    std::vector<const double*> vofk,
    HContainer<std::complex<double>>* hR);

void cal_gint_rho(
    const std::vector<HContainer<double>*>& dm_vec,
    const int nspin,
    double **rho);

void cal_gint_tau(        
    const std::vector<HContainer<double>*>& dm_vec,
    const int nspin,
    double**tau);

void cal_gint_fvl(
    const int nspin,
    const std::vector<const double*>& vr_eff,
    const std::vector<HContainer<double>*>& dm_vec,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix* fvl,
    ModuleBase::matrix* svl);

void cal_gint_fvl_meta(
    const int nspin,
    const std::vector<const double*>& vr_eff,
    const std::vector<const double*>& vofk,
    const std::vector<HContainer<double>*>& dm_vec,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix* fvl,
    ModuleBase::matrix* svl);



} // namespace ModuleGint