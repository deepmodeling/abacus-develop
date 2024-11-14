#include "FORCE.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_base/memory.h"
#include "module_parameter/parameter.h"
template<>
elecstate::DensityMatrix<double, double> Force_LCAO<double>::cal_edm(const elecstate::ElecState* pelec,
    const psi::Psi<double>& psi,
    const elecstate::DensityMatrix<double, double>& dm,
    const K_Vectors& kv,
    const Parallel_Orbitals& pv,
    const int& nspin, 
    const int& nbands,
    const UnitCell& ucell,
    Record_adj& ra) const
{
    ModuleBase::matrix wg_ekb;
    wg_ekb.create(nspin, nbands);

    for(int is=0; is<nspin; is++)
    {
        for(int ib=0; ib<nbands; ib++)
        {
            wg_ekb(is,ib) = pelec->wg(is,ib) * pelec->ekb(is, ib);
        }
    }

    // construct a DensityMatrix for Gamma-Only
    elecstate::DensityMatrix<double, double> edm(&pv, nspin);
    
#ifdef __PEXSI
    if (PARAM.inp.ks_solver == "pexsi")
    {
        auto pes = dynamic_cast<const elecstate::ElecStateLCAO<double>*>(pelec);
        for (int ik = 0; ik < nspin; ik++)
        {
            edm.set_DMK_pointer(ik, pes->get_DM()->pexsi_EDM[ik]);
        }
        
    }
    else
#endif
    {
        elecstate::cal_dm_psi(edm.get_paraV_pointer(), wg_ekb, psi, edm);
    }
    return edm;
}

template<>
elecstate::DensityMatrix<std::complex<double>, double> Force_LCAO<std::complex<double>>::cal_edm(const elecstate::ElecState* pelec,
    const psi::Psi<std::complex<double>>& psi,
    const elecstate::DensityMatrix<std::complex<double>, double>& dm,
    const K_Vectors& kv,
    const Parallel_Orbitals& pv,
    const int& nspin, 
    const int& nbands,
    const UnitCell& ucell,
    Record_adj& ra) const
{

    // construct a DensityMatrix object
    const int nspin_dm = nspin == 2 ? 2 : 1;
    elecstate::DensityMatrix<std::complex<double>, double> edm(&pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm);

    //--------------------------------------------
    // calculate the energy density matrix here.
    //--------------------------------------------

    ModuleBase::matrix wg_ekb;
    wg_ekb.create(kv.get_nks(), nbands);
    ModuleBase::Memory::record("Force::wg_ekb", sizeof(double) * kv.get_nks() * nbands);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
    for (int ik = 0; ik < kv.get_nks(); ik++)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            wg_ekb(ik, ib) = pelec->wg(ik, ib) * pelec->ekb(ik, ib);
        }
    }

    // use the original formula (Hamiltonian matrix) to calculate energy density matrix
    if (dm.EDMK.size())
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            edm.set_DMK_pointer(ik, dm.EDMK[ik].c);
        }
    }
    else
    {
        // cal_dm_psi
        elecstate::cal_dm_psi(edm.get_paraV_pointer(), wg_ekb, psi, edm);
    }

    // cal_dm_2d
    edm.init_DMR(ra, &ucell);
    edm.cal_DMR();
    return edm;
}