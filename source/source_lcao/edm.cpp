#include "edm.h"
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_estate/elecstate_lcao.h"
#include "source_base/memory_recorder.h"
#include "source_io/module_parameter/parameter.h"
template<>
module_dm::DensityMatrix<double, double> CalEDM<double>::cal_edm(const elecstate::ElecState* pelec,
    const psi::Psi<double>& psi,
    const module_dm::DensityMatrix<double, double>& dm,
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
    module_dm::DensityMatrix<double, double> edm(&pv, nspin);
    
#ifdef __PEXSI
    if (PARAM.inp.ks_solver == "pexsi")
    {
        // auto pes = dynamic_cast<const elecstate::ElecStateLCAO<double>*>(pelec);
        for (int ik = 0; ik < nspin; ik++)
        {
            edm.set_dmk_ptr(ik, dm.edm_pexsi[ik]);
        }
        
    }
    else
#endif
    {
        module_dm::dm_from_psi(&pv, wg_ekb, psi, edm);
    }
    edm.init_dmr(ra, &ucell);
    edm.cal_dmr(-1);
    return edm;
}

template<>
module_dm::DensityMatrix<std::complex<double>, double> CalEDM<std::complex<double>>::cal_edm(
    const elecstate::ElecState* pelec,
    const psi::Psi<std::complex<double>>& psi,
    const module_dm::DensityMatrix<std::complex<double>, double>& dm,
    const K_Vectors& kv,
    const Parallel_Orbitals& pv,
    const int& nspin, 
    const int& nbands,
    const UnitCell& ucell,
    Record_adj& ra) const
{

    // construct a DensityMatrix object
    // Pass the global physical nspin so that for SOC/noncollinear (nspin==4,
    // nspin_dm==1) cal_dmr selects the spin-resolved (Pauli) branch; otherwise
    // the overlap force/stress would use the real-projected (wrong) DMR.
    const int nspin_dm = nspin == 2 ? 2 : 1;
    module_dm::DensityMatrix<std::complex<double>, double> edm(&pv, nspin_dm, kv.kvec_d, kv.get_nks() / nspin_dm,
                                                               nspin);

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
    if (dm.edmk.size())
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            edm.set_dmk_ptr(ik, dm.edmk[ik].c);
        }
    }
    else
    {
        // dm_from_psi
        module_dm::dm_from_psi(&pv, wg_ekb, psi, edm);
    }

    // cal_dm_2d
    edm.init_dmr(ra, &ucell);
    edm.cal_dmr(-1);
    return edm;
}
