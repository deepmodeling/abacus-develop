#include "source_estate/setup_estate_pw.h"
#include "source_estate/elecstate_pw.h"
#include "source_estate/elecstate_pw_sdft.h"
#include "source_estate/elecstate_tools.h"

namespace elecstate
{

void setup_estate_pw(
    UnitCell& ucell,
    K_Vectors& kv,
    Structure_Factor& sf,
    elecstate::ElecState*& pelec,
    Charge& chr,
    pseudopot_cell_vl& locpp,
    pseudopot_cell_vnl& ppcell,
    VSep*& vsep_cell,
    ModulePW::PW_Basis_K* pw_wfc,
    ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod,
    ModulePW::PW_Basis_Big* pw_big,
    surchem& solvent,
    const Input_para& inp)
{
    ModuleBase::TITLE("elecstate", "setup_estate_pw");

    const bool is_gpu = (inp.device == "gpu");
    const bool is_single = (inp.precision == "single");

#if ((defined __CUDA) || (defined __ROCM))
    if (is_gpu)
    {
        if (is_single)
        {
            setup_estate_pw_impl<std::complex<float>, base_device::DEVICE_GPU>(
                ucell, kv, sf, pelec, chr, locpp, ppcell, vsep_cell,
                pw_wfc, pw_rho, pw_rhod, pw_big, solvent, inp);
        }
        else
        {
            setup_estate_pw_impl<std::complex<double>, base_device::DEVICE_GPU>(
                ucell, kv, sf, pelec, chr, locpp, ppcell, vsep_cell,
                pw_wfc, pw_rho, pw_rhod, pw_big, solvent, inp);
        }
    }
    else
#endif
    {
        if (is_single)
        {
            setup_estate_pw_impl<std::complex<float>, base_device::DEVICE_CPU>(
                ucell, kv, sf, pelec, chr, locpp, ppcell, vsep_cell,
                pw_wfc, pw_rho, pw_rhod, pw_big, solvent, inp);
        }
        else
        {
            setup_estate_pw_impl<std::complex<double>, base_device::DEVICE_CPU>(
                ucell, kv, sf, pelec, chr, locpp, ppcell, vsep_cell,
                pw_wfc, pw_rho, pw_rhod, pw_big, solvent, inp);
        }
    }
}

template <typename T, typename Device>
void setup_estate_pw_impl(
    UnitCell& ucell,
    K_Vectors& kv,
    Structure_Factor& sf,
    elecstate::ElecState*& pelec,
    Charge& chr,
    pseudopot_cell_vl& locpp,
    pseudopot_cell_vnl& ppcell,
    VSep*& vsep_cell,
    ModulePW::PW_Basis_K* pw_wfc,
    ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod,
    ModulePW::PW_Basis_Big* pw_big,
    surchem& solvent,
    const Input_para& inp)
{
    if (pelec == nullptr)
    {
        if (inp.esolver_type == "sdft")
        {
            pelec = new elecstate::ElecStatePW_SDFT<std::complex<double>, Device>(pw_wfc,
                &chr, &kv, &ucell, &ppcell, pw_rho, pw_big);
        }
        else
        {
            pelec = new elecstate::ElecStatePW<T, Device>(pw_wfc,
                &chr, &kv, &ucell, &ppcell, pw_rho, pw_big);
        }
    }

    if (PARAM.inp.dfthalf_type > 0)
    {
        vsep_cell = new VSep;
        vsep_cell->init_vsep(*pw_rhod, ucell.sep_cell);
    }

    if (pelec->pot == nullptr)
    {
        pelec->pot = new elecstate::Potential(pw_rhod,
              pw_rho, &ucell, &locpp.vloc, &sf,
              &solvent, &(pelec->f_en.etxc), &(pelec->f_en.vtxc), vsep_cell);
    }

    locpp.init_vloc(ucell, pw_rhod);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    ppcell.init(ucell, &sf, pw_wfc);
    ppcell.init_vnl(ucell, pw_rhod);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    if (inp.ocp)
    {
        elecstate::fixed_weights(inp.ocp_kb,
                                 inp.nbands,
                                 inp.nelec,
                                 pelec->klist,
                                 pelec->wg,
                                 pelec->skip_weights);
    }
}

void teardown_estate_pw(elecstate::ElecState*& pelec, VSep*& vsep_cell)
{
    ModuleBase::TITLE("elecstate", "teardown_estate_pw");

    if (vsep_cell != nullptr)
    {
        delete vsep_cell;
    }

    if (pelec != nullptr)
    {
        delete pelec;
        pelec = nullptr;
    }
}

template <typename T, typename Device>
void teardown_estate_pw_impl(elecstate::ElecState*& pelec, VSep*& vsep_cell)
{
    ModuleBase::TITLE("elecstate", "teardown_estate_pw_impl");

    if (vsep_cell != nullptr)
    {
        delete vsep_cell;
    }

    if (pelec != nullptr)
    {
        auto* pw_elec = dynamic_cast<elecstate::ElecStatePW<T, Device>*>(pelec);
        if (pw_elec)
        {
            delete pw_elec;
            pelec = nullptr;
        }
        else
        {
            ModuleBase::WARNING_QUIT("elecstate::teardown_estate_pw_impl", "Invalid ElecState type");
        }
    }
}

template void setup_estate_pw_impl<std::complex<float>, base_device::DEVICE_CPU>(
    UnitCell& ucell, K_Vectors& kv, Structure_Factor& sf,
    elecstate::ElecState*& pelec, Charge& chr,
    pseudopot_cell_vl& locpp, pseudopot_cell_vnl& ppcell, VSep*& vsep_cell,
    ModulePW::PW_Basis_K* pw_wfc, ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod, ModulePW::PW_Basis_Big* pw_big,
    surchem& solvent, const Input_para& inp);

template void setup_estate_pw_impl<std::complex<double>, base_device::DEVICE_CPU>(
    UnitCell& ucell, K_Vectors& kv, Structure_Factor& sf,
    elecstate::ElecState*& pelec, Charge& chr,
    pseudopot_cell_vl& locpp, pseudopot_cell_vnl& ppcell, VSep*& vsep_cell,
    ModulePW::PW_Basis_K* pw_wfc, ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod, ModulePW::PW_Basis_Big* pw_big,
    surchem& solvent, const Input_para& inp);

template void teardown_estate_pw_impl<std::complex<float>, base_device::DEVICE_CPU>(
    elecstate::ElecState*& pelec, VSep*& vsep_cell);

template void teardown_estate_pw_impl<std::complex<double>, base_device::DEVICE_CPU>(
    elecstate::ElecState*& pelec, VSep*& vsep_cell);

#if ((defined __CUDA) || (defined __ROCM))

template void setup_estate_pw_impl<std::complex<float>, base_device::DEVICE_GPU>(
    UnitCell& ucell, K_Vectors& kv, Structure_Factor& sf,
    elecstate::ElecState*& pelec, Charge& chr,
    pseudopot_cell_vl& locpp, pseudopot_cell_vnl& ppcell, VSep*& vsep_cell,
    ModulePW::PW_Basis_K* pw_wfc, ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod, ModulePW::PW_Basis_Big* pw_big,
    surchem& solvent, const Input_para& inp);

template void setup_estate_pw_impl<std::complex<double>, base_device::DEVICE_GPU>(
    UnitCell& ucell, K_Vectors& kv, Structure_Factor& sf,
    elecstate::ElecState*& pelec, Charge& chr,
    pseudopot_cell_vl& locpp, pseudopot_cell_vnl& ppcell, VSep*& vsep_cell,
    ModulePW::PW_Basis_K* pw_wfc, ModulePW::PW_Basis* pw_rho,
    ModulePW::PW_Basis* pw_rhod, ModulePW::PW_Basis_Big* pw_big,
    surchem& solvent, const Input_para& inp);

template void teardown_estate_pw_impl<std::complex<float>, base_device::DEVICE_GPU>(
    elecstate::ElecState*& pelec, VSep*& vsep_cell);

template void teardown_estate_pw_impl<std::complex<double>, base_device::DEVICE_GPU>(
    elecstate::ElecState*& pelec, VSep*& vsep_cell);

#endif

}
