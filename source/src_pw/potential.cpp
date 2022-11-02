#include "potential.h"

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/memory.h"
#include "../module_xc/xc_functional.h"
#include "global.h"
#include "math.h"
// new
#include "module_elecstate/potentials/efield.h"
#include "../module_surchem/surchem.h"
#include "module_elecstate/potentials/gatefield.h"
#include "H_Hartree_pw.h"
#ifdef __LCAO
#include "../src_lcao/ELEC_evolve.h"
#endif
#include "../module_base/timer.h"

Potential::Potential()
{
    vltot = nullptr;
    vr_eff1 = nullptr;
    vofk_eff1 = nullptr;
    this->out_pot = 0;
}

Potential::~Potential()
{
    delete[] vltot;
    delete[] vr_eff1;
    delete[] vofk_eff1;
#ifdef __CUDA
    cudaFree(d_vr_eff1);
#endif
}

void Potential::allocate(const int nrxx)
{
    ModuleBase::TITLE("Potential", "allocate");
    assert(nrxx >= 0);

    delete[] this->vltot;
    if (nrxx > 0)
        this->vltot = new double[nrxx];
    else
        this->vltot = nullptr;
    ModuleBase::Memory::record("Potential", "vltot", nrxx, "double");

    this->vr.create(GlobalV::NSPIN, nrxx);
    this->vr_eff.create(GlobalV::NSPIN, nrxx);
    ModuleBase::Memory::record("Potential", "vr", GlobalV::NSPIN * nrxx, "double");
    ModuleBase::Memory::record("Potential", "vr_eff", GlobalV::NSPIN * nrxx, "double");

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        this->vofk.create(GlobalV::NSPIN, nrxx);
        ModuleBase::Memory::record("Potential", "vofk", GlobalV::NSPIN * nrxx, "double");
        delete[] this->vofk_eff1;
        if(nrxx > 0)    this->vofk_eff1 = new double[nrxx];
        else            this->vofk_eff1 = nullptr;
        ModuleBase::Memory::record("Potential", "vofk_eff1", nrxx, "double");   
    }

    delete[] this->vr_eff1;
    if (nrxx > 0)
        this->vr_eff1 = new double[nrxx];
    else
        this->vr_eff1 = nullptr;
#ifdef __CUDA
    cudaMalloc((void **)&this->d_vr_eff1, nrxx * sizeof(double));
#endif
    ModuleBase::Memory::record("Potential", "vr_eff1", nrxx, "double");

    this->vnew.create(GlobalV::NSPIN, nrxx);
    ModuleBase::Memory::record("Potential", "vnew", GlobalV::NSPIN * nrxx, "double");

    if (GlobalV::imp_sol)
    {
        GlobalC::solvent_model.allocate(nrxx, GlobalV::NSPIN);
    }

    return;
}

//----------------------------------------------------------
//  Initializes the self consistent potential
//----------------------------------------------------------
void Potential::init_pot(const int &istep, // number of ionic steps
                         ModuleBase::ComplexMatrix &sf // structure factors
)
{
    ModuleBase::TITLE("Potential", "init_pot");
    ModuleBase::timer::tick("Potential", "init_pot");

    assert(istep >= 0);

    // total potential in real space
    this->vr_eff.zero_out();

    // the vltot should and must be zero here.
    ModuleBase::GlobalFunc::ZEROS(this->vltot, GlobalC::rhopw->nrxx);

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        this->vofk.zero_out();
    }

    //-------------------------------------------------------------------
    // (1) local pseudopotential + electric field (if any) in vltot
    //-------------------------------------------------------------------
    if (GlobalV::VION_IN_H)
    {
        this->set_local_pot(this->vltot, // 3D local pseudopotentials
                            GlobalC::ucell.ntype,
                            GlobalC::ppcell.vloc,
                            GlobalC::rhopw,
                            sf // structure factors
        );
    }
    else
    {
        for (int ir = 0; ir < GlobalC::rhopw->nrxx; ++ir)
        {
            this->vltot[ir] = 0.0;
        }
    }

    // zhengdy-soc, pauli matrix, just index 0 has vlocal term
    int nspin0 = GlobalV::NSPIN;

    if (GlobalV::NSPIN == 4)
    {
        nspin0 = 1;
    }

    for (int is = 0; is < nspin0; ++is)
    {
        for (int ir = 0; ir < GlobalC::rhopw->nrxx; ++ir)
        {
            this->vr_eff(is, ir) = this->vltot[ir];
        }
    }

    // core correction potential.
    GlobalC::CHR.set_rho_core(GlobalC::sf.strucFac);

    //--------------------------------------------------------------------
    // (2) other effective potentials need charge density,
    // choose charge density from ionic step 0.
    //--------------------------------------------------------------------
    if (istep == 0)
    {
        GlobalC::CHR.init_rho();
    }

    // renormalize the charge density
    GlobalC::CHR.renormalize_rho();

    //----------------------------------------------------------
    // (3) compute Hartree and XC potentials saves in vr
    //----------------------------------------------------------
    this->vr = this->v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);

    //----------------------------------------------------------
    // (4) total potentials
    //----------------------------------------------------------
#ifdef __LCAO
    if (ELEC_evolve::td_vext == 0)
    {
        this->set_vr_eff();
    }
    else
    {
        this->set_vrs_tddft(istep);
    }
#else
    this->set_vr_eff();
#endif

    // plots
    // figure::picture(this->vr_eff1,GlobalC::rhopw->nx,GlobalC::rhopw->ny,GlobalC::rhopw->nz);
    ModuleBase::timer::tick("Potential", "init_pot");
    return;
}

//==========================================================
// This routine computes the local potential in real space
//==========================================================
void Potential::set_local_pot(double *vl_pseudo, // store the local pseudopotential
                              const int &ntype, // number of atom types
                              ModuleBase::matrix &vloc, // local pseduopotentials
                              ModulePW::PW_Basis *rho_basis,
                              ModuleBase::ComplexMatrix &sf // structure factors
) const
{
    ModuleBase::TITLE("Potential", "set_local_pot");
    ModuleBase::timer::tick("Potential", "set_local_pot");

    std::complex<double> *vg = new std::complex<double>[rho_basis->npw];

    ModuleBase::GlobalFunc::ZEROS(vg, rho_basis->npw);

    for (int it = 0; it < ntype; it++)
    {
        for (int ig = 0; ig < rho_basis->npw; ig++)
        {
            vg[ig] += vloc(it, rho_basis->ig2igg[ig]) * sf(it, ig);
        }
    }

    GlobalC::UFFT.ToRealSpace(vg, vl_pseudo, rho_basis);

    if (GlobalV::EFIELD_FLAG && !GlobalV::DIP_COR_FLAG)
    {
        ModuleBase::matrix v_efield(GlobalV::NSPIN, GlobalC::rhopw->nrxx);
        v_efield = elecstate::Efield::add_efield(GlobalC::ucell, GlobalC::rhopw, GlobalV::NSPIN, GlobalC::CHR.rho, GlobalC::solvent_model);
        for (int ir = 0; ir < GlobalC::rhopw->nrxx; ++ir)
        {
            vl_pseudo[ir] += v_efield(0, ir);
        }
    }

    if( GlobalV::GATE_FLAG)
    {
        elecstate::Gatefield::add_gatefield(vl_pseudo, GlobalC::ucell, GlobalC::rhopw, true, true);
    }

    delete[] vg;

    // GlobalV::ofs_running <<" set local pseudopotential done." << std::endl;
    ModuleBase::timer::tick("Potential", "set_local_pot");
    return;
}

//==========================================================
// This routine computes the Hartree and Exchange and Correlation
// potential and energies which corresponds to a given charge density
// The XC potential is computed in real space, while the
// Hartree potential is computed in reciprocal space.
//==========================================================
ModuleBase::matrix Potential::v_of_rho(const double *const *const rho_in, const double *const rho_core_in)
{
    ModuleBase::TITLE("Potential", "v_of_rho");
    ModuleBase::timer::tick("Potential", "v_of_rho");

    ModuleBase::matrix v(GlobalV::NSPIN, GlobalC::rhopw->nrxx);

    //----------------------------------------------------------
    //  calculate the exchange-correlation potential
    //----------------------------------------------------------

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
#ifdef USE_LIBXC
        const std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> etxc_vtxc_v
            = XC_Functional::v_xc_meta(GlobalC::rhopw->nrxx,
                                       GlobalC::rhopw->nxyz,
                                       GlobalC::ucell.omega,
                                       rho_in,
                                       GlobalC::CHR.rho_core,
                                       GlobalC::CHR.kin_r);
        GlobalC::en.etxc = std::get<0>(etxc_vtxc_v);
        GlobalC::en.vtxc = std::get<1>(etxc_vtxc_v);
        v += std::get<2>(etxc_vtxc_v);
        vofk = std::get<3>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("v_of_rho", "to use mGGA, compile with LIBXC");
#endif
    }
    else
    {
        const std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v = XC_Functional::v_xc(GlobalC::rhopw->nrxx,
                                                                                               GlobalC::rhopw->nxyz,
                                                                                               GlobalC::ucell.omega,
                                                                                               rho_in,
                                                                                               GlobalC::CHR.rho_core);
        GlobalC::en.etxc = std::get<0>(etxc_vtxc_v);
        GlobalC::en.vtxc = std::get<1>(etxc_vtxc_v);
        v += std::get<2>(etxc_vtxc_v);
    }

    //----------------------------------------------------------
    //  calculate the Hartree potential
    //----------------------------------------------------------
    if (GlobalV::VH_IN_H)
    {
        v += H_Hartree_pw::v_hartree(GlobalC::ucell, GlobalC::rhopw, GlobalV::NSPIN, rho_in);
        if (GlobalV::imp_sol)
        {
            v += GlobalC::solvent_model.v_correction(GlobalC::ucell, GlobalC::rhopw, GlobalV::NSPIN, rho_in);
        }
    }

    //----------------------------------------------------------
    //  calculate the efield and dipole correction
    //----------------------------------------------------------
    if (GlobalV::EFIELD_FLAG && GlobalV::DIP_COR_FLAG)
    {
        v += elecstate::Efield::add_efield(GlobalC::ucell, GlobalC::rhopw, GlobalV::NSPIN, rho_in, GlobalC::solvent_model);
    }

    // test get ntot_reci
    // complex<double> *tmpn = new complex<double>[GlobalC::rhopw->npw];
    // ModuleBase::GlobalFunc::ZEROS(tmpn, GlobalC::rhopw->npw);
    // GlobalC::solvent_model.get_totn_reci(GlobalC::ucell, GlobalC::rhopw, tmpn);
    // delete[] tmpn;

    ModuleBase::timer::tick("Potential", "v_of_rho");
    return v;
} // end subroutine v_of_rho

//==========================================================
// set the effective potential vr_eff on the real space grid
// used in h_psi, adding the (spin dependent) scf (H+xc)
// part and the sum of all the local pseudopotential
// contributions.
//==========================================================
void Potential::set_vr_eff(void)
{
    ModuleBase::TITLE("Potential", "set_vr_eff");
    ModuleBase::timer::tick("Potential", "set_vr_eff");

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        //=================================================================
        // define the total local potential (external + scf) for each spin
        //=================================================================
        if (GlobalV::NSPIN == 4 && is > 0)
        {
            for (int i = 0; i < GlobalC::rhopw->nrxx; i++)
            {
                this->vr_eff(is, i) = this->vr(is, i);
            }
        }
        else
        {
            for (int i = 0; i < GlobalC::rhopw->nrxx; i++)
            {
                this->vr_eff(is, i) = this->vltot[i] + this->vr(is, i);
            }
        }
    }

    ModuleBase::timer::tick("Potential", "set_vr_eff");
    return;
}
