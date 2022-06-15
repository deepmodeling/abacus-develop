#include "potential.h"

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/memory.h"
#include "../module_xc/xc_functional.h"
#include "global.h"
#include "math.h"
// new
#include "../module_surchem/surchem.h"
#include "H_Hartree_pw.h"
#include "../module_surchem/efield.h"
#ifdef __LCAO
#include "../src_lcao/ELEC_evolve.h"
#endif
#include "../module_base/timer.h"

Potential::Potential()
{
    vltot = nullptr;
    vr_eff1 = nullptr;
    this->out_pot = 0;
}

Potential::~Potential()
{
    delete[] vltot;
    delete[] vr_eff1;
#ifdef __CUDA
    cudaFree(d_vr_eff1);
#endif
}

void Potential::allocate(const int nrxx)
{
    ModuleBase::TITLE("Potential", "allocate");
    assert(nrxx >= 0);

    delete[] this->vltot;
    if(nrxx > 0)    this->vltot = new double[nrxx];
    else            this->vltot = nullptr;
    ModuleBase::Memory::record("Potential", "vltot", nrxx, "double");

    this->vr.create(GlobalV::NSPIN, nrxx);
    this->vr_eff.create(GlobalV::NSPIN, nrxx);
    ModuleBase::Memory::record("Potential", "vr", GlobalV::NSPIN * nrxx, "double");
    ModuleBase::Memory::record("Potential", "vr_eff", GlobalV::NSPIN * nrxx, "double");

    if (XC_Functional::get_func_type() == 3)
    {
        this->vofk.create(GlobalV::NSPIN, nrxx);
        ModuleBase::Memory::record("Potential", "vofk", GlobalV::NSPIN * nrxx, "double");
    }

    delete[] this->vr_eff1;
    if(nrxx > 0)    this->vr_eff1 = new double[nrxx];
    else            this->vr_eff1 = nullptr;
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

    if (XC_Functional::get_func_type() == 3)
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
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "init_chg", init_chg);

        std::cout << " START POTENTIAL      : " << init_chg << std::endl;
        if (this->init_chg == "atomic") // mohan add 2007-10-17
        {
        start_from_atomic:
            GlobalC::CHR.atomic_rho(GlobalV::NSPIN, GlobalC::CHR.rho, GlobalC::rhopw);
        }
        else if (this->init_chg == "file")
        {
            GlobalV::ofs_running << " try to start potential from file : ";
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                std::stringstream ssc;
                ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_CHG";
                GlobalV::ofs_running << ssc.str() << std::endl;
                // mohan update 2012-02-10
                if (GlobalC::CHR.read_rho(is, ssc.str(), GlobalC::CHR.rho[is]))
                {
                    GlobalV::ofs_running << " Read in the charge density: " << ssc.str() << std::endl;
                }
                else if (is > 0 && GlobalV::NSPIN == 4)
                {
                    // read only spin (up+down)
                    if (GlobalV::PRENSPIN == 1)
                    {
                        GlobalV::ofs_running << " Didn't read in the charge density but autoset it for spin " << is + 1
                                             << std::endl;
                        for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
                        {
                            GlobalC::CHR.rho[is][ir] = 0.0;
                        }
                    }
                    //
                    else if (GlobalV::PRENSPIN == 2)
                    { // read up and down , then rearrange them.
                        if (is == 1)
                        {
                            ModuleBase::WARNING_QUIT("potential::init_pot", "Incomplete charge density file!");
                        }
                        else if (is == 2)
                        {
                            GlobalV::ofs_running << " Didn't read in the charge density but would rearrange it later. "
                                                 << std::endl;
                        }
                        else if (is == 3)
                        {
                            GlobalV::ofs_running << " rearrange charge density " << std::endl;
                            for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
                            {
                                GlobalC::CHR.rho[3][ir] = GlobalC::CHR.rho[0][ir] - GlobalC::CHR.rho[1][ir];
                                GlobalC::CHR.rho[0][ir] = GlobalC::CHR.rho[0][ir] + GlobalC::CHR.rho[1][ir];
                                GlobalC::CHR.rho[1][ir] = 0.0;
                                GlobalC::CHR.rho[2][ir] = 0.0;
                            }
                        }
                    }
                    else
                    {
                        ModuleBase::WARNING_QUIT("potential::init_pot", "Incomplete charge density file!");
                    }
                }
                else
                {
                    GlobalV::ofs_running << " Start charge density from atomic charge density." << std::endl;
                    goto start_from_atomic;
                }
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("potential::init_pot", "init_chg is wrong!");
        }

        // Peize Lin add 2020.04.04
        if (GlobalC::restart.info_load.load_charge && !GlobalC::restart.info_load.load_charge_finish)
        {
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                GlobalC::restart.load_disk("charge", is);
            }
            GlobalC::restart.info_load.load_charge_finish = true;
        }
    }
    else
    {
        // the extrapolation part moves to ions.cpp.
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
                              ModulePW::PW_Basis* rho_basis,
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
        v_efield = Efield::add_efield(GlobalC::ucell, GlobalC::rhopw, GlobalV::NSPIN, GlobalC::CHR.rho);
        for (int ir = 0; ir < GlobalC::rhopw->nrxx; ++ir)
        {
            vl_pseudo[ir] += v_efield(0, ir);
        }
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

    if (XC_Functional::get_func_type() == 3)
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

            // test energy outside
            cout << "energy Outside: " << endl;
            GlobalC::solvent_model.cal_Ael(GlobalC::ucell, GlobalC::rhopw);
            GlobalC::solvent_model.cal_Acav(GlobalC::ucell, GlobalC::rhopw);
        }
    }

    //----------------------------------------------------------
    //  calculate the efield and dipole correction
    //----------------------------------------------------------
    if (GlobalV::EFIELD_FLAG && GlobalV::DIP_COR_FLAG)
    {
        v += Efield::add_efield(GlobalC::ucell, GlobalC::rhopw, GlobalV::NSPIN, rho_in);
    }


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

// ----------------------------------------------------------------------
void Potential::newd(void)
{
    ModuleBase::TITLE("Potential", "newd");

    // distringuish non-local pseudopotential in REAL or RECIPROCAL space.
    // if in real space, call new_r
    // if in reciprocal space, call new_g

    // new g:
    //----------------------------------------------------------------------
    //  This routine computes the integral of the effective potential with
    //  the Q function and adds it to the bare ionic D term which is used
    //  to compute the non-local term in the US scheme.

    // no ultrasoft potentials: use bare coefficients for projectors
    // if( spin_orbital) ....
    // else if(noncolin) ....
    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {
        const int it = GlobalC::ucell.iat2it[iat];
        const int nht = GlobalC::ucell.atoms[it].nh;
        // nht: number of beta functions per atom type
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            for (int ih = 0; ih < nht; ih++)
            {
                for (int jh = ih; jh < nht; jh++)
                {
                    if (GlobalV::LSPINORB)
                    {
                        GlobalC::ppcell.deeq_nc(is, iat, ih, jh) = GlobalC::ppcell.dvan_so(is, it, ih, jh);
                        GlobalC::ppcell.deeq_nc(is, iat, jh, ih) = GlobalC::ppcell.dvan_so(is, it, jh, ih);
                    }
                    else if (GlobalV::NSPIN == 4)
                    {
                        if (is == 0)
                        {
                            GlobalC::ppcell.deeq_nc(is, iat, ih, jh) = GlobalC::ppcell.dvan(it, ih, jh);
                            GlobalC::ppcell.deeq_nc(is, iat, jh, ih) = GlobalC::ppcell.dvan(it, ih, jh);
                        }
                        else if (is == 1)
                        {
                            GlobalC::ppcell.deeq_nc(is, iat, ih, jh) = std::complex<double>(0.0, 0.0);
                            GlobalC::ppcell.deeq_nc(is, iat, jh, ih) = std::complex<double>(0.0, 0.0);
                        }
                        else if (is == 2)
                        {
                            GlobalC::ppcell.deeq_nc(is, iat, ih, jh) = std::complex<double>(0.0, 0.0);
                            GlobalC::ppcell.deeq_nc(is, iat, jh, ih) = std::complex<double>(0.0, 0.0);
                        }
                        else if (is == 3)
                        {
                            GlobalC::ppcell.deeq_nc(is, iat, ih, jh) = GlobalC::ppcell.dvan(it, ih, jh);
                            GlobalC::ppcell.deeq_nc(is, iat, jh, ih) = GlobalC::ppcell.dvan(it, ih, jh);
                        }
                    }
                    else
                    {
                        GlobalC::ppcell.deeq(is, iat, ih, jh) = GlobalC::ppcell.dvan(it, ih, jh);
                        GlobalC::ppcell.deeq(is, iat, jh, ih) = GlobalC::ppcell.dvan(it, ih, jh);
                    }
                }
            }
        }
    }
#ifdef __CUDA
    cudaMemcpy(GlobalC::ppcell.d_deeq,
               GlobalC::ppcell.deeq.ptr,
               GlobalV::NSPIN * GlobalC::ucell.nat * GlobalC::ppcell.nhm * GlobalC::ppcell.nhm * sizeof(double),
               cudaMemcpyHostToDevice);
#endif
    return;
} // end subroutine newd
