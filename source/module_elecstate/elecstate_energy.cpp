#include <cmath>

#include "elecstate.h"
#include "elecstate_getters.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"
#ifdef USE_PAW
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#endif

namespace elecstate
{
/// @brief calculate band gap
void ElecState::cal_bandgap()
{
    if (this->ekb.nr == 0 || this->ekb.nc == 0)
    { // which means no homo and no lumo
        this->bandgap = 0.0;
        return;
    }
    int nbands = PARAM.inp.nbands;
    int nks = this->klist->get_nks();
    double homo = this->ekb(0, 0);
    double lumo = this->ekb(0, nbands - 1);
    for (int ib = 0; ib < nbands; ib++)
    {
        for (int ik = 0; ik < nks; ik++)
        {
            if (!(this->ekb(ik, ib) - this->eferm.ef > 1e-5) && homo < this->ekb(ik, ib))
            {
                homo = this->ekb(ik, ib);
            }
            if (this->ekb(ik, ib) - this->eferm.ef > 1e-5 && lumo > this->ekb(ik, ib))
            {
                lumo = this->ekb(ik, ib);
            }
        }
    }
    this->bandgap = lumo - homo;
}

/// @brief calculate spin up & down band gap
/// @todo add isk[ik] so as to discriminate different spins
void ElecState::cal_bandgap_updw()
{
    if (this->ekb.nr == 0 || this->ekb.nc == 0)
    { // which means no homo and no lumo
        this->bandgap_up = 0.0;
        this->bandgap_dw = 0.0;
        return;
    }
    int nbands = PARAM.inp.nbands;
    int nks = this->klist->get_nks();
    double homo_up = this->ekb(0, 0);
    double lumo_up = this->ekb(0, nbands - 1);
    double homo_dw = this->ekb(0, 0);
    double lumo_dw = this->ekb(0, nbands - 1);
    for (int ib = 0; ib < nbands; ib++)
    {
        for (int ik = 0; ik < nks; ik++)
        {
            if (this->klist->isk[ik] == 0)
            {
                if (!(this->ekb(ik, ib) - this->eferm.ef_up > 1e-5) && homo_up < this->ekb(ik, ib))
                {
                    homo_up = this->ekb(ik, ib);
                }
                if (this->ekb(ik, ib) - this->eferm.ef_up > 1e-5 && lumo_up > this->ekb(ik, ib))
                {
                    lumo_up = this->ekb(ik, ib);
                }
            }
            if (this->klist->isk[ik] == 1)
            {
                if (!(this->ekb(ik, ib) - this->eferm.ef_dw > 1e-5) && homo_dw < this->ekb(ik, ib))
                {
                    homo_dw = this->ekb(ik, ib);
                }
                if (this->ekb(ik, ib) - this->eferm.ef_dw > 1e-5 && lumo_dw > this->ekb(ik, ib))
                {
                    lumo_dw = this->ekb(ik, ib);
                }
            }
        }
    }
    this->bandgap_up = lumo_up - homo_up;
    this->bandgap_dw = lumo_dw - homo_dw;
}

/// @brief calculate deband
double ElecState::cal_delta_eband(const UnitCell& ucell) const
{
    // out potentials from potential mixing
    // total energy and band energy corrections
    double deband0 = 0.0;

    double deband_aux = 0.0;

    // only potential related with charge is used here for energy correction
    // on the fly calculate it here by v_effective - v_fixed
    const double* v_eff = this->pot->get_effective_v(0);
    const double* v_fixed = this->pot->get_fixed_v();
    const double* v_ofk = nullptr;
    const bool v_ofk_flag =(get_xc_func_type() == 3 
                            || get_xc_func_type() == 5);
#ifdef USE_PAW
    if(PARAM.inp.use_paw)
    {
        ModuleBase::matrix v_xc;
        const std::tuple<double, double, ModuleBase::matrix> etxc_vtxc_v
            = XC_Functional::v_xc(this->charge->nrxx, this->charge, &ucell);
        v_xc = std::get<2>(etxc_vtxc_v);

        for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
        {
            deband_aux -= this->charge->rho[0][ir] * v_xc(0,ir);
        }
        if (PARAM.inp.nspin == 2)
        {
            for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
            {
                deband_aux -= this->charge->rho[1][ir] * v_xc(1,ir);
            }
        }
    }
#endif

    if(!PARAM.inp.use_paw)
    {
        for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
        {
            deband_aux -= this->charge->rho[0][ir] * (v_eff[ir] - v_fixed[ir]);
        }
        if (v_ofk_flag)
        {
            v_ofk = this->pot->get_effective_vofk(0);
            // cause in the get_effective_vofk, the func will return nullptr
            if(v_ofk==nullptr && this->charge->rhopw->nrxx>0)
            {
                ModuleBase::WARNING_QUIT("ElecState::cal_delta_eband","v_ofk is nullptr");
            }
            for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
            {
                deband_aux -= this->charge->kin_r[0][ir] * v_ofk[ir];
            }
        }
        
        if (PARAM.inp.nspin == 2)
        {
            v_eff = this->pot->get_effective_v(1);
            for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
            {
                deband_aux -= this->charge->rho[1][ir] * (v_eff[ir] - v_fixed[ir]);
            }
            if (v_ofk_flag)
            {
                v_ofk = this->pot->get_effective_vofk(1);
                if(v_ofk==nullptr && this->charge->rhopw->nrxx>0)
                {
                    ModuleBase::WARNING_QUIT("ElecState::cal_delta_eband","v_ofk is nullptr");
                }
                for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
                {
                    deband_aux -= this->charge->kin_r[1][ir] * v_ofk[ir];
                }
            }
        }
        else if (PARAM.inp.nspin == 4)
        {
            for (int is = 1; is < 4; is++)
            {
                v_eff = this->pot->get_effective_v(is);
                for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
                {
                    deband_aux -= this->charge->rho[is][ir] * v_eff[ir];
                }
            }
        }
    }

#ifdef __MPI
    MPI_Allreduce(&deband_aux, &deband0, 1, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#else
    deband0 = deband_aux;
#endif

    deband0 *= this->omega / this->charge->rhopw->nxyz;

    // \int rho(r) v_{exx}(r) dr = 2 E_{exx}[rho]
    deband0 -= 2 * this->f_en.exx; // Peize Lin add 2017-10-16
    return deband0;
}

/// @brief calculate descf
double ElecState::cal_delta_escf() const
{
    ModuleBase::TITLE("energy", "delta_escf");
    double descf = 0.0;

    // now rho1 is "mixed" charge density
    // and rho1_save is "output" charge density
    // because in "deband" the energy is calculated from "output" charge density,
    // so here is the correction.
    // only potential related with charge is used here for energy correction
    // on the fly calculate it here by v_effective - v_fixed
    const double* v_eff = this->pot->get_effective_v(0);
    const double* v_fixed = this->pot->get_fixed_v();
    const double* v_ofk = nullptr;

    if (get_xc_func_type() == 3 || get_xc_func_type() == 5)
    {
        v_ofk = this->pot->get_effective_vofk(0);
    }
    for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
    {
        descf -= (this->charge->rho[0][ir] - this->charge->rho_save[0][ir]) * (v_eff[ir] - v_fixed[ir]);
        if (get_xc_func_type() == 3 || get_xc_func_type() == 5)
        {
            // cause in the get_effective_vofk, the func will return nullptr
            assert(v_ofk!=nullptr);
            descf -= (this->charge->kin_r[0][ir] - this->charge->kin_r_save[0][ir]) * v_ofk[ir];
        }
    }

    if (PARAM.inp.nspin == 2)
    {
        v_eff = this->pot->get_effective_v(1);
        if (get_xc_func_type() == 3 || get_xc_func_type() == 5)
        {
            v_ofk = this->pot->get_effective_vofk(1);
        }
        for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
        {
            descf -= (this->charge->rho[1][ir] - this->charge->rho_save[1][ir]) * (v_eff[ir] - v_fixed[ir]);
            if (get_xc_func_type() == 3 || get_xc_func_type() == 5)
            {
                descf -= (this->charge->kin_r[1][ir] - this->charge->kin_r_save[1][ir]) * v_ofk[ir];
            }
        }
    }
    if (PARAM.inp.nspin == 4)
    {
        for (int is = 1; is < 4; is++)
        {
            v_eff = this->pot->get_effective_v(is);
            for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
            {
                descf -= (this->charge->rho[is][ir] - this->charge->rho_save[is][ir]) * v_eff[ir];
            }
        }
    }

#ifdef __MPI
    Parallel_Reduce::reduce_pool(descf);
#endif

    assert(this->charge->rhopw->nxyz > 0);

    descf *= this->omega / this->charge->rhopw->nxyz;
    return descf;
}

/// @brief calculation if converged
void ElecState::cal_converged()
{
    // update etxc and vtxc
    // allocate vnew in get_vnew()
    this->pot->get_vnew(this->charge, this->vnew);
    this->vnew_exist = true;
    // vnew will be used in force_scc()

    // set descf to 0
    this->f_en.descf = 0.0;
}

/**
 * @brief calculate energies
 *
 * @param type: 1 means Harris-Foulkes functinoal; 
 * @param type: 2 means Kohn-Sham functional;
 */
void ElecState::cal_energies(const int type)
{
    //! Hartree energy
    this->f_en.hartree_energy = get_hartree_energy();

    //! energy from E-field
    this->f_en.efield = get_etot_efield();

    //! energy from gate-field
    this->f_en.gatefield = get_etot_gatefield();

    //! energy from implicit solvation model 
    if (PARAM.inp.imp_sol)
    {
        this->f_en.esol_el = get_solvent_model_Ael();
        this->f_en.esol_cav = get_solvent_model_Acav();
    }

    //! spin constrained energy
    if (PARAM.inp.sc_mag_switch)
    {
        this->f_en.escon = get_spin_constrain_energy();
    }

     // energy from DFT+U
    if (PARAM.inp.dft_plus_u)
    {
        this->f_en.edftu = get_dftu_energy();
    }

#ifdef __DEEPKS
    // energy from deepks
    if (PARAM.inp.deepks_scf)
    {
        this->f_en.edeepks_scf = get_deepks_E_delta() - get_deepks_E_delta_band();
    }
#endif

    if (type == 1) // Harris-Foulkes functional
    {
        this->f_en.calculate_harris();
    }
    else if (type == 2)// Kohn-Sham functional
    {
        this->f_en.calculate_etot();
    }
    else
    {
        ModuleBase::WARNING_QUIT("cal_energies", "The form of total energy functional is unknown!");
    }
}

} // namespace elecstate
