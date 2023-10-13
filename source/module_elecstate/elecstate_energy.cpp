#include <cmath>

#include "elecstate.h"
#include "elecstate_getters.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"

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
    int nbands = GlobalV::NBANDS;
    int nks = this->klist->nks;
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
    int nbands = GlobalV::NBANDS;
    int nks = this->klist->nks;
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
double ElecState::cal_delta_eband() const
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
    if (get_xc_func_type() == 3 || get_xc_func_type() == 5)
    {
        v_ofk = this->pot->get_effective_vofk(0);
    }

    for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
    {
        deband_aux -= this->charge->rho[0][ir] * (v_eff[ir] - v_fixed[ir]);
        if (get_xc_func_type() == 3 || get_xc_func_type() == 5)
        {
            deband_aux -= this->charge->kin_r[0][ir] * v_ofk[ir];
        }
    }

    if (GlobalV::NSPIN == 2)
    {
        v_eff = this->pot->get_effective_v(1);
        v_ofk = this->pot->get_effective_vofk(1);
        for (int ir = 0; ir < this->charge->rhopw->nrxx; ir++)
        {
            deband_aux -= this->charge->rho[1][ir] * (v_eff[ir] - v_fixed[ir]);
            if (get_xc_func_type() == 3 || get_xc_func_type() == 5)
            {
                deband_aux -= this->charge->kin_r[1][ir] * v_ofk[ir];
            }
        }
    }
    else if (GlobalV::NSPIN == 4)
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
            descf -= (this->charge->kin_r[0][ir] - this->charge->kin_r_save[0][ir]) * v_ofk[ir];
        }
    }

    if (GlobalV::NSPIN == 2)
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
    if (GlobalV::NSPIN == 4)
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
 * @param type: 1 means harris energy; 2 means etot
 */
void ElecState::cal_energies(const int type)
{
    this->f_en.hartree_energy = get_hartree_energy();
    this->f_en.efield = get_etot_efield();
    this->f_en.gatefield = get_etot_gatefield();
    if (GlobalV::imp_sol)
    {
        this->f_en.esol_el = get_solvent_model_Ael();
        this->f_en.esol_cav = get_solvent_model_Acav();
    }
#ifdef __LCAO
    if (GlobalV::sc_mag_switch)
    {
        this->f_en.escon = get_spin_constrain_energy();
    }
    if (GlobalV::dft_plus_u)
    {
        this->f_en.edftu = get_dftu_energy();
    }
#endif
#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        this->f_en.edeepks_scf = get_deepks_E_delta() - get_deepks_E_delta_band();
    }
#endif
    if (type == 1) // harris
    {
        this->f_en.calculate_harris();
    }
    else // etot
    {
        this->f_en.calculate_etot();
    }
}

} // namespace elecstate