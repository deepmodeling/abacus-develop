#include "esolver_of_tddft.h"

#include "source_io/module_parameter/parameter.h"
#include "source_io/cube_io.h"
#include "source_io/output_log.h"
#include "source_io/write_elecstat_pot.h"
//-----------temporary-------------------------
#include "source_base/global_function.h"
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_hamilt/module_ewald/H_Ewald_pw.h"
#include "source_pw/module_pwdft/global.h"
#include "source_io/print_info.h"
#include "source_estate/cal_ux.h"
//-----force-------------------
#include "source_pw/module_pwdft/forces.h"
//-----stress------------------
#include "source_pw/module_ofdft/of_stress_pw.h"

namespace ModuleESolver
{

ESolver_OF_TDDFT::ESolver_OF_TDDFT()
{
    this->classname = "ESolver_OF_TDDFT";
    /*this->pphi_td= new std::complex<double>*[PARAM.inp.nspin];

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        this->pphi_td[is]= new std::complex<double>[pw_rho->nrxx];
    }*/
}

ESolver_OF_TDDFT::~ESolver_OF_TDDFT()
{
    //delete psi_td;
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] this->pphi_td[is];
    }
    delete[] this->pphi_td;
}


void ESolver_OF_TDDFT::runner(UnitCell& ucell, const int istep)
{
    ModuleBase::timer::tick("ESolver_OF_TDDFT", "runner");
    // get Ewald energy, initial rho and phi if necessary
    this->before_opt(istep, ucell);
    this->iter_ = 0;

    bool conv_esolver = false; // this conv_esolver is added by mohan 20250302 
#ifdef __MPI
    this->iter_time = MPI_Wtime();
#else
    this->iter_time = std::chrono::system_clock::now();
#endif

    if (istep==0)
    {
        this->pphi_td= new std::complex<double>*[PARAM.inp.nspin];

        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            this->pphi_td[is]= new std::complex<double>[pw_rho->nrxx];
        }
    }

    if ((istep<1) && PARAM.inp.init_chg != "file")
    {
        while (true)
        {
            // once we get a new rho and phi, update potential
            this->update_potential(ucell);

            // calculate the energy of new rho and phi
            this->energy_llast_ = this->energy_last_;
            this->energy_last_ = this->energy_current_;
            this->energy_current_ = this->cal_energy();


            // check if the job is done
            if (this->check_exit(conv_esolver))
            {
                break;
            }

            // find the optimization direction and step lenghth theta according to the potential
            this->optimize(ucell);

            std::cout<<"optimize------"<<std::endl;

            // update the rho and phi based on the direction and theta
            this->update_rho();

            this->iter_++;

            ESolver_FP::iter_finish(ucell, istep, this->iter_, conv_esolver);
        }

        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ir = 0; ir < pw_rho->nrxx; ++ir)
            {
                pphi_td[is][ir]=pphi_[is][ir];
            }
        }
    }
    else
    {
         std::cout<<"propagate_psi -------"<<std::endl;
        this->propagate_psi(ucell,this->pphi_td,this->pw_rho);
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ir = 0; ir < pw_rho->nrxx; ++ir)
            {
                pphi_[is][ir]=abs(pphi_td[is][ir]);
            }
        }
        conv_esolver=true;
    }

    this->after_opt(istep, ucell, conv_esolver);

    ModuleBase::timer::tick("ESolver_OF_TDDFT", "runner");
}

/**
 * @brief Prepare to optimize the charge density,
 * update elecstate, kedf, and opts if needed
 * calculate ewald energy, initialize the rho, phi, theta
 *
 * @param istep
 * @param ucell
 */
void ESolver_OF_TDDFT::before_opt(const int istep, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_OF", "before_opt");
    ModuleBase::timer::tick("ESolver_OF", "before_opt");

    //! 1) call before_scf() of ESolver_FP
    ESolver_FP::before_scf(ucell, istep);

    if (ucell.cell_parameter_updated)
    {
        this->dV_ = ucell.omega / this->pw_rho->nxyz;

        // initialize elecstate, including potential
        this->init_elecstate(ucell);

        // Initialize KEDF
        this->kedf_manager_->init(PARAM.inp, this->pw_rho, this->dV_, this->nelec_[0]);

        // Initialize optimization methods
        this->init_opt();

        // Refresh the arrays
        delete this->psi_;
        delete this->psi_td;
        this->psi_ = new psi::Psi<double>(1, 
                                          PARAM.inp.nspin, 
                                          this->pw_rho->nrxx,
                                          this->pw_rho->nrxx,
                                          true);
        /*this->psi_td = new psi::Psi<std::complex<double>>(1, 
                                          PARAM.inp.nspin, 
                                          this->pw_rho->nrxx,
                                          this->pw_rho->nrxx,
                                          true);*/
        //this->pphi_td= new std::complex<double>*[PARAM.inp.nspin];

        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            this->pphi_[is] = this->psi_->get_pointer(is);
            //this->pphi_td[is] = this->psi_td->get_pointer(is);
            //this->pphi_td[is]= new std::complex<double>[pw_rho->nrxx];
        }

        delete this->ptemp_rho_;
        this->ptemp_rho_ = new Charge();
        this->ptemp_rho_->set_rhopw(this->pw_rho);
        this->ptemp_rho_->allocate(PARAM.inp.nspin);

        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            delete[] this->pdLdphi_[is];
            delete[] this->pdEdphi_[is];
            delete[] this->pdirect_[is];
            delete[] this->precip_dir_[is];
            this->pdLdphi_[is] = new double[this->pw_rho->nrxx];
            this->pdEdphi_[is] = new double[this->pw_rho->nrxx];
            this->pdirect_[is] = new double[this->pw_rho->nrxx];
            this->precip_dir_[is] = new std::complex<double>[pw_rho->npw];
        }
    }

    this->pelec->init_scf(istep, ucell, Pgrid, sf.strucFac, locpp.numeric, ucell.symm);

    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, this->chr, this->pw_rho, ucell.symm);
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        if (PARAM.inp.init_chg != "file")
        {
            for (int ibs = 0; ibs < this->pw_rho->nrxx; ++ibs)
            {
                // Here we initialize rho to be uniform,
                // because the rho got by pot.init_pot -> Charge::atomic_rho may contain minus elements.
                this->chr.rho[is][ibs] = this->nelec_[is] / this->pelec->omega;
                this->pphi_[is][ibs] = sqrt(this->chr.rho[is][ibs]);                                                 
            }
        }
        else
        {
            for (int ibs = 0; ibs < this->pw_rho->nrxx; ++ibs)
            {
                this->pphi_[is][ibs] = sqrt(this->chr.rho[is][ibs]);
            }
        }
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        this->pelec->eferm.set_efval(is, 0);
        this->theta_[is] = 0.;
        ModuleBase::GlobalFunc::ZEROS(this->pdLdphi_[is], this->pw_rho->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdEdphi_[is], this->pw_rho->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdirect_[is], this->pw_rho->nrxx);
    }
    if (PARAM.inp.nspin == 1)
    {
        this->theta_[0] = 0.2;
    }

    ModuleBase::timer::tick("ESolver_OF", "before_opt");
}

void ESolver_OF_TDDFT::get_Hpsi(UnitCell& ucell, const std::complex<double>* const* psi_, ModulePW::PW_Basis* pw_rho, std::complex<double>** Hpsi)
{
    // update rho
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            this->chr.rho[is][ir] = abs(psi_[is][ir])*abs(psi_[is][ir]);
        }
    }

    this->pelec->pot->update_from_charge(&this->chr, &ucell); // Hartree + XC + external
    this->get_tf_potential(this->chr.rho,pw_rho ,this->pelec->pot->get_effective_v()); // TF potential
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        const double* vr_eff = this->pelec->pot->get_effective_v(is);
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            Hpsi[is][ir] = vr_eff[ir]*psi_[is][ir];
        }
    }
    this->get_vw_potential_phi(psi_, pw_rho, Hpsi);
}

void ESolver_OF_TDDFT::get_tf_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot)
{
    const double c_tf_
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2;
    if (PARAM.inp.nspin == 1)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rpot(0, ir) += 5.0 / 3.0 * c_tf_ * std::pow(prho[0][ir], 2. / 3.);
        }
    }
    else if (PARAM.inp.nspin == 2)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ir = 0; ir < pw_rho->nrxx; ++ir)
            {
                rpot(is, ir) += 5.0 / 3.0 * c_tf_ * std::pow(2. * prho[is][ir], 2. / 3.);
            }
        }
    }
}

void ESolver_OF_TDDFT::get_vw_potential_phi(const std::complex<double>* const* pphi, ModulePW::PW_Basis* pw_rho, std::complex<double>** Hpsi)
{
    std::complex<double>** rLapPhi = new std::complex<double>*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        rLapPhi[is] = new std::complex<double>[pw_rho->nrxx];
    }
    std::complex<double>** recipPhi = new std::complex<double>*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        recipPhi[is] = new std::complex<double>[pw_rho->npw];

        pw_rho->real2recip(pphi[is], recipPhi[is]);
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            recipPhi[is][ik] *= pw_rho->gg[ik] * pw_rho->tpiba2;
        }
        pw_rho->recip2real(recipPhi[is], rLapPhi[is]);
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            Hpsi[is][ik]+=rLapPhi[is][ik];
        }
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        delete[] recipPhi[is];
        delete[] rLapPhi[is];
    }
    delete[] recipPhi;
    delete[] rLapPhi;
}

void ESolver_OF_TDDFT::get_CD_potential(const std::complex<double>* const* psi_, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot)
{
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        //recipCurrent = new std::complex<double>[pw_rho->npw];
        //delete[] recipCurrent;
    }
}

void ESolver_OF_TDDFT::propagate_psi(UnitCell& ucell, std::complex<double>** pphi_, ModulePW::PW_Basis* pw_rho)
{
    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagte_psi");

    std::complex<double> imag(0.0,1.0);
    double dt=PARAM.inp.mdp.md_dt;
    std::complex<double>** K1 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** K2 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** K3 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** K4 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** psi1 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** psi2 = new std::complex<double>*[PARAM.inp.nspin];
    std::complex<double>** psi3 = new std::complex<double>*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        K1[is] = new std::complex<double>[pw_rho->nrxx];
        K2[is] = new std::complex<double>[pw_rho->nrxx];
        K3[is] = new std::complex<double>[pw_rho->nrxx];
        K4[is] = new std::complex<double>[pw_rho->nrxx];
        psi1[is] = new std::complex<double>[pw_rho->nrxx];
        psi2[is] = new std::complex<double>[pw_rho->nrxx];
        psi3[is] = new std::complex<double>[pw_rho->nrxx];
    }
    get_Hpsi(ucell,pphi_,pw_rho,K1);
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K1[is][ir]=-1.0*K1[is][ir]*dt*imag;
            psi1[is][ir]=pphi_[is][ir]+0.5*K1[is][ir];
        }
    }
    get_Hpsi(ucell,psi1,pw_rho,K2);
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K2[is][ir]=-1.0*K2[is][ir]*dt*imag;
            psi2[is][ir]=pphi_[is][ir]+0.5*K2[is][ir];
        }
    }
    get_Hpsi(ucell,psi2,pw_rho,K3);
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K3[is][ir]=-1.0*K3[is][ir]*dt*imag;
            psi3[is][ir]=pphi_[is][ir]+K3[is][ir];
        }
    }
    get_Hpsi(ucell,psi3,pw_rho,K4);
    for (int is = 0; is < PARAM.inp.nspin; ++is){
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            K4[is][ir]=-1.0*K4[is][ir]*dt*imag;
            pphi_[is][ir]+=1.0/6.0*(K1[is][ir]+2.0*K2[is][ir]+2.0*K3[is][ir]+K4[is][ir]);
        }
    }

    for (int is = 0; is < PARAM.inp.nspin; ++is) {
        delete[] K1[is];
        delete[] K2[is];
        delete[] K3[is];
        delete[] K4[is];
        delete[] psi1[is];
        delete[] psi2[is];
        delete[] psi3[is];
    }
    delete[] K1;
    delete[] K2;
    delete[] K3;
    delete[] K4;
    delete[] psi1;
    delete[] psi2;
    delete[] psi3;
    ModuleBase::timer::tick("ESolver_OF_TDDFT", "propagte_psi");
}

} // namespace ModuleESolver
