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
    this->evolve_ofdft=new Evolve_OFDFT();
}

ESolver_OF_TDDFT::~ESolver_OF_TDDFT()
{
    delete this->evolve_ofdft;
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
        this->phi_td.resize(PARAM.inp.nspin*this->pw_rho->nrxx);
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

            // update the rho and phi based on the direction and theta
            this->update_rho();

            this->iter_++;

            ESolver_FP::iter_finish(ucell, istep, this->iter_, conv_esolver);
        }

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
            {
                phi_td[is*this->pw_rho->nrxx+ir]=pphi_[is][ir];
            }
        }
    }
    else if ((istep<1) && PARAM.inp.init_chg == "file")
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
            {
                phi_td[is*this->pw_rho->nrxx+ir]=pphi_[is][ir];
            }
        } 
        conv_esolver=true;
    }
    else
    {
        this->evolve_ofdft->propagate_psi(this->pelec, this->chr, ucell, this->phi_td, this->pw_rho);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
            {
                pphi_[is][ir]=std::abs(phi_td[is*this->pw_rho->nrxx+ir]);
            }
        }
        conv_esolver=true;
    }

    this->after_opt(istep, ucell, conv_esolver);

    ModuleBase::timer::tick("ESolver_OF_TDDFT", "runner");
}

} // namespace ModuleESolver
