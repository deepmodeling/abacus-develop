#include <vector>

#include "charge.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/libm/libm.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_elecstate/magnetism.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_io/rho_io.h"

void Charge::init_rho(elecstate::efermi& eferm_iout, const ModuleBase::ComplexMatrix& strucFac)
{
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "init_chg", GlobalV::init_chg);

    std::cout << " START CHARGE      : " << GlobalV::init_chg << std::endl;
    if (GlobalV::init_chg == "atomic") // mohan add 2007-10-17
    {
        this->atomic_rho(GlobalV::NSPIN, GlobalC::ucell.omega, rho, strucFac, GlobalC::ucell);
    }
    else if (GlobalV::init_chg == "file")
    {
        GlobalV::ofs_running << " try to read charge from file : ";
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            std::stringstream ssc;
            ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_CHG.cube";
            GlobalV::ofs_running << ssc.str() << std::endl;
            double& ef_tmp = eferm_iout.get_ef(is);
            if (ModuleIO::read_rho(
#ifdef __MPI
                    &(GlobalC::Pgrid),
#endif
                    is,
                    GlobalV::NSPIN,
                    ssc.str(),
                    this->rho[is],
                    this->rhopw->nx,
                    this->rhopw->ny,
                    this->rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell),
                    this->prenspin))
            {
                GlobalV::ofs_running << " Read in the charge density: " << ssc.str() << std::endl;
            }
            else if (is > 0)
            {
                if (prenspin == 1)
                {
                    GlobalV::ofs_running << " Didn't read in the charge density but autoset it for spin " << is + 1
                                         << std::endl;
                    for (int ir = 0; ir < this->rhopw->nrxx; ir++)
                    {
                        this->rho[is][ir] = 0.0;
                    }
                }
                //
                else if (prenspin == 2)
                { // read up and down , then rearrange them.
                    if (is == 1)
                    {
                        ModuleBase::WARNING_QUIT("Charge::init_rho", "Incomplete charge density file!");
                    }
                    else if (is == 2)
                    {
                        GlobalV::ofs_running << " Didn't read in the charge density but would rearrange it later. "
                                             << std::endl;
                    }
                    else if (is == 3)
                    {
                        GlobalV::ofs_running << " rearrange charge density " << std::endl;
                        for (int ir = 0; ir < this->rhopw->nrxx; ir++)
                        {
                            this->rho[3][ir] = this->rho[0][ir] - this->rho[1][ir];
                            this->rho[0][ir] = this->rho[0][ir] + this->rho[1][ir];
                            this->rho[1][ir] = 0.0;
                            this->rho[2][ir] = 0.0;
                        }
                    }
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT(
                    "init_rho",
                    "!!! Couldn't find the charge file !!! The default directory \n of SPIN1_CHG.cube is OUT.suffix, "
                    "or you must set read_file_dir \n to a specific directory. ");
            }
        }

        if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
        {
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                std::stringstream ssc;
                ssc << GlobalV::global_readin_dir << "SPIN" << is + 1 << "_TAU.cube";
                GlobalV::ofs_running << " try to read kinetic energy density from file : " << ssc.str() << std::endl;
                // mohan update 2012-02-10, sunliang update 2023-03-09
                if (ModuleIO::read_rho(
#ifdef __MPI
                        &(GlobalC::Pgrid),
#endif
                        is,
                        GlobalV::NSPIN,
                        ssc.str(),
                        this->kin_r[is],
                        this->rhopw->nx,
                        this->rhopw->ny,
                        this->rhopw->nz,
                        eferm_iout.ef,
                        &(GlobalC::ucell),
                        this->prenspin))
                {
                    GlobalV::ofs_running << " Read in the kinetic energy density: " << ssc.str() << std::endl;
                }
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("Charge::init_rho", "init_chg is wrong!");
    }

    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_load.load_charge && !GlobalC::restart.info_load.load_charge_finish)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            GlobalC::restart.load_disk("charge", is, this->nrxx, rho);
        }
        GlobalC::restart.info_load.load_charge_finish = true;
    }
}

//==========================================================
// computes the core charge on the real space 3D mesh.
//==========================================================
void Charge::set_rho_core(
    const ModuleBase::ComplexMatrix &structure_factor
)
{
    ModuleBase::TITLE("Charge","set_rho_core");
    ModuleBase::timer::tick("Charge","set_rho_core");

    // double eps = 1.e-10;
    //----------------------------------------------------------
    // LOCAL VARIABLES :
    // counter on mesh points
    // counter on atomic types
    // counter on g vectors
    //----------------------------------------------------------
    // int ir = 0;
    // int it = 0;
    // int ig = 0;

    bool bl = false;
    for (int it = 0; it<GlobalC::ucell.ntype; it++)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
            bl = true;
            break;
        }
    }

    if (!bl)
    {
        ModuleBase::GlobalFunc::ZEROS( this->rho_core, this->rhopw->nrxx);
    	ModuleBase::timer::tick("Charge","set_rho_core");
        return;
    }

    double *rhocg = new double[this->rhopw->ngg];
    ModuleBase::GlobalFunc::ZEROS(rhocg, this->rhopw->ngg );

	// three dimension.
    std::complex<double> *vg = new std::complex<double>[this->rhopw->npw];	

    for (int it = 0; it < GlobalC::ucell.ntype;it++)
    {
        if (GlobalC::ucell.atoms[it].ncpp.nlcc)
        {
//----------------------------------------------------------
// EXPLAIN : drhoc compute the radial fourier transform for
// each shell of g vec
//----------------------------------------------------------
            this->non_linear_core_correction(
                GlobalC::ppcell.numeric,
                GlobalC::ucell.atoms[it].ncpp.msh,
                GlobalC::ucell.atoms[it].ncpp.r,
                GlobalC::ucell.atoms[it].ncpp.rab,
                GlobalC::ucell.atoms[it].ncpp.rho_atc,
                rhocg);
//----------------------------------------------------------
// EXPLAIN : multiply by the structure factor and sum
//----------------------------------------------------------
            for (int ig = 0; ig < this->rhopw->npw ; ig++)
            {
                vg[ig] += structure_factor(it, ig) * rhocg[this->rhopw->ig2igg[ig]];
            }
        }
    }

	// for tmp use.
	for(int ig=0; ig< this->rhopw->npw; ig++)
	{
		this->rhog_core[ig] = vg[ig];
	}

    this->rhopw->recip2real(vg, this->rho_core);

    // test on the charge and computation of the core energy
    double rhoima = 0.0;
    double rhoneg = 0.0;
    for (int ir = 0; ir < this->rhopw->nrxx; ir++)
    {
        rhoneg += std::min(0.0, this->rhopw->ft.get_auxr_data<double>()[ir].real());
        rhoima += std::abs(this->rhopw->ft.get_auxr_data<double>()[ir].imag());
        // NOTE: Core charge is computed in reciprocal space and brought to real
        // space by FFT. For non smooth core charges (or insufficient cut-off)
        // this may result in negative values in some grid points.
        // Up to October 1999 the core charge was forced to be positive definite.
        // This induces an error in the force, and probably stress, calculation if
        // the number of grid points where the core charge would be otherwise neg
        // is large. The error disappears for sufficiently high cut-off, but may be
        // rather large and it is better to leave the core charge as it is.
        // If you insist to have it positive definite (with the possible problems
        // mentioned above) uncomment the following lines.  SdG, Oct 15 1999
    }

	// mohan fix bug 2011-04-03
	Parallel_Reduce::reduce_double_pool( rhoneg );
	Parallel_Reduce::reduce_double_pool( rhoima );

	// mohan changed 2010-2-2, make this same as in atomic_rho.
	// still lack something......
    rhoneg /= this->rhopw->nxyz * GlobalC::ucell.omega;
    rhoima /= this->rhopw->nxyz * GlobalC::ucell.omega;

    // calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
    // The term was present in previous versions of the code but it shouldn't
    delete [] rhocg;
    delete [] vg;
    ModuleBase::timer::tick("Charge","set_rho_core");
    return;
} // end subroutine set_rhoc

void Charge::non_linear_core_correction
(
    const bool &numeric,
    const int mesh,
    const double *r,
    const double *rab,
    const double *rhoc,
    double *rhocg) const
{
    ModuleBase::TITLE("charge","drhoc");

	// use labmda instead of repeating codes
	const auto kernel = [&](int num_threads, int thread_id)
	{

	double gx = 0.0;
    double rhocg1 = 0.0;
    double *aux;

    // here we compute the fourier transform is the charge in numeric form
    if (numeric)
    {
        aux = new double [mesh];
        // G=0 term

        int igl0 = 0;
        if (this->rhopw->gg_uniq [0] < 1.0e-8)
        {
			// single thread term
			if (thread_id == 0)
			{
				for (int ir = 0;ir < mesh; ir++)
				{
					aux [ir] = r [ir] * r [ir] * rhoc [ir];
				}
				ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
				//rhocg [1] = fpi * rhocg1 / omega;
				rhocg [0] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;//mohan modify 2008-01-19
			}
            igl0 = 1;
        }

		int igl_beg, igl_end;
		// exclude igl0
		ModuleBase::TASK_DIST_1D(num_threads, thread_id, this->rhopw->ngg - igl0, igl_beg, igl_end);
		igl_beg += igl0;
		igl_end += igl_beg;

        // G <> 0 term
        for (int igl = igl_beg; igl < igl_end;igl++) 
        {
            gx = sqrt(this->rhopw->gg_uniq[igl] * GlobalC::ucell.tpiba2);
            ModuleBase::Sphbes::Spherical_Bessel(mesh, r, gx, 0, aux);
            for (int ir = 0;ir < mesh; ir++) 
            {
                aux [ir] = r[ir] * r[ir] * rhoc [ir] * aux [ir];
            } //  enddo
            ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
            rhocg [igl] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;
        } //  enddo
        delete [] aux;
    }
    else
    {
        // here the case where the charge is in analytic form,
        // check old version before 2008-12-9
    }

	}; // end kernel

	// do not use omp parallel when this function is already in parallel block
	// 
	// it is called in parallel block in Forces::cal_force_cc,
	// but not in other funtcion such as Stress_Func::stress_cc.
	ModuleBase::TRY_OMP_PARALLEL(kernel);

    return;
}