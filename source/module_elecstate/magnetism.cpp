#include "magnetism.h"

#include "source_base/parallel_reduce.h"
#include "module_parameter/parameter.h"
//#include "module_elecstate/module_charge/charge.h"

Magnetism::Magnetism()
{
    tot_mag = 0.0;
    abs_mag = 0.0;
    std::fill(tot_mag_nc, tot_mag_nc + 3, 0.0);
    std::fill(ux_, ux_ + 3, 0.0);
}

Magnetism::~Magnetism()
{
    delete[] start_mag;
}

void Magnetism::compute_mag(const double& omega,
		const int& nrxx, 
		const int& nxyz, 
		const double* const * rho, 
		double* nelec_spin)
{
    assert(nxyz>0);

    if (PARAM.inp.nspin==2)
    {
        this->tot_mag = 0.00;
        this->abs_mag = 0.00;

        for (int ir=0; ir<nrxx; ir++)
        {
            double diff = rho[0][ir] - rho[1][ir];
            this->tot_mag += diff;
            this->abs_mag += std::abs(diff);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(this->tot_mag);
        Parallel_Reduce::reduce_pool(this->abs_mag);
#endif
        this->tot_mag *= omega / nxyz;
        this->abs_mag *= omega / nxyz;

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total magnetism (Bohr mag/cell)",this->tot_mag);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Absolute magnetism (Bohr mag/cell)",this->abs_mag);
		
		//update number of electrons for each spin
		//if TWO_EFERMI, no need to update
		if(!PARAM.globalv.two_fermi)
		{
			nelec_spin[0] = (PARAM.inp.nelec + this->tot_mag) / 2;
			nelec_spin[1] = (PARAM.inp.nelec - this->tot_mag) / 2;
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Electron number for spin up", nelec_spin[0]);
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Electron number for spin down", nelec_spin[1]);
		}
    }

	// noncolliear :
	else if(PARAM.inp.nspin==4)
	{
		for(int i=0;i<3;i++) 
		{
			this->tot_mag_nc[i] = 0.00;
		}

		this->abs_mag = 0.00;
		for (int ir=0; ir<nrxx; ir++)
		{
			double diff = sqrt(pow(rho[1][ir], 2) + pow(rho[2][ir], 2) +pow(rho[3][ir], 2));
 
			for(int i=0;i<3;i++) 
			{
				this->tot_mag_nc[i] += rho[i+1][ir];
			}
			this->abs_mag += std::abs(diff);
		}
#ifdef __MPI
        Parallel_Reduce::reduce_pool(this->tot_mag_nc, 3);
        Parallel_Reduce::reduce_pool(this->abs_mag);
#endif
		for(int i=0;i<3;i++) 
		{
			this->tot_mag_nc[i] *= omega/ nxyz;
		}

		this->abs_mag *= omega/ nxyz;

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total magnetism along x (Bohr mag/cell)",
                                   this->tot_mag_nc[0]);

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"                      y (Bohr mag/cell)",
                                   this->tot_mag_nc[1]);

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"                      z (Bohr mag/cell)",
                                   this->tot_mag_nc[2]);

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Absolute magnetism (Bohr mag/cell)",this->abs_mag);
	}

    return;
}


bool Magnetism::judge_parallel(const double a[3], const ModuleBase::Vector3<double> &b)
{
   bool jp=false;

   double cross=0.0;

   cross = pow((a[1]*b.z-a[2]*b.y),2) 
	   + pow((a[2]*b.x-a[0]*b.z),2) 
	   + pow((a[0]*b.y-a[1]*b.x),2);

   jp = (fabs(cross)<1e-6);
   return jp;
}
