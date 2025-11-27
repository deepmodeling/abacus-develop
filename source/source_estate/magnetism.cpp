#include "magnetism.h"

#include "source_base/parallel_reduce.h"
#include "source_io/module_parameter/parameter.h"
//#include "source_estate/module_charge/charge.h"

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
    assert(omega>0.0);
    assert(nxyz>0);

    const double fac = omega / nxyz;

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
        this->tot_mag *= fac;
        this->abs_mag *= fac;

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
			const double r1 = rho[1][ir];
			const double r2 = rho[2][ir];
			const double r3 = rho[3][ir];
			double diff = sqrt(r1*r1 + r2*r2 + r3*r3);
 
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
			this->tot_mag_nc[i] *= fac;
            // mohan add 2025-06-21
			if( std::abs(this->tot_mag_nc[i]) < 1.0e-16)
			{
				this->tot_mag_nc[i] = 0.0;
			}
		}

		this->abs_mag *= fac;

        // mohan update 2025-06-21
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total magnetism (Bohr mag/cell)",
                                   this->tot_mag_nc[0], this->tot_mag_nc[1], this->tot_mag_nc[2]);

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Absolute magnetism (Bohr mag/cell)",this->abs_mag);
	}

    return;
}


bool Magnetism::judge_parallel(const double a[3], const ModuleBase::Vector3<double> &b)
{
   bool jp=false;

   double cross=0.0;

   const double c1 = a[1]*b.z - a[2]*b.y;
   const double c2 = a[2]*b.x - a[0]*b.z;
   const double c3 = a[0]*b.y - a[1]*b.x;
   cross = c1*c1 + c2*c2 + c3*c3;

   jp = (fabs(cross)<1e-6);
   return jp;
}
