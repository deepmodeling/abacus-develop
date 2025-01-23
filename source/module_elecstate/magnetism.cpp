#include "magnetism.h"

#include "module_base/parallel_reduce.h"
#include "module_parameter/parameter.h"

Magnetism::Magnetism()
{
    tot_magnetization = 0.0;
    abs_magnetization = 0.0;
    std::fill(tot_magnetization_nc, tot_magnetization_nc + 3, 0.0);
    std::fill(ux_, ux_ + 3, 0.0);
}

Magnetism::~Magnetism()
{
    delete[] start_magnetization;
}

void Magnetism::compute_magnetization(const double& omega,
                                      const int& nrxx, 
                                      const int& nxyz, 
                                      const double* const * rho, 
                                      double* nelec_spin)
{
    if (PARAM.inp.nspin==2)
    {
        this->tot_magnetization = 0.00;
        this->abs_magnetization = 0.00;

        for (int ir=0; ir<nrxx; ir++)
        {
            double diff = rho[0][ir] - rho[1][ir];
            this->tot_magnetization += diff;
            this->abs_magnetization += std::abs(diff);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(this->tot_magnetization);
        Parallel_Reduce::reduce_pool(this->abs_magnetization);
#endif
        this->tot_magnetization *= omega / nxyz;
        this->abs_magnetization *= omega / nxyz;

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"total magnetism (Bohr mag/cell)",this->tot_magnetization);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
		
		//update number of electrons for each spin
		//if TWO_EFERMI, no need to update
		if(!PARAM.globalv.two_fermi)
		{
			nelec_spin[0] = (PARAM.inp.nelec + this->tot_magnetization) / 2;
			nelec_spin[1] = (PARAM.inp.nelec - this->tot_magnetization) / 2;
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nelec for spin up", nelec_spin[0]);
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nelec for spin down", nelec_spin[1]);
		}
    }

	// noncolliear :
	else if(PARAM.inp.nspin==4)
	{
		for(int i=0;i<3;i++) {this->tot_magnetization_nc[i] = 0.00;
}
		this->abs_magnetization = 0.00;
		for (int ir=0; ir<nrxx; ir++)
		{
			double diff = sqrt(pow(rho[1][ir], 2) + pow(rho[2][ir], 2) +pow(rho[3][ir], 2));
 
			for(int i=0;i<3;i++) {this->tot_magnetization_nc[i] += rho[i+1][ir];
}
			this->abs_magnetization += std::abs(diff);
		}
#ifdef __MPI
        Parallel_Reduce::reduce_pool(this->tot_magnetization_nc, 3);
        Parallel_Reduce::reduce_pool(this->abs_magnetization);
#endif
		for(int i=0;i<3;i++) {this->tot_magnetization_nc[i] *= omega/ nxyz;
}
		this->abs_magnetization *= omega/ nxyz;
		GlobalV::ofs_running<<"total magnetism (Bohr mag/cell)"<<'\t'<<this->tot_magnetization_nc[0]<<'\t'<<this->tot_magnetization_nc[1]<<'\t'<<this->tot_magnetization_nc[2]<<'\n';
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
	}

    return;
}

bool Magnetism::judge_parallel(double a[3], ModuleBase::Vector3<double> b)
{
   bool jp=false;
   double cross;
   cross = pow((a[1]*b.z-a[2]*b.y),2) +  pow((a[2]*b.x-a[0]*b.z),2) + pow((a[0]*b.y-a[1]*b.x),2);
   jp = (fabs(cross)<1e-6);
   return jp;
}
