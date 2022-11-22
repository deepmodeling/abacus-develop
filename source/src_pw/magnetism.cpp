#include "magnetism.h"
#ifndef __CELL
#include "global.h"
#endif
#include "../src_parallel/parallel_reduce.h"

Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}

Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

#ifndef __CELL
void Magnetism::compute_magnetization(const Charge* const chr)
{
    if (GlobalV::NSPIN==2)
    {
        this->tot_magnetization = 0.00;
        this->abs_magnetization = 0.00;

        for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
        {
            double diff = chr->rho[0][ir] - chr->rho[1][ir];
            this->tot_magnetization += diff;
            this->abs_magnetization += abs(diff);
        }

        Parallel_Reduce::reduce_double_pool( this->tot_magnetization );
        Parallel_Reduce::reduce_double_pool( this->abs_magnetization );
        this->tot_magnetization *= GlobalC::ucell.omega / GlobalC::rhopw->nxyz;
        this->abs_magnetization *= GlobalC::ucell.omega / GlobalC::rhopw->nxyz;

		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"total magnetism (Bohr mag/cell)",this->tot_magnetization);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
		
		if(GlobalV::TWO_EFERMI)
		{
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nelup",get_nelup());
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"neldw",get_neldw());
		}
		else
		{
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nelec",GlobalV::nelec);
		}
    }

	// noncolliear :
	else if(GlobalV::NSPIN==4)
	{
		for(int i=0;i<3;i++)this->tot_magnetization_nc[i] = 0.00;
		this->abs_magnetization = 0.00;
		for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
		{
			double diff = sqrt(pow(chr->rho[1][ir], 2) + pow(chr->rho[2][ir], 2) +pow(chr->rho[3][ir], 2));
 
			for(int i=0;i<3;i++)this->tot_magnetization_nc[i] += chr->rho[i+1][ir];
			this->abs_magnetization += abs(diff);
		}
		Parallel_Reduce::reduce_double_pool( this->tot_magnetization_nc, 3 );
		Parallel_Reduce::reduce_double_pool( this->abs_magnetization );

		for(int i=0;i<3;i++)this->tot_magnetization_nc[i] *= GlobalC::ucell.omega / GlobalC::rhopw->nxyz;
		this->abs_magnetization *= GlobalC::ucell.omega / GlobalC::rhopw->nxyz;
		GlobalV::ofs_running<<"total magnetism (Bohr mag/cell)"<<'\t'<<this->tot_magnetization_nc[0]<<'\t'<<this->tot_magnetization_nc[1]<<'\t'<<this->tot_magnetization_nc[2]<<'\n';
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
	}

    return;
}


double Magnetism::get_nelup(void)
{
	double nelup = 0.0;
	if(GlobalV::TWO_EFERMI)
	{
//===============================================================
//  this type of electrons are used as "fixed" magnetization.
//===============================================================
		nelup = 0.5 * GlobalV::nelec + 0.5 * tot_magnetization;
	}
	else
	{
		nelup = 0.5 * GlobalV::nelec;
	}
    return nelup;
}


double Magnetism::get_neldw(void)
{
	double neldw = 0.0;
	if(GlobalV::TWO_EFERMI)
	{
//===============================================================
//  this type of electrons are used as "fixed" magnetization.
//===============================================================
		neldw = 0.5 * GlobalV::nelec - 0.5 * tot_magnetization;
	}
	else
	{
		neldw = 0.5 * GlobalV::nelec;
	}
    return neldw ;
}

bool Magnetism::judge_parallel(double a[3], ModuleBase::Vector3<double> b)
{
   bool jp=false;
   double cross;
   cross = pow((a[1]*b.z-a[2]*b.y),2) +  pow((a[2]*b.x-a[0]*b.z),2) + pow((a[0]*b.y-a[1]*b.x),2);
   jp = (fabs(cross)<1e-6);
   return jp;
}
#endif