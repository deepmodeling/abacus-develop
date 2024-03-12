#include "parallel_kpoints.h"
#include "module_base/parallel_global.h"
#include "module_base/parallel_common.h"

Parallel_Kpoints::Parallel_Kpoints()
{
    nks_pool = nullptr;
    startk_pool = nullptr;
    whichpool = nullptr;
}

Parallel_Kpoints::~Parallel_Kpoints()
{
    delete[] nks_pool;
    delete[] startk_pool;
    delete[] whichpool;
}


// the kpoints here are reduced after symmetry applied.
void Parallel_Kpoints::kinfo(int &nkstot)
{
#ifdef __MPI
    Parallel_Common::bcast_int(nkstot);
    this->get_nks_pool(nkstot);
    this->get_startk_pool(nkstot);
    this->get_whichpool(nkstot);
#endif
    return;
}

#ifdef __MPI
void Parallel_Kpoints::get_whichpool(const int &nkstot)
{
    delete[] whichpool;
    this->whichpool = new int[nkstot];
    ModuleBase::GlobalFunc::ZEROS(whichpool, nkstot);

	//std::cout << " calculate : whichpool" << std::endl;
	//std::cout << " nkstot is " << nkstot << std::endl;

    for (int i=0; i<GlobalV::KPAR; i++)
    {
        for (int ik=0; ik< this->nks_pool[i]; ik++)
        {
            const int k_now = ik + startk_pool[i];
            this->whichpool[k_now] = i;
            //GlobalV::ofs_running << "\n whichpool[" << k_now <<"] = " << whichpool[k_now];
        }
    }

    return;
}

void Parallel_Kpoints::get_nks_pool(const int &nkstot)
{
    delete[] nks_pool;
    this->nks_pool = new int[GlobalV::KPAR];
    ModuleBase::GlobalFunc::ZEROS(nks_pool, GlobalV::KPAR);

    const int nks_ave = nkstot/GlobalV::KPAR;
    const int remain = nkstot%GlobalV::KPAR;

    //GlobalV::ofs_running << "\n nkstot = " << nkstot;
    //GlobalV::ofs_running << "\n GlobalV::KPAR = " << GlobalV::KPAR;
    //GlobalV::ofs_running << "\n nks_ave = " << nks_ave;

    for (int i=0; i<GlobalV::KPAR; i++)
    {
        this->nks_pool[i] = nks_ave;
        if (i<remain)
        {
            nks_pool[i]++;
        }
        //GlobalV::ofs_running << "\n nks_pool[i] = " << nks_pool[i];
    }
    return;
}

void Parallel_Kpoints::get_startk_pool(const int &nkstot)
{
    delete[] startk_pool;
    startk_pool = new int[GlobalV::KPAR];
    //const int remain = nkstot%GlobalV::KPAR;

    startk_pool[0] = 0;
    for (int i=1; i<GlobalV::KPAR; i++)
    {
        startk_pool[i] = startk_pool[i-1] + nks_pool[i-1];
        //GlobalV::ofs_running << "\n startk_pool[i] = " << startk_pool[i];
    }
    return;
}
#endif


void Parallel_Kpoints::pool_collection(double &value, const double *wk, const int &ik)
{
#ifdef __MPI

    if (GlobalV::RANK_IN_POOL==0)
    {
        const int ik_now = ik - this->startk_pool[GlobalV::MY_POOL];
        //GlobalV::ofs_running << "\n\n ik=" << ik << " ik_now=" << ik_now;
        if (this->whichpool[ik] == GlobalV::MY_POOL)
        {
            if (GlobalV::MY_POOL > 0)
            {
                value = wk[ik_now];
                MPI_Send(&value, 1, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD);
                //GlobalV::ofs_running << "\n send wk[" << ik << "]=" << value;
            }
            else
            {
                value = wk[ik_now];
                //GlobalV::ofs_running << "\n wk[" << ik << "]=" << value;
            }
        }
        else
        {
            if (GlobalV::MY_RANK==0)
            {
                MPI_Status ierror;
                const int iproc = this->startpro_pool[ this->whichpool[ik] ];
                MPI_Recv(&value, 1, MPI_DOUBLE, iproc, ik, MPI_COMM_WORLD,&ierror);
                //GlobalV::ofs_running << "\n receive wk[" << ik << "]=" << value << " from proc=" << iproc;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    value = wk[ik];
#endif
    return;
}


void Parallel_Kpoints::pool_collection(double *valuea, double *valueb, const ModuleBase::realArray &a, const ModuleBase::realArray &b, const int &ik)
{
    const int dim2 = a.getBound2();
    const int dim3 = a.getBound3();
    const int dim4 = a.getBound4();
    assert( a.getBound2() == b.getBound2() );
    assert( a.getBound3() == b.getBound3() );
    assert( a.getBound4() == b.getBound4() );
    const int dim = dim2 * dim3 * dim4;
#ifdef __MPI
    const int ik_now = ik - this->startk_pool[GlobalV::MY_POOL];
    const int begin = ik_now * dim2 * dim3 * dim4;
    double* pa = &a.ptr[begin];
    double* pb = &b.ptr[begin];

    const int pool = this->whichpool[ik];

	GlobalV::ofs_running << "\n ik=" << ik;

    if (GlobalV::RANK_IN_POOL==0)
    {
        if (GlobalV::MY_POOL==0)
        {
            if (pool==0)
            {
                for (int i=0; i<dim ; i++)
                {
                    valuea[i] = *pa;
                    valueb[i] = *pb;
                    ++pa;
                    ++pb;
                }
            }
            else
            {
				GlobalV::ofs_running << " receive data.";
                MPI_Status ierror;
                MPI_Recv(valuea, dim, MPI_DOUBLE, this->startpro_pool[pool], ik*2+0, MPI_COMM_WORLD,&ierror);
                MPI_Recv(valueb, dim, MPI_DOUBLE, this->startpro_pool[pool], ik*2+1, MPI_COMM_WORLD,&ierror);
            }
        }
        else
        {
            if (GlobalV::MY_POOL == pool)
            {
				GlobalV::ofs_running << " send data.";
                MPI_Send(pa, dim, MPI_DOUBLE, 0, ik*2+0, MPI_COMM_WORLD);
                MPI_Send(pb, dim, MPI_DOUBLE, 0, ik*2+1, MPI_COMM_WORLD);
            }
        }
    }
	else
	{
		GlobalV::ofs_running << "\n do nothing.";
	}
    MPI_Barrier(MPI_COMM_WORLD);

    /*


    	if(this->whichpool[ik] == GlobalV::MY_POOL)
    	{
    		if(GlobalV::MY_POOL > 0 && GlobalV::RANK_IN_POOL == 0)
    		{
    			// data transfer ends.
    			MPI_Send(pa, dim, MPI_DOUBLE, 0, ik*2,   MPI_COMM_WORLD);
    			MPI_Send(pb, dim, MPI_DOUBLE, 0, ik*2+1, MPI_COMM_WORLD);
    		}
    		else if(GlobalV::MY_POOL == 0 && GlobalV::MY_RANK == 0)
    		{
    //			std::cout << "\n ik = " << ik << std::endl;
    			// data transfer begin.
    			for(int i=0; i<dim; i++)
    			{
    				valuea[i] = *pa;
    				valueb[i] = *pb;
    				++pa;
    				++pb;
    			}
    			// data transfer ends.
    		}
    	}
    	else
    	{
    		if(GlobalV::MY_RANK==0)
    		{
    			MPI_Status* ierror;
    			const int iproc = this->startpro_pool[ this->whichpool[ik] ];
    			MPI_Recv(valuea, dim, MPI_DOUBLE, iproc, ik*2,   MPI_COMM_WORLD,ierror);
    			MPI_Recv(valueb, dim, MPI_DOUBLE, iproc, ik*2+1, MPI_COMM_WORLD,ierror);
    		}
    	}
    	*/
#else
    // data transfer ends.
    const int begin = ik * dim2 * dim3 * dim4;
    double* pa = &a.ptr[begin];
    double* pb = &b.ptr[begin];
    for (int i=0; i<dim; i++)
    {
        valuea[i] = *pa;
        valueb[i] = *pb;
        ++pa;
        ++pb;
    }
    // data transfer ends.
#endif
    return;
}


void Parallel_Kpoints::pool_collection(std::complex<double> *value, const ModuleBase::ComplexArray &w, const int &ik)
{
    const int dim2 = w.getBound2();
    const int dim3 = w.getBound3();
    const int dim4 = w.getBound4();
    const int dim = dim2 * dim3 * dim4;
#ifdef __MPI
    const int ik_now = ik - this->startk_pool[GlobalV::MY_POOL];
    const int begin = ik_now * dim2 * dim3 * dim4;
    std::complex<double>* p = &w.ptr[begin];
    //temprary restrict kpar=1 for NSPIN=2 case for generating_orbitals
    int pool = 0; 
    if(GlobalV::NSPIN != 2) pool = this->whichpool[ik];

	GlobalV::ofs_running << "\n ik=" << ik;

    if (GlobalV::RANK_IN_POOL==0)
    {
        if (GlobalV::MY_POOL==0)
        {
            if (pool==0)
            {
                for (int i=0; i<dim ; i++)
                {
                    value[i] = *p;
                    ++p;
                }
            }
            else
            {
				GlobalV::ofs_running << " receive data.";
                MPI_Status ierror;
                MPI_Recv(value, dim, MPI_DOUBLE, this->startpro_pool[pool], ik*2+0, MPI_COMM_WORLD,&ierror);
            }
        }
        else
        {
            if (GlobalV::MY_POOL == pool)
            {
				GlobalV::ofs_running << " send data.";
                MPI_Send(p, dim, MPI_DOUBLE, 0, ik*2+0, MPI_COMM_WORLD);
            }
        }
    }
	else
	{
		GlobalV::ofs_running << "\n do nothing.";
	}
    MPI_Barrier(MPI_COMM_WORLD);

#else
    // data transfer ends.
    const int begin = ik * dim2 * dim3 * dim4;
    std::complex<double> * p = &w.ptr[begin];
    for (int i=0; i<dim; i++)
    {
        value[i] = *p;
        ++p;
    }
    // data transfer ends.
#endif
}
