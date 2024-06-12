#include "parallel_kpoints.h"
#include "module_base/parallel_global.h"
#include "module_base/parallel_common.h"

Parallel_Kpoints::Parallel_Kpoints()
{
}

Parallel_Kpoints::~Parallel_Kpoints()
{
}


// the kpoints here are reduced after symmetry applied.
void Parallel_Kpoints::kinfo(int &nkstot)
{
#ifdef __MPI
    this->kpar = GlobalV::KPAR;    // number of pools
    this->my_pool = GlobalV::MY_POOL;
    this->rank_in_pool = GlobalV::RANK_IN_POOL;
    this->nproc = GlobalV::NPROC;
    this->nspin = GlobalV::NSPIN;

    Parallel_Common::bcast_int(nkstot);
    this->get_nks_pool(nkstot);     // assign k-points to each pool
    this->get_startk_pool(nkstot);  // get the start k-point index for each pool
    this->get_whichpool(nkstot);    // get the pool index for each k-point

    this->set_startpro_pool();      // get the start processor index for each pool
    
    this->nkstot_np = nkstot;       
    this->nks_np = this->nks_pool[this->my_pool]; // number of k-points in this pool
#else
    this->kpar = 1;
    this->my_pool = 0;
    this->rank_in_pool = 0;
    this->nproc = 1;
    this->nspin = GlobalV::NSPIN;
    this->nkstot_np = nkstot;
    this->nks_np = nkstot;
#endif
    return;
}

#ifdef __MPI
void Parallel_Kpoints::get_whichpool(const int &nkstot)
{
    this->whichpool.resize(nkstot, 0);

	//std::cout << " calculate : whichpool" << std::endl;
	//std::cout << " nkstot is " << nkstot << std::endl;

    for (int i=0; i<this->kpar; i++)
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
    nks_pool.resize(this->kpar, 0);

    const int nks_ave = nkstot/this->kpar;
    const int remain = nkstot%this->kpar;

    //GlobalV::ofs_running << "\n nkstot = " << nkstot;
    //GlobalV::ofs_running << "\n this->kpar = " << this->kpar;
    //GlobalV::ofs_running << "\n nks_ave = " << nks_ave;

    for (int i=0; i<this->kpar; i++)
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
    startk_pool.resize(this->kpar, 0);
    //const int remain = nkstot%this->kpar;

    startk_pool[0] = 0;
    for (int i=1; i<this->kpar; i++)
    {
        startk_pool[i] = startk_pool[i-1] + nks_pool[i-1];
        //GlobalV::ofs_running << "\n startk_pool[i] = " << startk_pool[i];
    }
    return;
}

void Parallel_Kpoints::set_startpro_pool(void)
{
    startpro_pool.resize(this->kpar, 0);

    const int nproc_ave = this->nproc/this->kpar;
    const int remain = this->nproc%this->kpar;

    startpro_pool[0] = 0;
    for (int i=1; i<this->kpar; i++)
    {
        startpro_pool[i] = startpro_pool[i-1] + nproc_ave;
        if (i-1<remain)
        {
            startpro_pool[i]++;
        }
        //GlobalV::ofs_running << "\n startpro_pool[i] = " << startpro_pool[i];
    }
    return;
}

// gather kpoints from all processor pools, only need to be called by the first processor of each pool.
void Parallel_Kpoints::gatherkvec(const std::vector<ModuleBase::Vector3<double>>& vec_local,
                                  std::vector<ModuleBase::Vector3<double>>& vec_global) const
{
    vec_global.resize(this->nkstot_np, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    for (int i = 0; i < this->nks_np; ++i)
    {
        if (this->rank_in_pool==0)
        {
            vec_global[i + startk_pool[this->my_pool]] = vec_local[i];
        }
        //vec_global[i + startk_pool[GlobalV::MY_POOL]] = vec_local[i] / double(GlobalV::NPROC_IN_POOL);
    }

    MPI_Allreduce(MPI_IN_PLACE, &vec_global[0], 3 * this->nkstot_np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return;
}
#endif


void Parallel_Kpoints::pool_collection(double &value, const double *wk, const int &ik)
{
#ifdef __MPI
     const int ik_now = ik - this->startk_pool[this->my_pool];
    //GlobalV::ofs_running << "\n\n ik=" << ik << " ik_now=" << ik_now;

    const int pool = this->whichpool[ik];
    
    if (this->rank_in_pool==0)
    {
        if (this->my_pool==0)
        {
            if (pool==0)
            {
                value = wk[ik_now];
            }
            else
            {
                GlobalV::ofs_running << " receive data.";
                MPI_Status ierror;
                MPI_Recv(&value, 1, MPI_DOUBLE, this->startpro_pool[pool], ik, MPI_COMM_WORLD,&ierror);
            }
        }
        else
        {
            if (this->my_pool == pool)
            {
                GlobalV::ofs_running << " send data.";
                MPI_Send(&wk[ik_now], 1, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        GlobalV::ofs_running << "\n do nothing.";
    }

    MPI_Barrier(MPI_COMM_WORLD);
#else
    value = wk[ik];
#endif
    return;
}


void Parallel_Kpoints::pool_collection(double *value_re, double *value_im, const ModuleBase::realArray &re, const ModuleBase::realArray &im, const int &ik)
{
    const int dim2 = re.getBound2();
    const int dim3 = re.getBound3();
    const int dim4 = re.getBound4();
    assert( re.getBound2() == im.getBound2() );
    assert( re.getBound3() == im.getBound3() );
    assert( re.getBound4() == im.getBound4() );
    const int dim = dim2 * dim3 * dim4;
    pool_collection_aux(value_re, re, dim, ik);
    pool_collection_aux(value_im, im, dim, ik);
    return;
}


void Parallel_Kpoints::pool_collection(std::complex<double> *value, const ModuleBase::ComplexArray &w, const int &ik)
{
    const int dim2 = w.getBound2();
    const int dim3 = w.getBound3();
    const int dim4 = w.getBound4();
    const int dim = dim2 * dim3 * dim4;
    pool_collection_aux(value, w, dim, ik);
}

template <class T, class V> void Parallel_Kpoints::pool_collection_aux(T *value, const V &w, const int& dim, const int &ik)
{
#ifdef __MPI
    const int ik_now = ik - this->startk_pool[this->my_pool];
    const int begin = ik_now * dim;
    T* p = &w.ptr[begin];
    //temprary restrict kpar=1 for NSPIN=2 case for generating_orbitals
    int pool = 0;
    if(this->nspin != 2) pool = this->whichpool[ik];

	GlobalV::ofs_running << "\n ik=" << ik;

    if (this->rank_in_pool==0)
    {
        if (this->my_pool==0)
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
            if (this->my_pool == pool)
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
    const int begin = ik * dim;
    T * p = &w.ptr[begin];
    for (int i=0; i<dim; i++)
    {
        value[i] = *p;
        ++p;
    }
    // data transfer ends.
#endif
}