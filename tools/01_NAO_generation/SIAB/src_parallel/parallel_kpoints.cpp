#include "parallel_kpoints.h"
#include "parallel_global.h"
#include "parallel_common.h"

Parallel_Kpoints::Parallel_Kpoints()
{
	nks_pool = new int[1];
	startk_pool = new int[1];
	whichpool = new int[1];
}
	
Parallel_Kpoints::~Parallel_Kpoints()
{
	delete[] nks_pool;
	delete[] startk_pool;
	delete[] whichpool;
}

void Parallel_Kpoints::init(void)
{
	//TITLE("Parallel_Kpoints","init");
#ifdef __MPI
//----------------------------------------------------------
// CALL Function : divide_pools 
//----------------------------------------------------------
	this->divide_pools();

// for test

	if(MY_RANK==0)
	{		
		cout << "\n    " << setw(8) << "MY_RANK"
			<< setw(8) << "MY_POOL"
			<< setw(13) << "RANK_IN_POOL"
			<< setw(6) << "NPROC"
			<< setw(6) << "KPAR"
			<< setw(14) << "NPROC_IN_POOL" << endl;
	}

		ofs_running << "\n    " << setw(8) << "MY_RANK"
			<< setw(8) << "MY_POOL"
			<< setw(13) << "RANK_IN_POOL"
			<< setw(6) << "NPROC"
			<< setw(6) << "KPAR"
			<< setw(14) << "NPROC_IN_POOL" << endl;
	
	for(int i=0;i<NPROC;i++)
	{	
		if(MY_RANK == i)
		{
			cout << " I'm" << setw(8) << MY_RANK
			<< setw(8) << MY_POOL
			<< setw(13) << RANK_IN_POOL
			<< setw(6) << NPROC
			<< setw(6) << KPAR
			<< setw(14) << NPROC_IN_POOL << endl;

			ofs_running << " I'm" << setw(8) << MY_RANK
			<< setw(8) << MY_POOL
			<< setw(13) << RANK_IN_POOL
			<< setw(6) << NPROC
			<< setw(6) << KPAR
			<< setw(14) << NPROC_IN_POOL << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if(MY_RANK != 0 )
	{
		cout.rdbuf(NULL);
	}

	return;
#endif
}

void Parallel_Kpoints::divide_pools(void)
{

#ifdef __MPI
	//cout<<"\n ==> mpi_split()"<<endl;
	int i=0;
	int j=0;
	if(NPROC<KPAR)
	{
		cout<<"\n NPROC=" << NPROC << " KPAR=" << KPAR;
		cout<<"Error : Too many pools !"<<endl;
		exit(0);
	}
	//if(kv.nkstot<KPAR) cout<<"Error !"<<endl;
	
	// (1) per process in each pool
	NPROC_IN_POOL = NPROC/KPAR;
	if(MY_RANK < (NPROC%KPAR)*(NPROC_IN_POOL+1)) 
	{
		NPROC_IN_POOL++;
	}

	// (2) To know how many process in pool j.
	nproc_pool = new int[KPAR];
	ZEROS(nproc_pool, KPAR);
	for(i=0;i<NPROC;i++)
	{
		j = i%KPAR;
		nproc_pool[j]++;
	}

	// (3) To know start proc index in each pool.
	startpro_pool = new int[KPAR];
	ZEROS(startpro_pool, KPAR);
	for(i=1;i<KPAR;i++) 
	{
		startpro_pool[i]=startpro_pool[i-1]+nproc_pool[i-1];
	}

	// use 'MY_RANK' to know 'MY_POOL'.
	for(i=0;i<KPAR;i++) 
	{
		if(MY_RANK >= startpro_pool[i]) 
		{
			MY_POOL=i;
		}
	}

	int key = 1;
	RANK_IN_POOL = MY_RANK-startpro_pool[MY_POOL];
	
	//========================================================
	// MPI_Comm_Split: Creates new communicators based on 
	// colors(2nd parameter) and keys(3rd parameter)	
	// Note: The color must be non-negative or MPI_UNDEFINED.
	//========================================================	
	MPI_Comm_split(MPI_COMM_WORLD,MY_POOL,key,&POOL_WORLD);
#endif
	
	return;
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

void Parallel_Kpoints::get_whichpool(const int &nkstot)
{
	delete[] whichpool;
	this->whichpool = new int[nkstot];
	ZEROS(whichpool, nkstot);

	for(int i=0; i<KPAR; i++)
	{
		for(int ik=0; ik< this->nks_pool[i]; ik++)
		{
			const int k_now = ik + startk_pool[i];
			this->whichpool[k_now] = i; 
			//ofs_running << "\n whichpool[" << k_now <<"] = " << whichpool[k_now];
		}
	}
	
	return;
}

void Parallel_Kpoints::get_nks_pool(const int &nkstot)
{
	delete[] nks_pool;
	this->nks_pool = new int[KPAR];
	ZEROS(nks_pool, KPAR);
	 
	const int nks_ave = nkstot/KPAR;
	const int remain = nkstot%KPAR;

	//ofs_running << "\n nkstot = " << nkstot;
	//ofs_running << "\n KPAR = " << KPAR;
	//ofs_running << "\n nks_ave = " << nks_ave;
	
	for(int i=0; i<KPAR; i++)
	{
		this->nks_pool[i] = nks_ave;
		if(i<remain)
		{
			nks_pool[i]++;
		}
		//ofs_running << "\n nks_pool[i] = " << nks_pool[i];
	}
	return;
}

void Parallel_Kpoints::get_startk_pool(const int &nkstot)
{
	delete[] startk_pool;
	startk_pool = new int[KPAR]; 
	const int remain = nkstot%KPAR;
	
	startk_pool[0] = 0;
	for(int i=1; i<KPAR; i++)
	{
		startk_pool[i] = startk_pool[i-1] + nks_pool[i-1];
		//ofs_running << "\n startk_pool[i] = " << startk_pool[i];
	}
	return;
}


void Parallel_Kpoints::pool_collection(double &value, const double *wk, const int &ik)
{
#ifdef __MPI
	const int ik_now = ik - this->startk_pool[MY_POOL];
	MPI_Status* ierror;
	if(this->whichpool[ik] == MY_POOL)
	{
		if(MY_POOL > 0 && RANK_IN_POOL == 0)
		{
			value = wk[ik_now];
			MPI_Send(&value, 1, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD);
		}
		else if(MY_POOL == 0 && MY_RANK == 0)
		{
			value = wk[ik_now];
		}
	}
	else
	{
		if(MY_RANK==0)
		{
			const int iproc = this->startpro_pool[ this->whichpool[ik] ];
			MPI_Recv(&value, 1, MPI_DOUBLE, iproc, ik, MPI_COMM_WORLD,ierror);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
#else
	value = wk[ik];
#endif
	return;
}


void Parallel_Kpoints::pool_collection(double *value, const realArray& overlap, const int &ik)
{
	const int dim2 = overlap.getBound2();
	const int dim3 = overlap.getBound3();
	const int dim4 = overlap.getBound4();
	const int dim = dim2 * dim3 * dim4;

#ifdef __MPI
	const int ik_now = ik - this->startk_pool[MY_POOL];
	const int begin = ik_now * dim2 * dim3 * dim4;
	double* pointer = &overlap.ptr[begin];
	
	MPI_Status* ierror;
	if(this->whichpool[ik] == MY_POOL)
	{
		if(MY_POOL > 0 && RANK_IN_POOL == 0)
		{
			// data transfer begin.
			for(int i=0; i<dim; i++)
			{
				value[i] = *pointer;
				++pointer;
			}
			// data transfer ends.
			MPI_Send(value, dim, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD);
		}
		else if(MY_POOL == 0 && MY_RANK == 0)
		{
//			cout << "\n ik = " << ik << endl;
			// data transfer begin.
			for(int i=0; i<dim; i++)
			{
				value[i] = *pointer;
				++pointer;
			}
			// data transfer ends.
		}	
	}
	else
	{
		if(MY_RANK==0)
		{
			const int iproc = this->startpro_pool[ this->whichpool[ik] ];
			MPI_Recv(value, dim, MPI_DOUBLE, iproc, ik, MPI_COMM_WORLD,ierror);
		}
	}
#else
	// data transfer ends.
	const int begin = ik * dim2 * dim3 * dim4;
	double* pointer = &overlap.ptr[begin];
	for(int i=0; i<dim; i++)
	{
		value[i] = *pointer;
		++pointer;
	}		
	// data transfer ends.
#endif
	return;
}
