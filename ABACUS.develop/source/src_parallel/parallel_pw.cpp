#include "../src_pw/global.h"	// only alat,Ecut,wfac,Nx,Ny,Nz,kv,out
#include "../src_pw/tools.h"
#include "parallel_pw.h"

Parallel_PW::Parallel_PW()
{
//-------------------------------------
// SURVIVE TIME: 
//-------------------------------------
	this->ngm_i = 0;
	this->ngm_i_record = 0;
	this->allocate_igl2g = false;

	this->npps = new int[1];
	this->st_start = new int[1];
	this->nst_per = new int[1];
	this->npw_per = new int[1];

	this->ig_l2g = new int[1];
	
	this->index = new bool[1];
	this->index_ip = new int[1];

	this->isind = new int[1];
	this->ismap = new int[1];

    this->ngm_i_number = new int[1];
    this->ngm_i2 = 0;
    this->ngm_i_record2 = 0;
    this->allocate_igl2g_final_scf = false; //LiuXh add 20180619
    this->ngm_i_final_scf = 0;
    this->ngm_i_record_final_scf = 0;
}

Parallel_PW::~Parallel_PW()
{
	delete[] npps;
	delete[] st_start;
	delete[] nst_per;
	delete[] npw_per;
	
	delete [] ig_l2g;
	
	for(int i=0;i<this->n1;i++)
	{
		delete[] st[i];
	}
	
	delete[] st;
	delete[] st_i;
	delete[] st_j;
	delete[] st_n;
	
	delete[] index;
	delete[] index_ip;
	
	delete[] isind;
	delete[] ismap;
	
delete[] ngm_i_number;
}

void Parallel_PW::init(const double &gcut_in, const int &n1_in, const int &n2_in, const int &n3_in,
		const int &bz, const int &nproc_in,const int &rank_in)
{
	if(test_pw) TITLE("Parallel_PW","init");

	this->gcut=gcut_in;
	this->n1=n1_in;
	this->n2=n2_in;
	this->n3=n3_in;

	// nproc_use == nproc_p ( number of processes in this pool ).
	this->nproc_use = nproc_in;
	this->rank_use = rank_in;

	//==========================================================
	// npps : number of xy planes of each process in this pool. 
	//==========================================================
	delete[] this->npps;
	this->npps = new int[nproc_use];
	ZEROS(npps, nproc_use);

	const int nb3 = this->n3 / bz;
	
	for (int ip = 0;ip < nproc_use;ip++)
	{
		npps[ip] = nb3/nproc_use;
		npps[ip] *= bz;//mohan update 2011-04-22
		if(ip < nb3%nproc_use) this->npps[ip]+=bz;
		
		if(test_pw)
		{
			ofs_running << " processor="<<ip<<" npps="<<npps[ip] << endl;
		}
	}

	delete[] this->ngm_i_number;
	this->ngm_i_number = new int[10];
	ZEROS(ngm_i_number, 10);

	return;
}

void Parallel_PW::columns_map(void)
{
	if(test_pw) TITLE("Parallel_PW","columns_map");
	timer::tick("Parallel_PW","columns_map");
		
	int ibox[3];
	
	// set the center at origin point.
	ibox[0] = int(n1 / 2.0) + 1;
	ibox[1] = int(n2 / 2.0) + 1;
	ibox[2] = int(n3 / 2.0) + 1;

	//----------------------------------------------
	// st : number of sticks(n1*n2)
	// first dim: x
	// second dim: y
	// FUNCTION: count the number of nonzero element
	// on each grid(x,y)
	//----------------------------------------------
	this->st = new int *[n1];
	for (int i = 0;i < n1;i++)
	{
		this->st[i] = new int[n2];
		for(int j=0; j<n2; j++)
		{
			st[i][j] = 0;
		}
	}

	// get the number of plane wave within 'gcut'
	int ng = 0;
	// x is the slowest.
	for (int i = -ibox[0]; i <= ibox[0]; i++)
	{
		// y is always in the middle.
		for (int j = -ibox[1]; j <= ibox[1]; j++)
		{
			// z is the fastest.
			for (int k = -ibox[2]; k <= ibox[2]; k++)
			{
				Vector3<double> f(i,j,k);
				// g2= |f|^2 in the unit of (2Pi/lat0)^2
				double g2 = f * (ucell.GGT * f);

				// gcut is from input.
				if (g2 <= this->gcut)
				{
					// get the index of grid(x,y)
					int m1 = i;
					int m2 = j;

					if (i < 0) m1 = i + n1;
					if (j < 0) m2 = j + n2;

					assert(m1<n1);
					assert(m2<n2);//mohan 2009-11-9

					this->st[m1][m2]++;

					//cout<<setw(12)<<m1+m2*this->n1<<setw(12)<<this->st[m1][m2]++<<endl;
					ng++;
					// ng: plane wave number if spheri.
				}
			}
		}
	}

	// calculate the number of sticks.
	this->nst = 0;
	for (int i = 0;i < n1;i++)
	{
		for (int j = 0;j < n2;j++)
		{
			if (st[i][j] != 0)
			{
				this->nst++;
			}
		}
	}

	OUT(ofs_running,"number of plane waves",ng);
	OUT(ofs_running,"number of sticks", nst);

	timer::tick("Parallel_PW","columns_map");
	return;
}

void Parallel_PW::restore_st(void)
{
	if(test_pw) TITLE("Parallel_PW","restore_st");
	
	// assert the number of total sticks has
	// been calculated.
	assert(this->nst>0);

	if(this->nst < nproc_use)
	{
		WARNING_QUIT("Parallel_PW::restore_st","nst(number of sticks in FFT xy plane) < nproc_use, too many processes used.");
	}

	//---------------------------------------
	// record the information of each stick,
	// (1) location: st_i, st_j
	// (2) how many nonzero elements: st_n
	//---------------------------------------
	delete[] st_i;
	delete[] st_j;
	delete[] st_n;
	this->st_i = new int[nst];
	this->st_j = new int[nst];
	this->st_n = new int[nst];

	int k = 0;
	// First time we get the 'nst',
	// Here we search again in order
	// to restore the sticks.
	for (int i = 0;i < n1;i++)
	{
		for (int j = 0;j < n2;j++)
		{
			if (st[i][j] != 0)
			{
				this->st_i[k] = i; //i2
				this->st_j[k] = j; //i1 y
				this->st_n[k] = this->st[i][j]; 

				//cout<<setw(12)<<i<<setw(12)<<j<<setw(12)<<st_n[k]<<endl;
				//cout<<setw(12)<<i+j*this->n1<<setw(12)<<st_n[k]<<endl;
				k++;
			}
		}
	}
	return;
}


void Parallel_PW::columns_and_pw_distribution(void)
{
	if(test_pw) TITLE("Parallel_PW","columns_and_pw_distribution");

	// number of sticks of each processor in one pool.
	delete[] nst_per;
	this->nst_per = new int[nproc_use];
	ZEROS(nst_per, nproc_use);
	// number of planewaves of each processor in one pool.
	delete[] npw_per;
	this->npw_per = new int[nproc_use];
	ZEROS(npw_per, nproc_use);

	//******************************************************
	// Two choices  
	// (1): sticks found marked , 
	// (2): sort sticks 
	// Here we choose first ,which is easier to accompolish 
	// sticks found mark 
	//******************************************************
	delete[] index;
	this->index = new bool[this->nst];
	ZEROS(index,nst);

	delete[] index_ip;
	this->index_ip = new int[this->n1*this->n2];
	for (int i = 0;i < n1*n2 ;i++)
	{
		this->index_ip[i] = -1;
		//because processor begin from 0, so -1 is illegal number.
	}

	return;
}

void Parallel_PW::max_pw_column(int &pw, int &max_i, int &max_j)
{
//	cout<<"\n ==> Parallel_PW.max_pw_column()";
	// search in all sticks
	int use=0;
	for (int i = 0;i < this->nst;i++)
	{
		// Find out which column has not been used 
		if (!this->index[i])
		{
			// find out the only column which have most planewaves
			// in all the columns that have not been used ,
			// the column located at (max_i,max_j)
			// the planewave number is pw
			if (this->st_n[i] > pw)
			{
				max_i = this->st_i[i];
				max_j = this->st_j[i];
				pw    = this->st_n[i];
				use = i;
			}
		}
	}
	
	/*
	cout<<"\n use = "<<use
		<<" index = "<<index[use]
		<<" max_i = "<<max_i
		<<" max_j = "<<max_j
		<<" pw = "<<pw;
	*/
	// means the index stick has been used 
	this->index[use] = 1;

	//cout<<setw(24)<<"End max_pw_column()"<<endl;
	return;
}

void Parallel_PW::fft_dlay_set(void)
{
	int i = 0;

	int ip = 0;//index of processors
	delete[] st_start;
	this->st_start = new int[nproc_use];
	ZEROS(st_start, nproc_use);

	// calculate start stick index in each
	// processor in this pool.
	for (ip = 1;ip < nproc_use;ip++)
	{
		for (i = 0;i < ip;i++)
		{
			this->st_start[ip] += this->nst_per[i];
		}
//		cout<<"\n cpu="<<ip<<" st_start="<<this->st_start[ip];
	}

	int nxy = this->n1 * this->n2;

	delete[] isind;
	delete[] ismap;
	this->isind = new int[nxy];
	this->ismap = new int[this->nst];
	ZEROS(ismap, nst);

	for (i = 0;i < nxy;i++)
	{
		this->isind[i] = -1;
		//because isind[i]=0,means i point has stick,
		//so -1 has illegal meaning.
	}

//	cout << setw(12) << "isind dim" << setw(12) << nxy << endl;
//	cout << setw(12) << "ismap dim" << setw(12) << this->nst << endl;

	int *st_move = new int[nproc_use];
	ZEROS(st_move,nproc_use);
	
	for (i = 0;i < nxy;i++)
	{
		//	find out grid i belongs to which process
		//	index_ip[i] has been calculated before.
		ip = this->index_ip[i];

		if (ip >= 0)
		{
			//cout<<setw(12)<<i<<setw(12)<<ip<<setw(12)<<nst[ip]<<endl;
			int is = st_move[ip] + this->st_start[ip] ;
			this->ismap[ is ] = i;
			
			if (ip == rank_use)
			{
				//====================================
				// In plane x-y
				// each i correspond to one column
				// find out how many columns in this 
				// 'rank_use' process.
				//====================================
				this->isind[i] = st_move[ip];
			}
			// increase number of columns
			st_move[ip]++;	
		}
	}

	// Don't delete this test below,it's useful !
		/*
		cout<<"========================================"<<endl;
		cout<<setw(12)<<"col"<<setw(12)<<"nxy"<<endl;
		for(i=0;i<this->nst;i++)
		{
			cout<<setw(12)<<i<<setw(12)<<ismap[i]<<endl;
		}
		cout<<"========================================"<<endl;
		cout<<setw(12)<<"nxy"<<setw(12)<<"col(per)"<<endl;
		for(i=0;i<nxy;i++)
		{
			if(this->isind[i]>=0)
			{
				cout<<setw(12)<<i<<setw(12)<<isind[i]<<endl;
			}
		}
		*/	
	delete[] st_move;
	return;
}

void Parallel_PW::fft_map(
	int *ig2fft, // ig2fftc for charge grid, ig2fftw for wave function grid.
	const int ngm, // ngmc for charge grid, ngmw for wave function grid.
	const int &ngmc_g_in, // cutgg_num
        int ggchg_time
)
{
	if(test_pw) TITLE("Parallel_PW","fft_map");

	// if already found all the plane waves in this process, return.
	if(ngm_i == ngm) return;

	if(!allocate_igl2g)
	{
		delete[] ig_l2g;
		this->ig_l2g = new int[ngm];
		allocate_igl2g = true;
	}

	for (int ig = 0;ig < ngmc_g_in ;ig++)
	{
		int m1 = (int)pw.gdirect_global[ig].x;
		int m2 = (int)pw.gdirect_global[ig].y;

		if (m1 < 0){m1 += this->n1;}
		if (m2 < 0){m2 += this->n2;}

		int mc = m2 + m1 * this->n2; // mohan 2009-11-09

		if (this->isind[mc] >= 0)
		{
			if (pw.gg_global[ig] <= this->gcut)
			{
				this->ig_l2g[ngm_i] = ig;
				if(this->gcut == pw.ggchg) // for charge, not for wave functions.
				{
					// set 
					pw.gdirect[ngm_i] = pw.gdirect_global[ig];
					pw.gcar[ngm_i] = pw.gdirect[ngm_i]*ucell.G;
					pw.gg[ngm_i] = pw.gg_global[ig];
				}
				++ngm_i;
			}
		}
	}

        ngm_i_number[ggchg_time-1] = ngm_i; //LiuXh add 20180619

	/*
	if (ngm_i != this->npw_per[rank_use])
	{
		cout << setw(12) << "ngm_i = " << setw(12) << ngm_i << endl;
		cout << setw(12) << "npw_per = " << this->npw_per[rank_use] << endl;
		WARNING_QUIT("Parallel_PW::fft_map","ngm_i !=npw_per[rank_use]");
	}
	*/
	
//	cout<<setw(12)<<"i"<<setw(12)<<"x"<<setw(12)<<"y"<<setw(12)<<"z"<<endl;
	for (int ig=ngm_i_record; ig<ngm_i; ig++)
	{
		//  Mohan fix the bug 2008-4-3 17:19
		//	ig_l2g[i] is the index of ig_global
		//	it tells the position in ig_global
		//  for each plane wave in local processor
		//  The local plane wave number is ngm_i = pw.ngmc
		//  The total plane wave number is pw.ngmc_g
		int x = int(pw.gdirect_global[this->ig_l2g[ig]].x);
		int y = int(pw.gdirect_global[this->ig_l2g[ig]].y);
		int z = int(pw.gdirect_global[this->ig_l2g[ig]].z);

		if (x < 0) x += this->n1;
		if (y < 0) y += this->n2;
		if (z < 0) z += this->n3;

		const int index_now = y+x*n2;
		
	//	cout<<setw(12)<<x+y*n1<<setw(12)<<this->index_ip[x+y*this->n1]<<endl;
		if(this->isind[index_now]==-1)
		{
			cerr<<"I don't know this grid !"<<endl;
			cout<<setw(12)<<"xy"
				<<setw(12)<<"isind"<<setw(12)<<"ig2fft"<<endl;
			cout<<setw(12)<<index_now
				<<setw(12)<<this->isind[index_now]<<setw(12)<<ig2fft[ig]<<endl;
			WARNING_QUIT("Parallel_PW::fft_map","isind == -1 !");
		}
		ig2fft[ig] = z + this->isind[index_now] * this->n3;
		//cout << "\n i=" << i << " ig2fft=" << ig2fft[i];
	}

	ngm_i_record = ngm_i;
//	cout << "\n ngm_i_record = " << ngm_i_record;
	//cout << setw(12) << "ngm_i" << setw(12) << ngm_i << endl;

	return;
}


void Parallel_PW::print_data(ofstream &print)const
{
	print<<setw(12)<<"gcut"<<setw(12)<<this->gcut<<endl;
	print<<setw(12)<<"nst"<<setw(12)<<this->nst<<endl;
	print<<setw(12)<<"n1"<<setw(12)<<this->n1<<endl;
	print<<setw(12)<<"n2"<<setw(12)<<this->n2<<endl;
	print<<setw(12)<<"n3"<<setw(12)<<this->n3<<endl;

	print<<setw(12)<<"Processor"<<setw(12)<<"Planes"<<setw(12)<<"Columns"
		 <<setw(12)<<"st_start"<<setw(12)<<"PWs"<<endl;
	for(int i=0;i<nproc_use;i++)
	{
		print<<setw(12)<<i<<setw(12)<<this->npps[i]<<setw(12)<<this->nst_per[i]
			 <<setw(12)<<this->st_start[i]<<setw(12)<<this->npw_per[i]<<endl;
	}


	print<<"\n";
	print<<"====== ismap ======="<<endl;
	print<<setw(12)<<"Columns"<<setw(12)<<"nxy"<<endl;
	for(int i=0;i<this->nst;i++)
	{
		if(i>0)
		{
			if(this->ismap[i]<this->ismap[i-1])
			{
				print<<"--------------------------"<<endl;
			}
		}
		print<<setw(12)<<i<<setw(12)<<this->ismap[i]<<endl;
	}

	print<<"\n";
	print<<"====== isind ======"<<endl;
	print<<setw(12)<<"nxy"<<setw(12)<<"Columns"<<endl;
	for(int i=0;i<this->n1*this->n2;i++)
	{	
		if(this->isind[i]>=0)
		{
			print<<setw(12)<<i<<setw(12)<<this->isind[i]<<endl;
		}
	}

	return;
}

//LiuXh add a new function here,
//20180515
void Parallel_PW::fft_map_after_vc(
        int *ig2fft, // ig2fftc for charge grid, ig2fftw for wave function grid.
        const int ngm, // ngmc for charge grid, ngmw for wave function grid.
        const int &ngmc_g_in, // cutgg_num
        int ggchg_time
)
{
    if(test_pw) TITLE("Parallel_PW","fft_map_after_vc");

    if(!allocate_igl2g)
    {
        delete[] ig_l2g;
        this->ig_l2g = new int[ngm];
        allocate_igl2g = true;
    }

    if(this->gcut == pw.ggchg) ngm_i2 = ngm_i_number[ggchg_time-1];
    for (int ig = ngm_i_record2;ig < ngm_i2;ig++)
    {
        if(this->gcut == pw.ggchg)
        {
            pw.gcar[ig] = pw.gdirect[ig]*ucell.G;
            pw.gg[ig] = pw.get_NormG_cartesian(ig);
        }
    }


    if(this->gcut == pw.ggchg) ngm_i_record2 = ngm_i2;

    if(ggchg_time == 10) {ngm_i2 = 0; ngm_i_record2 = 0;}

    return;
}


void Parallel_PW::fft_map_final_scf(
	int *ig2fft, // ig2fftc for charge grid, ig2fftw for wave function grid.
	const int ngm, // ngmc for charge grid, ngmw for wave function grid.
	const int &ngmc_g_in // cutgg_num
)
{
	if(test_pw) TITLE("Parallel_PW","fft_map_final_scf");

	// if already found all the plane waves in this process, return.
	if(ngm_i_final_scf == ngm) return;

	if(!allocate_igl2g_final_scf)
	{
		delete[] ig_l2g;
		this->ig_l2g = new int[ngm];
		allocate_igl2g_final_scf = true;
	}

	for (int ig = 0;ig < ngmc_g_in ;ig++)
	{
		int m1 = (int)pw.gdirect_global[ig].x;
		int m2 = (int)pw.gdirect_global[ig].y;

		if (m1 < 0){m1 += this->n1;}
		if (m2 < 0){m2 += this->n2;}

		int mc = m2 + m1 * this->n2; // mohan 2009-11-09

		if (this->isind[mc] >= 0)
		{
			if (pw.gg_global[ig] <= this->gcut)
			{
				this->ig_l2g[ngm_i_final_scf] = ig;
				if(this->gcut == pw.ggchg) // for charge, not for wave functions.
				{
					// set 
					pw.gdirect[ngm_i_final_scf] = pw.gdirect_global[ig];
					pw.gcar[ngm_i_final_scf] = pw.gdirect[ngm_i_final_scf]*ucell.G;
					pw.gg[ngm_i_final_scf] = pw.gg_global[ig];
				}
				++ngm_i_final_scf;
			}
		}
	}

	for (int ig=ngm_i_record_final_scf; ig<ngm_i_final_scf; ig++)
	{
		//  Mohan fix the bug 2008-4-3 17:19
		//	ig_l2g[i] is the index of ig_global
		//	it tells the position in ig_global
		//  for each plane wave in local processor
		//  The local plane wave number is ngm_i = pw.ngmc
		//  The total plane wave number is pw.ngmc_g
		int x = int(pw.gdirect_global[this->ig_l2g[ig]].x);
		int y = int(pw.gdirect_global[this->ig_l2g[ig]].y);
		int z = int(pw.gdirect_global[this->ig_l2g[ig]].z);

		if (x < 0) x += this->n1;
		if (y < 0) y += this->n2;
		if (z < 0) z += this->n3;

		const int index_now = y+x*n2;
		
	//	cout<<setw(12)<<x+y*n1<<setw(12)<<this->index_ip[x+y*this->n1]<<endl;
		if(this->isind[index_now]==-1)
		{
			cerr<<"I don't know this grid !"<<endl;
			cout<<setw(12)<<"xy"
				<<setw(12)<<"isind"<<setw(12)<<"ig2fft"<<endl;
			cout<<setw(12)<<index_now
				<<setw(12)<<this->isind[index_now]<<setw(12)<<ig2fft[ig]<<endl;
			WARNING_QUIT("Parallel_PW::fft_map","isind == -1 !");
		}
		ig2fft[ig] = z + this->isind[index_now] * this->n3;
		//cout << "\n i=" << i << " ig2fft=" << ig2fft[i];
	}

	ngm_i_record_final_scf = ngm_i_final_scf;
//	cout << "\n ngm_i_record = " << ngm_i_record;
	//cout << setw(12) << "ngm_i" << setw(12) << ngm_i << endl;

	return;
}
