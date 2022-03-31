#include "ORB_control.h"
#include "ORB_gen_tables.h"
#include "../module_base/timer.h"
#include "../src_parallel/parallel_common.h"
#include "../src_io/wf_local.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/memory.h"

//#include "build_st_pw.h"

ORB_control::ORB_control(
    const bool& gamma_only_in,
    const int& nlocal_in,
    const int& nbands_in,
    const int& nspin_in,
    const int& dsize_in,
    const int& nb2d_in,
    const int& dcolor_in,
    const int& drank_in,
    const int& myrank_in,
    const std::string& calculation_in,
    const std::string& ks_solver_in) :
    gamma_only(gamma_only_in),
    nlocal(nlocal_in),
    nbands(nbands_in),
    nspin(nspin_in),
    dsize(dsize_in),
    nb2d(nb2d_in),
    dcolor(dcolor_in),
    drank(drank_in),
    myrank(myrank_in),
    calculation(calculation_in),
    ks_solver(ks_solver_in),
    setup_2d(true)
{
    this->ParaV.nspin = nspin_in;
}

ORB_control::ORB_control() :
    setup_2d(false)
{}
ORB_control::~ORB_control()
{}

void ORB_control::read_orb_first(
	std::ofstream &ofs_in, 
	LCAO_Orbitals &orb,
	const int &ntype, // mohan add 2021-04-26
	const int &lmax, // mohan add 2021-04-26 
	const double &lcao_ecut_in, // mohan add 2021-04-16
	const double &lcao_dk_in, // mohan add 2021-04-16
	const double &lcao_dr_in, // mohan add 2021-04-16
	const double &lcao_rmax_in, // mohan add 2021-04-16
	const bool &deepks_setorb,
	const int &out_mat_r,
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank // mohan add 2021-04-26
) 
{
    ModuleBase::TITLE("ORB_control","read_orb_first");
	ModuleBase::timer::tick("ORB_control","read_orb_first");
    
	/////////////////////////////////////////////////////////////////
	/// (1) FUNCTION : use 'info' to generate 'Numerical Orbital'
	///
	/// (1) RESULT : We have 'Numerical Orbital' for calculate S-table and T-table.
	/////////////////////////////////////////////////////////////////

	// mohan add 2021-04-16
	assert(ntype>0);
	assert(lmax>=0);
	assert(lcao_ecut_in>0.0);
	assert(lcao_dk_in>0.0);
	assert(lcao_dr_in>0.0);
	assert(lcao_rmax_in>0.0);

	// mohan add 2021-04-16
	orb.ecutwfc = lcao_ecut_in;
	orb.dk = lcao_dk_in;
	orb.dR = lcao_dr_in;
	orb.Rmax = lcao_rmax_in;
	
    orb.Read_Orbitals(
		ofs_in,
		ntype, 
		lmax, 
		deepks_setorb, 
		out_mat_r, 
		force_flag,
		my_rank);

	ModuleBase::timer::tick("ORB_control","read_orb_first");
	return;
}

void ORB_control::set_orb_tables(
	std::ofstream &ofs_in,
	ORB_gen_tables &OGT, 
	LCAO_Orbitals &orb,
	const double &lat0,
	const bool &deepks_setorb,
	const int &Lmax_exx,
	const int &nprojmax, 
	const int* nproj,
	const Numerical_Nonlocal* beta_) 
{
    ModuleBase::TITLE("ORB_control","set_orb_tables");
	ModuleBase::timer::tick("ORB_control","set_orb_tables");


#ifdef __NORMAL

#else
	if(calculation=="test")
	{
		ModuleBase::timer::tick("ORB_control","set_orb_tables");
		return;
	}
#endif


    ///////////////////////////////////////////////////////////////////
    /// (2) FUNCTION : Generate Gaunt_Coefficients and S-table using OGT.init
	/// 	   Must have 'Numerical Orbital' infomation
	///
	/// (2) RESULT : we have tabulated S table for use.
    ///////////////////////////////////////////////////////////////////
    const int job0 = 3;
    /// job0 :
    /// 1. generate overlap table
    /// 2. generate kinetic table
    /// 3. generate overlap & kinetic table
    OGT.gen_tables(ofs_in, job0, orb, Lmax_exx, deepks_setorb, nprojmax, nproj, beta_);
    // init lat0, in order to interpolated value from this table.

	assert(lat0>0.0);
    OGT.set_unit(lat0);

	ModuleBase::timer::tick("ORB_control","set_orb_tables");
    return;
}

void ORB_control::clear_after_ions(
	ORB_gen_tables &OGT, 
	LCAO_Orbitals &orb,
	const bool &deepks_setorb,
	const int* nproj_)
{
    ModuleBase::TITLE("ORB_control","clear_after_ions");
    OGT.MOT.Destroy_Table(orb);
    OGT.tbeta.Destroy_Table_Beta(orb.get_ntype(), orb.Phi, nproj_);
    
	//caoyu add 2021-03-18
    if (deepks_setorb) 
	{
        OGT.talpha.Destroy_Table_Alpha(orb);
    }
    return;
}


void ORB_control::setup_2d_division(std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("ORB_control","setup_2d_division");
    ofs_running << "\n SETUP THE DIVISION OF H/S MATRIX" << std::endl;
    
    // (1) calculate nrow, ncol, nloc.
    if (ks_solver=="genelpa" || ks_solver=="hpseps" || ks_solver=="scalpack" 
        || ks_solver=="selinv" || ks_solver=="scalapack_gvx")
    {
        ofs_running << " divide the H&S matrix using 2D block algorithms." << std::endl;
#ifdef __MPI
        // storage form of H and S matrices on each processor
        // is determined in 'divide_HS_2d' subroutine
        this->divide_HS_2d(DIAG_WORLD, ofs_running, ofs_warning);
#else
        ModuleBase::WARNING_QUIT("LCAO_Matrix::init","diago method is not ready.");
#endif
	}
	else
	{
		// the full matrix
		this->ParaV.nloc = nlocal * nlocal;
	}

	// (2) set the trace, then we can calculate the nnr.
	// for 2d: calculate po.nloc first, then trace_loc_row and trace_loc_col
	// for O(N): calculate the three together.
	this->set_trace(ofs_running);
}


void ORB_control::set_parameters(std::ofstream& ofs_running,
        std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("ORB_control","set_parameters");

    Parallel_Orbitals* pv = &this->ParaV;
    // set loc_size
	if(gamma_only)//xiaohui add 2014-12-21
	{
		pv->loc_size=nbands/dsize;

		// mohan add 2012-03-29
		if(pv->loc_size==0)
		{
			ofs_warning << " loc_size=0" << " in proc " << myrank+1 << std::endl;
			ModuleBase::WARNING_QUIT("ORB_control::set_parameters","NLOCAL < DSIZE");
		}

		if (drank<nbands%dsize) pv->loc_size+=1;
		if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"local size",pv->loc_size);

		// set loc_sizes
		delete[] pv->loc_sizes;
		pv->loc_sizes = new int[dsize];
		ModuleBase::GlobalFunc::ZEROS(pv->loc_sizes, dsize);

		pv->lastband_in_proc = 0;
		pv->lastband_number = 0;
		int count_bands = 0;
		for (int i=0; i<dsize; i++)
		{
			if (i<nbands%dsize)
			{
				// mohan modify 2010-07-05
				pv->loc_sizes[i]=nbands/dsize+1;
			}
			else
			{
				pv->loc_sizes[i]=nbands/dsize;
			}
			count_bands += pv->loc_sizes[i];
			if (count_bands >= nbands)
			{
				pv->lastband_in_proc = i;
				pv->lastband_number = nbands - (count_bands - pv->loc_sizes[i]);
				break;
			}
		}
	}
	else
	{
		pv->loc_size=nlocal/dsize;

		// mohan add 2012-03-29
		if(pv->loc_size==0)
		{
			ofs_warning << " loc_size=0" << " in proc " << myrank+1 << std::endl;
			ModuleBase::WARNING_QUIT("ORB_control::set_parameters","NLOCAL < DSIZE");
		}

		if (drank<nlocal%dsize) 
		{
            pv->loc_size += 1;
        }
		if(pv->testpb) ModuleBase::GlobalFunc::OUT(ofs_running,"local size",pv->loc_size);

		// set loc_sizes
		delete[] pv->loc_sizes;
		pv->loc_sizes = new int[dsize];
		ModuleBase::GlobalFunc::ZEROS(pv->loc_sizes, dsize);

		pv->lastband_in_proc = 0;
		pv->lastband_number = 0;
		int count_bands = 0;
		for (int i=0; i<dsize; i++)
		{
			if (i<nlocal%dsize)
			{
				// mohan modify 2010-07-05
				pv->loc_sizes[i]=nlocal/dsize+1;
			}
			else
			{
				pv->loc_sizes[i]=nlocal/dsize;
			}
			count_bands += pv->loc_sizes[i];
			if (count_bands >= nbands)
			{
				pv->lastband_in_proc = i;
				pv->lastband_number = nbands - (count_bands - pv->loc_sizes[i]);
				break;
			}
		}
	}//xiaohui add 2014-12-21

    if (ks_solver=="hpseps") //LiuXh add 2021-09-06, clear memory, Z_LOC only used in hpseps solver
    {
	    pv->Z_LOC = new double*[nspin];
	    for(int is=0; is<nspin; is++)
	    {
		    pv->Z_LOC[is] = new double[pv->loc_size * nlocal];
		    ModuleBase::GlobalFunc::ZEROS(pv->Z_LOC[is], pv->loc_size * nlocal);
	    }
	    pv->alloc_Z_LOC = true;//xiaohui add 2014-12-22
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"lastband_in_proc", pv->lastband_in_proc);
    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"lastband_number", pv->lastband_number);

    return;
}


#ifdef __MPI
// creat the 'comm_2D' stratege.
void ORB_control::mpi_creat_cart(MPI_Comm* comm_2D,
    int prow, int pcol, std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ORB_control","mpi_creat_cart");
    // the matrix is divided as ( dim[0] * dim[1] )
    int dim[2];
    int period[2]={1,1};
    int reorder=0;
    dim[0]=prow;
    dim[1]=pcol;

    if(this->ParaV.testpb) ofs_running << " dim = " << dim[0] << " * " << dim[1] << std::endl;

    MPI_Cart_create(DIAG_WORLD,2,dim,period,reorder,comm_2D);
    return;
}
#endif

#ifdef __MPI
void ORB_control::mat_2d(MPI_Comm vu,
    const int &M_A,
    const int &N_A,
    const int &nb,
    LocalMatrix &LM,
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
    ModuleBase::TITLE("ORB_control", "mat_2d");
    
    Parallel_Orbitals* pv = &this->ParaV;
    
    int dim[2];
    int period[2];
    int coord[2];
    int i,j,k,end_id;
    int block;

    // (0) every processor get it's id on the 2D comm
    // : ( coord[0], coord[1] )
    MPI_Cart_get(vu,2,dim,period,coord);

    // (1.1) how many blocks at least
    // eg. M_A = 6400, nb = 64;
    // so block = 10;
    block=M_A/nb;

    // (1.2) If data remain, add 1.
    if (block*nb<M_A)
    {
        block++;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"Total Row Blocks Number",block);

	// mohan add 2010-09-12
	if(dim[0]>block)
	{
		ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
		ofs_warning << " but, the number of row blocks is " << block << std::endl;
		ModuleBase::WARNING_QUIT("ORB_control::mat_2d","some processor has no row blocks, try a smaller 'nb2d' parameter.");
	}

    // (2.1) row_b : how many blocks for this processor. (at least)
    LM.row_b=block/dim[0];

    // (2.2) row_b : how many blocks in this processor.
    // if there are blocks remain, some processors add 1.
    if (coord[0]<block%dim[0])
    {
        LM.row_b++;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"Local Row Block Number",LM.row_b);

    // (3) end_id indicates the last block belong to
    // which processor.
    if (block%dim[0]==0)
    {
        end_id=dim[0]-1;
    }
    else
    {
        end_id=block%dim[0]-1;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"Ending Row Block in processor",end_id);

    // (4) row_num : how many rows in this processors :
    // the one owns the last block is different.
    if (coord[0]==end_id)
    {
        LM.row_num=(LM.row_b-1)*nb+(M_A-(block-1)*nb);
    }
    else
    {
        LM.row_num=LM.row_b*nb;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"Local rows (including nb)",LM.row_num);

    // (5) row_set, it's a global index :
    // save explicitly : every row in this processor
    // belongs to which row in the global matrix.
    delete[] LM.row_set;
    LM.row_set= new int[LM.row_num];
    j=0;
    for (i=0; i<LM.row_b; i++)
    {
        for (k=0; k<nb&&(coord[0]*nb+i*nb*dim[0]+k<M_A); k++,j++)
        {
            LM.row_set[j]=coord[0]*nb+i*nb*dim[0]+k;
           // ofs_running << " j=" << j << " row_set=" << LM.row_set[j] << std::endl;
        }
    }

    // the same procedures for columns.
    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"Total Col Blocks Number",block);

	if(dim[1]>block)
	{
		ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
		ofs_warning << " but, the number of column blocks is " << block << std::endl;
		ModuleBase::WARNING_QUIT("ORB_control::mat_2d","some processor has no column blocks.");
	}

    LM.col_b=block/dim[1];
    if (coord[1]<block%dim[1])
    {
        LM.col_b++;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"Local Row Block Number",LM.col_b);

    if (block%dim[1]==0)
    {
        end_id=dim[1]-1;
    }
    else
    {
        end_id=block%dim[1]-1;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"Ending Row Block in processor",end_id);

    if (coord[1]==end_id)
    {
        LM.col_num=(LM.col_b-1)*nb+(M_A-(block-1)*nb);
    }
    else
    {
        LM.col_num=LM.col_b*nb;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"Local columns (including nb)",LM.row_num);

    delete[] LM.col_set;
    LM.col_set = new int[LM.col_num];

    j=0;
    for (i=0; i<LM.col_b; i++)
    {
        for (k=0; k<nb&&(coord[1]*nb+i*nb*dim[1]+k<M_A); k++,j++)
        {
            LM.col_set[j]=coord[1]*nb+i*nb*dim[1]+k;
        }
    }
    LM.col_pos=0;
    LM.row_pos = 0;

    // for wavefuncton , calculate nbands_loc
    block=N_A/nb;
    if (block*nb<N_A)
    {
        block++;
    }
    if(dim[1]>block)
	{
		ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
		ofs_warning << " but, the number of bands-row-block is " << block << std::endl;
		ModuleBase::WARNING_QUIT("ORB_control::mat_2d","some processor has no bands-row-blocks.");
    }
    int col_b_bands = block / dim[1];
    if (coord[1] < block % dim[1])
    {
        col_b_bands++;
    }
    if (block%dim[1]==0)
    {
        end_id=dim[1]-1;
    }
    else
    {
        end_id=block%dim[1]-1;
    }
    if (coord[1]==end_id)
    {
        pv->ncol_bands=(col_b_bands-1)*nb+(N_A-(block-1)*nb);
    }
    else
    {
        pv->ncol_bands=col_b_bands*nb;
    }
    pv->nloc_wfc = pv->ncol_bands * LM.row_num;

    return;
}
#endif


#ifdef __MPI
// A : contains total matrix element in processor.
void ORB_control::data_distribution(
    MPI_Comm comm_2D,
    const std::string &file,
    const int &n,
    const int &nb,
    double *A,
    const LocalMatrix &LM)
{
    ModuleBase::TITLE("ORB_control", "data_distribution");
    Parallel_Orbitals* pv = &this->ParaV;
    MPI_Comm comm_row;
    MPI_Comm comm_col;
    MPI_Status status;

    int dim[2];
    int period[2];
    int coord[2];
    MPI_Cart_get(comm_2D,2,dim,period,coord);

    if(pv->testpb) ofs_running << "\n dim = " << dim[0] << " * " << dim[1] << std::endl;
    if(pv->testpb) ofs_running << " coord = ( " << coord[0] << " , " << coord[1] << ")." << std::endl;
    if(pv->testpb) ofs_running << " n = " << n << std::endl;

    mpi_sub_col(comm_2D,&comm_col);
    mpi_sub_row(comm_2D,&comm_row);

    // total number of processors
    const int myid = coord[0]*dim[1]+coord[1];

    // the matrix is n * n
    double* ele_val = new double[n];
    double* val = new double[n];
    int* sends = new int[dim[1]];
    int* fpt = new int[dim[1]];
    int* snd = new int[dim[1]];
    int* temp = new int[dim[1]];

    ModuleBase::GlobalFunc::ZEROS(ele_val, n);
    ModuleBase::GlobalFunc::ZEROS(val, n);
    ModuleBase::GlobalFunc::ZEROS(sends, dim[1]);
    ModuleBase::GlobalFunc::ZEROS(fpt, dim[1]);
    ModuleBase::GlobalFunc::ZEROS(snd, dim[1]);
    ModuleBase::GlobalFunc::ZEROS(temp, dim[1]);

    // the columes of matrix is divided by 'dim[1]' 'rows of processors'.
    // collect all information of each 'rows of processors'
    // collection data is saved in 'sends'
    snd[coord[1]] = LM.col_num;
    MPI_Allgather(&snd[coord[1]],1,MPI_INT,sends,1,MPI_INT,comm_row);

    // fpt : start column index after applied 'mat_2d' reorder algorithms
    // to matrix.
    fpt[0] = 0;
    for (int i=1; i<dim[1]; i++)
    {
        fpt[i]=fpt[i-1]+sends[i-1];
//      ofs_running << " col_pro = " << i << " start_col = " << fpt[i] << std::endl;
    }

//    ofs_running << "\n myid = " << myid << std::endl;

    int cur_i = 0;

    int iacol;
    int iarow;
    int ai;
    int aj;
    int tag = 0;

    bool find = true;
    if (myid==0)
    {
        FILE *fp;
        fp=fopen(file.c_str(),"rb");
        if (fp==NULL)
        {
            std::cout << " Can't find file : " << file << std::endl;
            find = false;
        }
        else
        {
            ofs_running << " Open file : " << file << std::endl;
            int dim = 0;
            fread(&dim,sizeof(int),1,fp);
            if (dim!=n)
            {
                find = false;
            }
            ofs_running << " Read in dimension = " << dim << std::endl;
        }
        int nrow = 0;
        while (nrow<n && !feof(fp))
        {
            ModuleBase::GlobalFunc::ZEROS(ele_val, n);
            ModuleBase::GlobalFunc::ZEROS(val, n);

            // read om one row elements.
//            ofs_running << "\n nrow = " << nrow << std::endl;

            for (int i=nrow; i<n; i++)
            {
                //if ((i-nrow)%8==0)ofs_running << std::endl;
                fread(&ele_val[i],sizeof(double),1,fp);
                //			ofs_running << " " << ele_val[i];
            }

            // start position of col_pro.
            for (int i=0; i<dim[1]; i++)
            {
                temp[i] = fpt[i];
            }

            for (int k=0; k<n; k++)
            {
                // calculate iarow and iacol.
                // belong to which col_pro.
                indxg2p(comm_2D,nb,nrow,k,&iarow,&iacol);
                val[temp[iacol]]=ele_val[k];
                temp[iacol]++;
            }

            indxg2l(nrow,0,nb,dim[0],dim[1],&ai,&aj);
            indxg2p(comm_2D,nb,nrow,0,&iarow,&iacol);

            const int incx = 1;
            if (iarow==0&&iacol==0)
            {
                BlasConnector::copy(LM.col_num,val,incx,&A[ai*LM.col_num],incx);
                for (int i=1; i<dim[1]; i++)
                {
//					ofs_running << " send to processor " << iarow*dim[1]+i << std::endl;
                    MPI_Send(&val[fpt[i]],sends[i],MPI_DOUBLE,iarow*dim[1]+i,tag,DIAG_WORLD);
                }
            }
            else
            {
                for (int i=0; i<dim[1]; i++)
                {
//					ofs_running << " else, send to processor " << iarow*dim[1]+i << std::endl;
                    MPI_Send(&val[fpt[i]],sends[i],MPI_DOUBLE,iarow*dim[1]+i,tag,DIAG_WORLD);
                }
            }
            nrow++;
        }// end read in nrow

        fclose(fp);
    }
    else
    {
        for (int j=0; j<LM.row_num; j++)
        {
//			ofs_running << " receive row = " <<  j << std::endl;
            MPI_Recv(&A[j*LM.col_num],LM.col_num,MPI_DOUBLE,0,tag,DIAG_WORLD,&status);
        }
    }

    /*
    for (int i=0; i<LM.row_num; i++)
    {
        ofs_running << "\n\n Row = " << i << std::endl;
        for (int j=0; j<LM.col_num; j++)
        {
            if (j%8==0) ofs_running << std::endl;
            ofs_running << " " << A[j*LM.col_num+i];
        }
    }
    */

    MPI_Comm_free(&comm_row);
    MPI_Comm_free(&comm_col);

    delete[] ele_val;
    delete[] val;
    delete[] sends;
    delete[] fpt;
    delete[] snd;
    delete[] temp;

#ifdef __MPI
    Parallel_Common::bcast_bool(find);
#endif

    //ModuleBase::GlobalFunc::OUT(ofs_running,"Find the H/S file",find);

    if (!find)
    {
        ModuleBase::WARNING_QUIT("ORB_control::data_distribution","Can't find the H/S file");
    }

    return;
}
#endif

/* this function is not called
#ifdef __MPI
#include "src_pdiag/pdgseps.h"
#include "src_pdiag/pzgseps.h"
void ORB_control::readin(
    const std::string& fa,
	const std::string &fb,
	const int &nlocal_tot,
	double *eigen,
	double *eigvr)
{
    ModuleBase::TITLE("Pdiag_Double","readin");

    Parallel_Orbitals* pv = &this->ParaV;
    int coord[2];
    int dim[2];
    int period[2];
    int i,j,tmp1,tmp2;
    int k,loc_size,loc_pos;
    double time1,time2;

    MPI_Comm comm=DIAG_HPSEPS_WORLD,comm_2D,comm_col,comm_row,newcomm;

    dim[0]=(int)sqrt((double)dsize);

    while (dsize%dim[0]!=0)
    {
        dim[0]=dim[0]-1;
    }
    dim[1]=dsize/dim[0];

    // call mpi_creat_cart
    this->mpi_creat_cart(&pv->comm_2D,dim[0],dim[1]);

    // call mat_2d
    this->mat_2d(pv->comm_2D, nlocal_tot,nlocal_tot,pv->nb,pv->MatrixInfo);

    pv->loc_size=nlocal_tot/dsize;
    if (drank<nlocal_tot%dsize) loc_size=loc_size+1;

    ofs_running << " loc_size = " << loc_size;

    //Distribute the matrix
    const long nloc = pv->MatrixInfo.col_num * pv->MatrixInfo.row_num;

    double *A = new double[nloc];
    double *B = new double[nloc];
    double *Z = new double[loc_size*nlocal_tot];
    ModuleBase::GlobalFunc::ZEROS(A, nloc);
    ModuleBase::GlobalFunc::ZEROS(B, nloc);
    ModuleBase::GlobalFunc::ZEROS(Z, loc_size * nlocal_tot);

    ofs_running << "\n Data distribution of H." << std::endl;
    this->data_distribution(pv->comm_2D,fa,nlocal_tot,pv->nb,A,pv->MatrixInfo);
    ofs_running << "\n Data distribution of S." << std::endl;
    this->data_distribution(pv->comm_2D,fb,nlocal_tot,pv->nb,B,pv->MatrixInfo);

    time1=MPI_Wtime();
    // call pdgseps
    char uplo = 'U';
    pdgseps(pv->comm_2D,nlocal_tot,pv->nb,A,B,Z,eigen,pv->MatrixInfo,uplo,pv->loc_size,loc_pos);
    time2=MPI_Wtime();
    ModuleBase::GlobalFunc::OUT(ofs_running,"time1",time1);
    ModuleBase::GlobalFunc::OUT(ofs_running,"time2",time2);

    //this->gath_eig(comm,n,eigvr,Z);

    ofs_running << "\n " << std::setw(6) << "Band" << std::setw(25) << "Ry" << std::setw(25) << " eV" << std::endl;
    for(int i=0; i<nlocal_tot; i++)
    {
        ofs_running << " " << std::setw(6) << i << std::setw(25) << eigen[i] << std::setw(25)<< eigen[i] * 13.6058 << std::endl;
    }

    delete[] A;
    delete[] B;
    delete[] Z;
}
#endif
*/