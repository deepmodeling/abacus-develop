#include "parallel_orbitals.h"
#include "../module_base/memory.h"
#include "module_orbital/ORB_control.h"
#ifdef __MPI
extern "C"
{
    #include "../module_base/blacs_connector.h"
	#include "../module_base/scalapack_connector.h"
}
#endif
Parallel_Orbitals::Parallel_Orbitals()
{
    loc_sizes = new int[1];
    trace_loc_row = new int[1];
    trace_loc_col = new int[1];

    testpb = 0;//mohan add 2011-03-16
	alloc_Z_LOC = false; //xiaohui add 2014-12-22
    // default value of nb is 1,
	// but can change to larger value from input.
    nb = 1;
	MatrixInfo.row_set = new int[1];
    MatrixInfo.col_set = new int[1];

    //in multi-k, 2D-block-division variables for FT (R<->k)
    nnr = 1;
    nlocdim = new int[1];	
	nlocstart = new int[1];

}

Parallel_Orbitals::~Parallel_Orbitals()
{
    delete[] trace_loc_row;
    delete[] trace_loc_col;
    delete[] loc_sizes;
    
    if (alloc_Z_LOC)//xiaohui add 2014-12-22
	{
		for(int is=0; is<this->nspin; is++)
		{
			delete[] Z_LOC[is];
		}
		delete[] Z_LOC;
	}
    delete[] MatrixInfo.row_set;
    delete[] MatrixInfo.col_set;
    
    delete[] nlocdim;
	delete[] nlocstart;
    
}

bool Parallel_Orbitals::in_this_processor(const int &iw1_all, const int &iw2_all) const
{
    if (trace_loc_row[iw1_all] == -1) return false;
    else if (trace_loc_col[iw2_all] == -1) return false;
    return true;
}

void ORB_control::set_trace(std::ofstream& ofs_running)
{
    ModuleBase::TITLE("ORB_control","set_trace");
    assert(nlocal > 0);
    
    Parallel_Orbitals* pv = &this->ParaV;

    delete[] pv->trace_loc_row;
    delete[] pv->trace_loc_col;

    ModuleBase::GlobalFunc::OUT(ofs_running,"trace_loc_row dimension",nlocal);
    ModuleBase::GlobalFunc::OUT(ofs_running,"trace_loc_col dimension",nlocal);

    pv->trace_loc_row = new int[nlocal];
    pv->trace_loc_col = new int[nlocal];
    // mohan update 2011-04-07
    for(int i=0; i<nlocal; i++)
    {
        pv->trace_loc_row[i] = -1;
        pv->trace_loc_col[i] = -1;
    }

    ModuleBase::Memory::record("ORB_control","trace_loc_row",nlocal,"int");
    ModuleBase::Memory::record("ORB_control","trace_loc_col",nlocal,"int");

    if(ks_solver=="lapack"
    || ks_solver=="cg"
    || ks_solver=="dav") //xiaohui add 2013-09-02
	{
		std::cout << " common settings for trace_loc_row and trace_loc_col " << std::endl;
		for (int i=0; i<nlocal; i++)
		{
			pv->trace_loc_row[i] = i;
			pv->trace_loc_col[i] = i;
		}
		pv->nrow = nlocal;
		pv->ncol = nlocal;
	}
#ifdef __MPI
    else if(ks_solver=="scalpack" || ks_solver=="genelpa" || ks_solver=="hpseps" 
		|| ks_solver=="selinv" || ks_solver=="scalapack_gvx") //xiaohui add 2013-09-02
    {
        // ofs_running << " nrow=" << nrow << std::endl;
        for (int irow=0; irow< pv->nrow; irow++)
        {
            int global_row = pv->MatrixInfo.row_set[irow];
            pv->trace_loc_row[global_row] = irow;
			// ofs_running << " global_row=" << global_row 
			// << " trace_loc_row=" << pv->trace_loc_row[global_row] << std::endl;
        }

        // ofs_running << " ncol=" << ncol << std::endl;
        for (int icol=0; icol< pv->ncol; icol++)
        {
            int global_col = pv->MatrixInfo.col_set[icol];
            pv->trace_loc_col[global_col] = icol;
			// ofs_running << " global_col=" << global_col 
			// << " trace_loc_col=" << pv->trace_loc_row[global_col] << std::endl;
        }
    }
#endif
    else 
    {
        std::cout << " Parallel Orbial, DIAGO_TYPE = " << ks_solver << std::endl;
        ModuleBase::WARNING_QUIT("ORB_control::set_trace","Check ks_solver.");
    }

    //---------------------------
    // print the trace for test.
    //---------------------------
    /*
    ofs_running << " " << std::setw(10) << "GlobalRow" << std::setw(10) << "LocalRow" << std::endl;
    for(int i=0; i<nlocal; i++)
    {
        ofs_running << " " << std::setw(10) << i << std::setw(10) << pv->trace_loc_row[i] << std::endl;

    }

    ofs_running << " " << std::setw(10) << "GlobalCol" << std::setw(10) << "LocalCol" << std::endl;
    for(int j=0; j<nlocal; j++)
    {
        ofs_running << " " << std::setw(10) << j << std::setw(10) << trace_loc_col[j] << std::endl;
    }
    */

    return;
}

#ifdef __MPI
inline int cart2blacs(
	MPI_Comm comm_2D,
	int nprows,
	int npcols,
    int Nlocal,
    int Nbands, 
    int nblk,
	int lld,
    int* desc,
    int* desc_wfc)
{
    int my_blacs_ctxt;
    int myprow, mypcol;
    int *usermap=new int[nprows*npcols];
    int info=0;
    for(int i=0; i<nprows; ++i)
    {
        for(int j=0; j<npcols; ++j)
        {
            int pcoord[2]={i, j};
            MPI_Cart_rank(comm_2D, pcoord, &usermap[i+j*nprows]);
        }
    }
    MPI_Fint comm_2D_f = MPI_Comm_c2f(comm_2D);
    Cblacs_get(comm_2D_f, 0, &my_blacs_ctxt);
    Cblacs_gridmap(&my_blacs_ctxt, usermap, nprows, nprows, npcols);
    Cblacs_gridinfo(my_blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);
    delete[] usermap;
    int ISRC=0;
    descinit_(desc, &Nlocal, &Nlocal, &nblk, &nblk, &ISRC, &ISRC, &my_blacs_ctxt, &lld, &info);
    descinit_(desc_wfc, &Nlocal, &Nbands, &nblk, &nblk, &ISRC, &ISRC, &my_blacs_ctxt, &lld, &info);

    return my_blacs_ctxt;
}
#endif

void ORB_control::divide_HS_2d
(
#ifdef __MPI
	MPI_Comm DIAG_WORLD,
#endif
    std::ofstream& ofs_running,
    std::ofstream& ofs_warning)
{
	ModuleBase::TITLE("ORB_control","divide_HS_2d");
	assert(nlocal>0);
    assert(dsize > 0);
    Parallel_Orbitals* pv = &this->ParaV;

#ifdef __MPI
	DIAG_HPSEPS_WORLD=DIAG_WORLD;
#endif

	if(dcolor!=0) return; // mohan add 2012-01-13

	// get the 2D index of computer.
	pv->dim0 = (int)sqrt((double)dsize); //mohan update 2012/01/13
	//while (GlobalV::NPROC_IN_POOL%dim0!=0)
	while (dsize%pv->dim0!=0)
	{
		pv->dim0 = pv->dim0 - 1;
	}
	assert(pv->dim0 > 0);
	pv->dim1=dsize/pv->dim0;

	if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"dim0",pv->dim0);
	if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"dim1",pv->dim1);

#ifdef __MPI
	// mohan add 2011-04-16
	if(nb2d==0)
	{
		if(nlocal>0) pv->nb = 1;
		if(nlocal>500) pv->nb = 32;
		if(nlocal>1000) pv->nb = 64;
	}
	else if(nb2d>0)
	{
		pv->nb = nb2d; // mohan add 2010-06-28
	}
	ModuleBase::GlobalFunc::OUT(ofs_running,"nb2d", pv->nb);

    this->set_parameters(ofs_running, ofs_warning);

	// call mpi_creat_cart
	this->mpi_creat_cart(&pv->comm_2D,pv->dim0,pv->dim1, ofs_running);

	// call mat_2d
    int try_nb = this->mat_2d(pv->comm_2D, nlocal, nbands, pv->nb,
        pv->MatrixInfo, ofs_running, ofs_warning);
    if(try_nb==1)
    {
        ofs_running<<" parameter nb2d is too large: nb2d = "<<pv->nb<<std::endl;
        ofs_running<<" reset nb2d to value 1, this set would make the program keep working but maybe get slower during diagonalization."<<std::endl;
        pv->nb = 1;
        try_nb = this->mat_2d(pv->comm_2D, nlocal, nbands, pv->nb,
        pv->MatrixInfo, ofs_running, ofs_warning);
    }

	// mohan add 2010-06-29
	pv->nrow = pv->MatrixInfo.row_num;
	pv->ncol = pv->MatrixInfo.col_num;
	pv->nloc = pv->MatrixInfo.col_num * pv->MatrixInfo.row_num;

	// init blacs context for genelpa
    if(ks_solver=="genelpa" || ks_solver=="scalapack_gvx")
    {
        pv->blacs_ctxt = cart2blacs(pv->comm_2D, pv->dim0, pv->dim1,
            nlocal, nbands, pv->nb, pv->nrow, pv->desc, pv->desc_wfc);
    }
#else // single processor used.
	pv->nb = nlocal;
	pv->nrow = nlocal;
	pv->ncol = nlocal;
	pv->nloc = nlocal * nlocal;
    this->set_parameters(ofs_running, ofs_warning);
    pv->MatrixInfo.row_b = 1;
	pv->MatrixInfo.row_num = nlocal;
	delete[] pv->MatrixInfo.row_set;
	pv->MatrixInfo.row_set = new int[nlocal];
	for(int i=0; i<nlocal; i++)
	{
		pv->MatrixInfo.row_set[i]=i;
	}
	pv->MatrixInfo.row_pos=0;

	pv->MatrixInfo.col_b = 1;
	pv->MatrixInfo.col_num = nlocal;
	delete[] pv->MatrixInfo.col_set;
	pv->MatrixInfo.col_set = new int[nlocal];
	for(int i=0; i<nlocal; i++)
	{
		pv->MatrixInfo.col_set[i]=i;
	}
	pv->MatrixInfo.col_pos=0;
#endif

	assert(pv->nloc>0);
	if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"MatrixInfo.row_num",pv->MatrixInfo.row_num);
	if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"MatrixInfo.col_num",pv->MatrixInfo.col_num);
	if(pv->testpb)ModuleBase::GlobalFunc::OUT(ofs_running,"nloc",pv->nloc);
	return;
}
