#include "pdiag_basic.h"
#include "src_parallel/parallel_common.h"
#include "src_pw/global.h"
#include "src_io/wf_local.h"
#include "src_global/lapack_connector.h"

Pdiag_Basic::Pdiag_Basic()
{
    loc_sizes = new int[1];
	testpb = 0;//mohan add 2011-03-16
	alloc_Z_LOC = false; //xiaohui add 2014-12-22
}

Pdiag_Basic::~Pdiag_Basic()
{
    delete[] loc_sizes;
	if(alloc_Z_LOC)//xiaohui add 2014-12-22
	{
		for(int is=0; is<NSPIN; is++)
		{
			delete[] Z_LOC[is];
		}
		delete[] Z_LOC;
	}
}

void Pdiag_Basic::set_parameters(void)
{
    TITLE("Pdiag_Basic","set_parameters");

    // set loc_size
	if(GAMMA_ONLY_LOCAL)//xiaohui add 2014-12-21
	{
		loc_size=NBANDS/DSIZE;

		// mohan add 2012-03-29
		if(loc_size==0)
		{
			ofs_warning << " loc_size=0" << " in proc " << MY_RANK+1 << endl;
			WARNING_QUIT("Pdiag_Basic::set_parameters","NLOCAL < DSIZE");
		}

		if (DRANK<NBANDS%DSIZE) loc_size=loc_size+1;
		if(testpb)OUT(ofs_running,"local size",loc_size);

		// set loc_sizes
		delete[] loc_sizes;
		loc_sizes = new int[DSIZE];
		ZEROS(loc_sizes, DSIZE);

		this->lastband_in_proc = 0;
		this->lastband_number = 0;
		int count_bands = 0;
		for (int i=0; i<DSIZE; i++)
		{
			if (i<NBANDS%DSIZE)
			{
				// mohan modify 2010-07-05
				loc_sizes[i]=NBANDS/DSIZE+1;
			}
			else
			{
				loc_sizes[i]=NBANDS/DSIZE;
			}
			count_bands += loc_sizes[i];
			if (count_bands >= NBANDS)
			{
				lastband_in_proc = i;
				lastband_number = NBANDS - (count_bands - loc_sizes[i]);
				break;
			}
		}
	}
	else
	{
		loc_size=NLOCAL/DSIZE;

		// mohan add 2012-03-29
		if(loc_size==0)
		{
			ofs_warning << " loc_size=0" << " in proc " << MY_RANK+1 << endl;
			WARNING_QUIT("Pdiag_Basic::set_parameters","NLOCAL < DSIZE");
		}

		if (DRANK<NLOCAL%DSIZE) 
		{
			loc_size=loc_size+1;
		}
		if(testpb) OUT(ofs_running,"local size",loc_size);

		// set loc_sizes
		delete[] loc_sizes;
		loc_sizes = new int[DSIZE];
		ZEROS(loc_sizes, DSIZE);

		this->lastband_in_proc = 0;
		this->lastband_number = 0;
		int count_bands = 0;
		for (int i=0; i<DSIZE; i++)
		{
			if (i<NLOCAL%DSIZE)
			{
				// mohan modify 2010-07-05
				loc_sizes[i]=NLOCAL/DSIZE+1;
			}
			else
			{
				loc_sizes[i]=NLOCAL/DSIZE;
			}
			count_bands += loc_sizes[i];
			if (count_bands >= NBANDS)
			{
				lastband_in_proc = i;
				lastband_number = NBANDS - (count_bands - loc_sizes[i]);
				break;
			}
		}
	}//xiaohui add 2014-12-21

	Z_LOC = new double*[NSPIN];
	for(int is=0; is<NSPIN; is++)
	{
		Z_LOC[is] = new double[loc_size * NLOCAL];
		ZEROS(Z_LOC[is], loc_size * NLOCAL);
	}
	alloc_Z_LOC = true;//xiaohui add 2014-12-22

    if(testpb)OUT(ofs_running,"lastband_in_proc",lastband_in_proc);
    if(testpb)OUT(ofs_running,"lastband_number",lastband_number);

    return;
}


#ifdef __MPI
// creat the 'comm_2D' stratege.
void Pdiag_Basic::mpi_creat_cart(MPI_Comm *comm_2D, int prow, int pcol)
{
    TITLE("Pdiag_Basic","mpi_creat_cart");
    // the matrix is divided as ( dim[0] * dim[1] )
    int dim[2];
    int period[2]={1,1};
    int reorder=0;
    dim[0]=prow;
    dim[1]=pcol;

    if(testpb)ofs_running << " dim = " << dim[0] << " * " << dim[1] << endl;

    MPI_Cart_create(DIAG_WORLD,2,dim,period,reorder,comm_2D);
    return;
}
#endif

#ifdef __MPI
void Pdiag_Basic::mat_2d(MPI_Comm vu,
                         const int &M_A,
                         const int &N_A,
                         const int &nb,
                         LocalMatrix &LM)
{
    TITLE("Pdiag_Basic","mat_2d");
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

    if(testpb)OUT(ofs_running,"Total Row Blocks Number",block);

	// mohan add 2010-09-12
	if(dim[0]>block)
	{
		ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << endl;
		ofs_warning << " but, the number of row blocks is " << block << endl;
		WARNING_QUIT("Pdiag_Basic::mat_2d","some processor has no row blocks, try a smaller 'nb2d' parameter.");
	}

    // (2.1) row_b : how many blocks for this processor. (at least)
    LM.row_b=block/dim[0];

    // (2.2) row_b : how many blocks in thie processor.
    // if there are blocks remain, some processors add 1.
    if (coord[0]<block%dim[0])
    {
        LM.row_b++;
    }

    if(testpb)OUT(ofs_running,"Local Row Block Number",LM.row_b);

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

    if(testpb)OUT(ofs_running,"Ending Row Block in processor",end_id);

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

    if(testpb)OUT(ofs_running,"Local rows (including nb)",LM.row_num);

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
           // ofs_running << " j=" << j << " row_set=" << LM.row_set[j] << endl;
        }
    }

    // the same procedures for columns.
    block=N_A/nb;
    if (block*nb<N_A)
    {
        block++;
    }

    if(testpb)OUT(ofs_running,"Total Col Blocks Number",block);

	if(dim[1]>block)
	{
		ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << endl;
		ofs_warning << " but, the number of column blocks is " << block << endl;
		WARNING_QUIT("Pdiag_Basic::mat_2d","some processor has no column blocks.");
	}

    LM.col_b=block/dim[1];
    if (coord[1]<block%dim[1])
    {
        LM.col_b++;
    }

    if(testpb)OUT(ofs_running,"Local Row Block Number",LM.col_b);

    if (block%dim[1]==0)
    {
        end_id=dim[1]-1;
    }
    else
    {
        end_id=block%dim[1]-1;
    }

    if(testpb)OUT(ofs_running,"Ending Row Block in processor",end_id);

    if (coord[1]==end_id)
    {
        LM.col_num=(LM.col_b-1)*nb+(N_A-(block-1)*nb);
    }
    else
    {
        LM.col_num=LM.col_b*nb;
    }

    if(testpb)OUT(ofs_running,"Local columns (including nb)",LM.row_num);

    delete[] LM.col_set;
    LM.col_set = new int[LM.col_num];

    j=0;
    for (i=0; i<LM.col_b; i++)
    {
        for (k=0; k<nb&&(coord[1]*nb+i*nb*dim[1]+k<N_A); k++,j++)
        {
            LM.col_set[j]=coord[1]*nb+i*nb*dim[1]+k;
        }
    }
    LM.col_pos=0;
    LM.row_pos=0;
    return;
}
#endif


#ifdef __MPI
// A : contains total matrix element in processor.
void Pdiag_Basic::data_distribution(
    MPI_Comm comm_2D,
    const string &file,
    const int &n,
    const int &nb,
    double *A,
    const LocalMatrix &LM)
{
    TITLE("Pdiag_Basic","data_distribution");
    MPI_Comm comm_row;
    MPI_Comm comm_col;
    MPI_Status status;

    int dim[2];
    int period[2];
    int coord[2];
    MPI_Cart_get(comm_2D,2,dim,period,coord);

    if(testpb) ofs_running << "\n dim = " << dim[0] << " * " << dim[1] << endl;
    if(testpb) ofs_running << " coord = ( " << coord[0] << " , " << coord[1] << ")." << endl;
    if(testpb) ofs_running << " n = " << n << endl;

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

    ZEROS(ele_val, n);
    ZEROS(val, n);
    ZEROS(sends, dim[1]);
    ZEROS(fpt, dim[1]);
    ZEROS(snd, dim[1]);
    ZEROS(temp, dim[1]);

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
//      ofs_running << " col_pro = " << i << " start_col = " << fpt[i] << endl;
    }

//    ofs_running << "\n myid = " << myid << endl;

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
            cout << " Can't find file : " << file << endl;
            find = false;
        }
        else
        {
            ofs_running << " Open file : " << file << endl;
            int dim = 0;
            fread(&dim,sizeof(int),1,fp);
            if (dim!=n)
            {
                find = false;
            }
            ofs_running << " Read in dimension = " << dim << endl;
        }
        int nrow = 0;
        while (nrow<n && !feof(fp))
        {
            ZEROS(ele_val, n);
            ZEROS(val, n);

            // read om one row elements.
//            ofs_running << "\n nrow = " << nrow << endl;

            for (int i=nrow; i<n; i++)
            {
                //if ((i-nrow)%8==0)ofs_running << endl;
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
                LapackConnector::copy(LM.col_num,val,incx,&A[ai*LM.col_num],incx);
                for (int i=1; i<dim[1]; i++)
                {
//					ofs_running << " send to processor " << iarow*dim[1]+i << endl;
                    MPI_Send(&val[fpt[i]],sends[i],MPI_DOUBLE,iarow*dim[1]+i,tag,DIAG_WORLD);
                }
            }
            else
            {
                for (int i=0; i<dim[1]; i++)
                {
//					ofs_running << " else, send to processor " << iarow*dim[1]+i << endl;
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
//			ofs_running << " receive row = " <<  j << endl;
            MPI_Recv(&A[j*LM.col_num],LM.col_num,MPI_DOUBLE,0,tag,DIAG_WORLD,&status);
        }
    }

    /*
    for (int i=0; i<LM.row_num; i++)
    {
        ofs_running << "\n\n Row = " << i << endl;
        for (int j=0; j<LM.col_num; j++)
        {
            if (j%8==0) ofs_running << endl;
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

    //OUT(ofs_running,"Find the H/S file",find);

    if (!find)
    {
        WARNING_QUIT("Pdiag_Basic::data_distribution","Can't find the H/S file");
    }

    return;
}
#endif

#ifdef __MPI
#include "src_pw/occupy.h"
void Pdiag_Basic::gath_eig(MPI_Comm comm,int n,double **wfc,double *Z)
{
    TITLE("Pdiag_Basic","gath_eig");
    time_t time_start = time(NULL);
//  ofs_running << " Start gath_eig Time : " << ctime(&time_start);

    int i,j,k;
    int nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);

    double **ctot;

	// mohan add 2010-07-03
	// the occupied bands are useless
	// for calculating charge density.
	if(DRANK>lastband_in_proc)
	{
		delete[] Z;
	}

	// first we need to collect all
	// the occupied bands.
	// NBANDS * NLOCAL	
	if(myid==0)
	{
		ctot = new double*[NBANDS];
    	for (int i=0; i<NBANDS; i++)
    	{
        	ctot[i] = new double[NLOCAL];
        	ZEROS(ctot[i], NLOCAL);
    	}
    	Memory::record("Pdiag_Basic","ctot",NBANDS*NLOCAL,"double");
	}

    k=0;
    if (myid==0)
    {
        // mohan add nbnd0 2010-07-02
        int nbnd0 = -1;
        if (NBANDS < loc_sizes[0])
        {
			// means all bands in this processor
			// is needed ( is occupied)
            nbnd0 = NBANDS;
        }
        else
        {
			// means this processor only save
			// part of NBANDS.
            nbnd0 = loc_sizes[0];
        }
        if(testpb) ofs_running << " nbnd in processor 0 is " << nbnd0 << endl;

//printf("from 0 to %d\n",nbnd0-1);
        for (i=0; i<nbnd0; i++)
        {
            for (j=0; j<NLOCAL; j++)
            {
				// change the order in processor 0.
				// the contribution from processor 0.
                ctot[k][j]=Z[j*loc_sizes[0]+i];
            }
            k++;
        }
		// Z is useless in processor 0 now.
		delete[] Z;
    }
    MPI_Barrier(comm);


    for (i=1; i<= this->lastband_in_proc; i++)
    {
        // mohan fix bug 2010-07-02
        // rows indicates the data structure of Z.
        // mpi_times indicates the data distribution
        // time, each time send a band.
        int rows = loc_sizes[i];
        int mpi_times;
        if (i==lastband_in_proc)
        {
            mpi_times = lastband_number;
        }
        else
        {
            mpi_times = loc_sizes[i];
        }
        if(testpb)ofs_running << " nbnd in processor " << i << " is " << mpi_times << endl;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;
                double *send = new double[n];
                int count = 0;

                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

				// send the data to processor 0.
                MPI_Send(send,n,MPI_DOUBLE,0,tag,comm);

                delete[] send;
            }
			// third part to delete Z;
			delete[] Z;
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                double *ctmp = new double[NLOCAL];
                ZEROS(ctmp, NLOCAL);
                int tag = j;
                
				// Processor 0 receive the data from other processors.
				MPI_Recv(ctmp,n,MPI_DOUBLE,i,tag,comm,&status);

                for (int m=0; m<NLOCAL; m++)
                {
                    ctot[k][m]=ctmp[m];
//					ofs_running << " receive Z=" << ctmp[m] << endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        //MPI_Barrier(comm);
    }
    if(testpb)OUT(ofs_running,"Final k",k);
/*
if(myid==0){
	double *vect=new double[NLOCAL*NBANDS];
	double *form=new double[NBANDS*NBANDS];
	int x,y;
	for(x=0;x<NBANDS;x++){
		for(y=0;y<NLOCAL;y++){
			vect[x*NLOCAL+y]=ctot[x][y];
		}
	}
		char chT='T';
		char chN='N';
		int	ne = NBANDS;
		int m1, n1, k1;
		double ONE=1.0,ZERO=0.0;
		m1 = ne;
		n1 = ne;
		k1 = NLOCAL;
		dgemm_(&chT, &chN, &m1, &n1, &k1, &ONE, vect, &k1, vect, &k1, &ZERO, form, &m1);
		double di=0.0,oth=0.0;
		for(x=0;x<ne;x++){
			for(y=0;y<ne;y++){
				if(x==y){
					di+=fabs(form[x*ne+y]);
				}else{
					oth+=fabs(form[x*ne+y]);
				}
			}
		}
		di-=ne;
		printf("\n\n\ndi=%.16lf\n\n\nother=%.16lf\n\n\n",di,oth);
		//assert(0>1);
}
*/
MPI_Barrier(comm);

	// mohan add 2010-09-10
	// output the wave function if required.
	// this is a bad position to output wave functions.
	// but it works!
    if(this->out_lowf)
	{
		// read is in ../src_algorithms/wf_local.cpp
		stringstream ss;
		ss << global_out_dir << "LOWF_GAMMA_S" << CURRENT_SPIN+1 << ".dat";
		// mohan add 2012-04-03, because we need the occupations for the
		// first iteration. 
		Occupy::calculate_weights();
		WF_Local::write_lowf( ss.str(), ctot );//mohan add 2010-09-09        
	}

	// mohan add 2010-09-10
	// distribution of local wave functions 
	// to each processor.
	// only used for GAMMA_ONLY_LOCAL
	WF_Local::distri_lowf( ctot, wfc);

	// mohan 2010-09-26
	// distribution of augmented wave functions.
	// mohan fix bug 2011-01-03
	// the third dimension of WFC_GAMMA_aug may be zero.
	// mohan fix bug 2011-03-03
	// when the third dimension is zero, we don't need to call
	// this function, however, as proc 0, it needs all the procssors
	// to give reports to it.
	//	cout << " block distri_lowf_aug" << endl;
	// mohan update 2021-02-12, delte BFIELD option
	WF_Local::distri_lowf_aug( ctot, LOWF.WFC_GAMMA_aug[CURRENT_SPIN]); 

	// clean staff.
	if(myid==0)
	{
    	for (int i=0; i<NBANDS; i++)
    	{
        	delete[] ctot[i];
    	}
    	delete[] ctot;
	}

    time_t time_end = time(NULL);
    //ofs_running << " End   gath_eig Time : " << ctime(&time_end);
	OUT_TIME("gather eigenvalues",time_start,time_end);
    return;
}

void Pdiag_Basic::gath_eig_complex(MPI_Comm comm,int n,complex<double> **cc,complex<double> *Z, const int &ik)
{
    TITLE("Pdiag_Basic","gath_eig_complex");
    time_t time_start = time(NULL);
    //ofs_running << " Start gath_eig_complex Time : " << ctime(&time_start);

    int i,j,k;
    int nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);

    complex<double> **ctot;

	// mohan add 2010-07-03
	// the occupied bands are useless
	// for calculating charge density.
	if(DRANK>lastband_in_proc)
	{
		delete[] Z;
	}

	// first we need to collect all
	// the occupied bands.
	// NBANDS * NLOCAL	
	if(DRANK==0)
	{
		ctot = new complex<double>*[NBANDS];
    	for (int i=0; i<NBANDS; i++)
    	{
        	ctot[i] = new complex<double>[NLOCAL];
        	ZEROS(ctot[i], NLOCAL);
    	}
    	Memory::record("Pdiag_Basic","ctot",NBANDS*NLOCAL,"cdouble");
	}

	k=0;
    if (myid==0)
    {
        // mohan add nbnd0 2010-07-02
        int nbnd0 = -1;
        if (NBANDS < loc_sizes[0])
        {
			// means all bands in this processor
			// is needed ( is occupied)
            nbnd0 = NBANDS;
        }
        else
        {
			// means this processor only save
			// part of NBANDS.
            nbnd0 = loc_sizes[0];
        }
        if(testpb)ofs_running << " nbnd in processor 0 is " << nbnd0 << endl;

        for (i=0; i<nbnd0; i++)
        {
            for (j=0; j<NLOCAL; j++)
            {
				// change the order in processor 0.
				// the contribution from processor 0.
                ctot[k][j]=Z[j*loc_sizes[0]+i];
            }
            k++;
        }
		// Z is useless in processor 0 now.
		delete[] Z;
    }
    MPI_Barrier(comm);

	for (i=1; i<= this->lastband_in_proc; i++)
    {
        // mohan fix bug 2010-07-02
        // rows indicates the data structure of Z.
        // mpi_times indicates the data distribution
        // time, each time send a band.
        int rows = loc_sizes[i];
        int mpi_times;
        if (i==lastband_in_proc)
        {
            mpi_times = lastband_number;
        }
        else
        {
            mpi_times = loc_sizes[i];
        }
        if(testpb)ofs_running << " nbnd in processor " << i << " is " << mpi_times << endl;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;
                complex<double> *send = new complex<double>[n];
                int count = 0;

                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

				// send the data to processor 0.
                MPI_Send(send,n,mpicomplex,0,tag,comm);

                delete[] send;
            }
			// third part to delete Z;
			delete[] Z;
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                complex<double> *ctmp = new complex<double>[NLOCAL];
                ZEROS(ctmp, NLOCAL);
                int tag = j;
                
				// Processor 0 receive the data from other processors.
				MPI_Recv(ctmp,n,mpicomplex,i,tag,comm,&status);

                for (int m=0; m<NLOCAL; m++)
                {
                    ctot[k][m]=ctmp[m];
//					ofs_running << " receive Z=" << ctmp[m] << endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        //MPI_Barrier(comm);
    }
    if(testpb)OUT(ofs_running,"Final k",k);

	// output the wave function if required.
	// this is a bad position to output wave functions.
	// but it works!
	stringstream ss;
	ss << global_out_dir << "LOWF_K_" << ik+1 << ".dat";
    if(this->out_lowf)
	{
//		cout << " write the wave functions" << endl;
		WF_Local::write_lowf_complex( ss.str(), ctot, ik );//mohan add 2010-09-09        
	}

	// mohan add 2010-09-10
	// distribution of local wave functions 
	// to each processor.
	WF_Local::distri_lowf_complex( ctot, cc);

	// mohan 2010-09-26
	// distribution of augmented wave functions.
	// mohan fix bug 2011-01-03
	// the third dimension of WFC_GAMMA_aug may be zero.
	// mohan fix bug 2011-03-03
	// when the third dimension is zero, we don't need to call
	// this function, however, as proc 0, it needs all the procssors
	// to give reports to it.
	// mohan add 2012-01-09
	// 
	// for complex
	// mohan update 2021-02-12, delete BFIELD option
	WF_Local::distri_lowf_aug_complex( ctot, LOWF.WFC_K_aug[ik]); //mohan add 2012-01-09 
	
	
	// clean staff.
	if(DRANK==0)
	{
    	for (int i=0; i<NBANDS; i++)
    	{
        	delete[] ctot[i];
    	}
    	delete[] ctot;
	}

    time_t time_end = time(NULL);
//    ofs_running << " End   gath_eig_complex Time : " << ctime(&time_end);

	OUT_TIME("gather eigenvalues",time_start,time_end);
    return;
}
#endif

#ifdef __MPI
void Pdiag_Basic::gath_full_eig(MPI_Comm comm,int n,double **c,double *Z)
{
    TITLE("Pdiag_Basic","gath_full_eig");

    time_t time_start = time(NULL);
    //ofs_running << " Start gath_full_eig Time : " << ctime(&time_start);

    int i,j,k,incx=1;
    int *loc_sizes,loc_size,nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);
    loc_sizes =(int*)malloc(sizeof(int)*nprocs);
    loc_size=n/nprocs;
    for (i=0; i<nprocs; i++)
    {
        if (i<n%nprocs)
        {
            loc_sizes[i]=loc_size+1;
        }
        else
        {
            loc_sizes[i]=loc_size;
        }
    }
    if (myid==0)
    {
        k=0;
        for (i=0; i<loc_sizes[0]; i++)
        {
            for (j=0; j<n; j++)
            {
                // i : column index;
                // j : row index.
//              c[k*n+j]=Z[i*n+j];
                c[k][j]=Z[j*loc_sizes[0]+i];
                //			ofs_running << " Z=" << Z[i*n+j] << endl;
            }
            k++;
        }
    }
    MPI_Barrier(comm);

    for (i=1; i<nprocs; i++)
    {
        const int rows = loc_sizes[i];
        const int mpi_times = rows;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;

                double *send = new double[n];

                int count = 0;
                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

                MPI_Send(send,n,MPI_DOUBLE,0,tag,comm);

                delete[] send;
            }
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                double *ctmp = new double[n];
                ZEROS(ctmp, n);
                int tag = j;
                MPI_Recv(ctmp,n,MPI_DOUBLE,i,tag,comm,&status);

                for (int m=0; m<n; m++)
                {
                    c[k][m]=ctmp[m];
//					ofs_running << " receive Z=" << ctmp[m] << endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        MPI_Barrier(comm);
    }

    for (int i=0; i<NLOCAL; i++)
    {
        Parallel_Common::bcast_double(c[i],NLOCAL);
    }

    time_t time_end = time(NULL);
//    ofs_running << " End   gath_full_eig Time : " << ctime(&time_end);

    return;
}

void Pdiag_Basic::gath_full_eig_complex(MPI_Comm comm,int n,complex<double> **c,complex<double> *Z)
{
    TITLE("Pdiag_Basic","gath_full_eig_complex");

    time_t time_start = time(NULL);
    //ofs_running << " Start gath_full_eig_complex Time : " << ctime(&time_start);

    int i,j,k,incx=1;
    int *loc_sizes,loc_size,nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);
    loc_sizes =(int*)malloc(sizeof(int)*nprocs);
    loc_size=n/nprocs;
    for (i=0; i<nprocs; i++)
    {
        if (i<n%nprocs)
        {
            loc_sizes[i]=loc_size+1;
        }
        else
        {
            loc_sizes[i]=loc_size;
        }
    }
    if (myid==0)
    {
        k=0;
        for (i=0; i<loc_sizes[0]; i++)
        {
            for (j=0; j<n; j++)
            {
                // i : column index;
                // j : row index.
//              c[k*n+j]=Z[i*n+j];
                c[k][j]=Z[j*loc_sizes[0]+i];
                //			ofs_running << " Z=" << Z[i*n+j] << endl;
            }
            k++;
        }
    }
    MPI_Barrier(comm);

	for (i=1; i<nprocs; i++)
	{
		const int rows = loc_sizes[i];
		const int mpi_times = rows;
		if (myid==i)
		{
			for (j=0; j<mpi_times; j++)
			{
				int tag = j;

				complex<double> *send = new complex<double>[n];

				int count = 0;
				for (int m=0; m<rows*n; m++)
				{
					if (m%rows==j)
					{
						send[count] = Z[m];
						++count;
					}
				}

				MPI_Send(send,n,mpicomplex,0,tag,comm);

				delete[] send;
			}
		}
		else if (myid==0)
		{
			int col=0;
			for (j=0; j<mpi_times; j++)
			{
				complex<double> *ctmp = new complex<double>[n];
				ZEROS(ctmp, n);
				int tag = j;
				MPI_Recv(ctmp,n,mpicomplex,i,tag,comm,&status);

				for (int m=0; m<n; m++)
				{
					c[k][m]=ctmp[m];
					//					ofs_running << " receive Z=" << ctmp[m] << endl;
				}
				k++;

				delete[] ctmp;
			}
		}
		MPI_Barrier(comm);
	}

	for (int i=0; i<NLOCAL; i++)
	{
		Parallel_Common::bcast_complex_double(c[i],NLOCAL);
	}

    time_t time_end = time(NULL);
    //ofs_running << " End   gath_full_eig Time : " << ctime(&time_end);
	OUT_TIME("gather full eigenvalues",time_start,time_end);

    return;
}
#endif
