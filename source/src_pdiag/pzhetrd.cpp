#include "pzhetrd.h"
#include "../module_base/lapack_connector.h"

void pzhetrd(MPI_Comm comm2D, LocalMatrix loc_A,int N,int NB,
             std::complex<double> *A,double *diag,double *off_diag,std::complex<double> *norm,char uplo)
{
/*
 *   PSEPS routine (version 1.0) --
 *   Computer Network Information Center, CAS. 
 *   December 15, 2005
 *   Purpose
 *    ========
 *    pzhetrd reduces a Hermitian matrix A to Hermitian  tridiagonal form T by 
 *    an unitary similarity transformation inparallel: Q'*A*Q = T.
 *    Arguments
 *    =========
 *    Notes: Each global data object is described by an associated struct variable 
 *           loc_A. This struct stores the information required to establish the mapping 
 *           between an object element and its corresponding process 
 *    uplo   (global input) char
 *            = 'U':  Upper triangles of sub( A ) and sub( B ) are stored;
 *            = 'L':  Lower triangles of sub( A ) and sub( B ) are stored.
 *    comm_2D (input) MPI_Comm
 *            MPI 2D communicator 
 *    N       (global input) INTEGER
 *            The number of columns and rows to be operated on matrices A��N >= 0. 
 *    NB      (input) INTEGER
 *            blocked size of 2D blocked cyclic matrix
 *    A       (local input/local output) double precision std::complex pointer,
 *            pointer into the local memory to an array of dimension (loc_A.row_num, loc_A.col_num).
 *            On entry, this array contains the local pieces of the  N-by-N Hermitian distributed 
 *            matrix sub( A ). If UPLO = 'U',the leading N-by-N upper triangular part of sub( A ) 
 *            contains the upper triangular part of the matrix.  If UPLO = 'L', the leading N-by-N
 *            lower triangular part of sub( A ) contains the lower triangular part of the matrix.
 *    loc_A   (local input) struct Loc_A
 *            This struct avaible stores the information required to establish the mapping between
 *            an object element and its corresponding process and memory location.
 *    diag    (local output) double precisoin array,
 *            The diagonal elements of the tridiagonal matrix T:
 *    off_diag (local output) DOUBLE PRECISION array, 
 *            The off-diagonal elements of the tridiagonal matrix T
 *    norm    (local output) double precision array,
 *             This array contains the scalar factors TAU of  the elementary reflectors. 
 */        
	ModuleBase::TITLE("Parallel","pzhetrd");
	std::complex<double>  W_R[loc_A.row_num][NB],U_R[loc_A.row_num][NB],
		W_C[loc_A.col_num][NB],U_C[loc_A.col_num][NB],
		W_C1[loc_A.col_num][NB],U_C1[loc_A.col_num][NB];
	std::complex<double>  w_r[loc_A.row_num],w_c[loc_A.col_num],u_r[loc_A.row_num],u_c[loc_A.col_num],
		x_r[loc_A.row_num],x_c[loc_A.col_num];
	std::complex<double>   beta,alpha;
	std::complex<double>   u[N],tw[N],tu[N],x[N],p[N],w[N],q[N],s,g;
	double t;
	int     bt;
	double z[N],dg[N];
	MPI_Comm comm_col,comm_row;
	int dim[2],period[2],coord[2];
	int iacol,iarow,loc_i,loc_j;
	int i,j,k,in,first_col,end_col;
	int size,root_id;
	int tmp,tmp1,tmp2,pos;
	int lda,ldb,ldc;
	int m,lm,ln,outi,outj;
	int bn,fl=0,l,indx=1,indy=1,incx=1,incy;
	char transa,transb;
	std::complex<double>* L = new std::complex<double>[loc_A.col_num*loc_A.row_num];
	std::complex<double>* L1 = new std::complex<double>[loc_A.col_num*loc_A.row_num];


	MPI_Cart_get(comm2D,2,dim,period,coord);
	mpi_sub_col(comm2D,&comm_col);
	mpi_sub_row(comm2D,&comm_row);

	int ha;
	bn=N/NB;
	if (bn*NB<N) bn++;

	for (i=loc_A.row_pos; i<loc_A.row_num; i++)
		for (j=loc_A.col_pos; j<loc_A.col_num; j++)
		{
			L1[i*loc_A.col_num+j]=L[i*loc_A.col_num+j]=std::complex<double> (0.0,0.0);
		}

	if (uplo=='U'||uplo=='u')
		for (bt=0; bt<bn; bt++)
		{
			//      printf("\nuplo=%c,bt=%d,bn=%d\n",uplo,bt,bn);
			for (i=0; i<loc_A.row_num; i++)
				for (j=0; j<NB; j++)
				{
					U_R[i][j]=W_R[i][j]=std::complex<double> (0.0,0.0);
				}

			for (i=0; i<loc_A.col_num; i++)
				for (j=0; j<NB; j++)
				{
					U_C[i][j]=W_C[i][j]=std::complex<double> (0.0,0.0);
					U_C1[i][j]=W_C1[i][j]=std::complex<double> (0.0,0.0);
				}
			//determain algorithm block size m
			first_col=bt*NB;
			end_col=(first_col+NB-1>N-2)?(N-2):(first_col+NB-1);
			m=end_col-first_col+1;
			indxg2p(comm2D,NB,first_col,first_col,&iarow,&iacol);          //diagonal element
			if (coord[0]==iarow) outi=loc_A.row_pos+m;
			else outi=loc_A.row_pos;
			if (coord[1]==iacol) outj=loc_A.col_pos+m;
			else outj=loc_A.col_pos;
			for (i=loc_A.row_pos; i<loc_A.row_num; i++)
			{
				tmp1=loc_A.row_set[i];
				for (j=loc_A.col_pos; j<loc_A.col_num; j++)
				{
					tmp2=loc_A.col_set[j];
					if (tmp1<=tmp2)
					{
						L[i*loc_A.col_num+j]=A[i*loc_A.col_num+j];
					}
					if (tmp1<tmp2)
					{
						L1[i*loc_A.col_num+j]=conj(A[i*loc_A.col_num+j]);
					}           
				}
			}

			//parallel translation in the panel block
			//printf("first_col=%d,end_col=%d",first_col,end_col);
			for (k=first_col; k<=end_col; k++)
			{
				for (i=k; i<N; i++)
				{ 
					x[i]=p[i]=u[i]=std::complex<double> (0.0,0.0);
				} 
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					x_r[i]=u_r[i]=w_r[i]=std::complex<double> (0.0,0.0);
				}
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					x_c[i]=u_c[i]=w_c[i]=std::complex<double> (0.0,0.0);
				}
				indxg2p(comm2D,NB,k,k,&iarow,&iacol);
				if (iacol==coord[1]) loc_A.col_pos++;
				//      generate u,tau
				if (iarow==coord[0])
				{
					for (i=loc_A.col_pos; i<loc_A.col_num; i++)
					{
						u_c[i]=conj(A[loc_A.row_pos*loc_A.col_num+i]);
					}
					MPI_Bcast(&u_c[loc_A.col_pos],loc_A.col_num-loc_A.col_pos,MPI_DOUBLE_COMPLEX,iarow,comm_col);
				}
				else
					MPI_Bcast(&u_c[loc_A.col_pos],loc_A.col_num-loc_A.col_pos,MPI_DOUBLE_COMPLEX,iarow,comm_col);
				indxg2p(comm2D,NB,k,k+1,&iarow,&iacol);
				indxg2l(k,k+1,NB,dim[0],dim[1],&loc_i,&loc_j);
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					tmp=loc_A.col_set[i];
					p[tmp]=u_c[i];
				}
				MPI_Allreduce(&p[k+1],&u[k+1],N-k-1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_row);
				alpha=u[k+1];
				lm=N-k-1;
				zlarfg_(&lm,&alpha,&u[k+2],&incx,&g);
				norm[k]=g;
				if (coord[1]==iacol)
				{
					u_c[loc_j]=std::complex<double> (1.0,0.0);
				}
				z[k]=alpha.real();
				u[k+1]=std::complex<double> (1.0,0.0);
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					tmp=loc_A.col_set[i];
					u_c[i]=u[tmp];
				}
				if (iarow==coord[0])
					for (i=loc_A.col_pos; i<loc_A.col_num; i++)
						A[loc_A.row_pos*loc_A.col_num+i]=u_c[i];
				if (iarow==coord[0]) loc_A.row_pos++;
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					u_r[i]=u[tmp];
				}
				if (k==N-2)  {fl=1;break;}

				for (i=k+1; i<N; i++)
				{
					x[i]=p[i]=std::complex<double> (0.0,0.0);
				}
				for (i=loc_A.row_pos; i<outi; i++)
				{
					tmp1=loc_A.row_set[i];
					for (j=loc_A.col_pos; j<loc_A.col_num; j++)
					{
						tmp2=loc_A.col_set[j];
						if (tmp1<=tmp2)
							L[i*loc_A.col_num+j]=A[i*loc_A.col_num+j];
						if (tmp1<tmp2)
						{
							L1[i*loc_A.col_num+j]=conj(A[i*loc_A.col_num+j]);
						}
					}
				}
				//local matrix-std::vector production Au
				transa='t';
				lm=loc_A.row_num-loc_A.row_pos;
				ln=loc_A.col_num-loc_A.col_pos;
				beta=std::complex<double> (1.0,0.0);
				lda=loc_A.col_num;
				pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
				incx=1;
				incy=1;
				zgemv_(&transa,&ln,&lm,&g,&L[pos],&lda,&u_c[loc_A.col_pos],&incx,&beta,
						&x_r[loc_A.row_pos],&incy);
				transa='n';
				zgemv_(&transa,&ln,&lm,&g,&L1[pos],&lda,&u_r[loc_A.row_pos],&incx,&beta,
						&x_c[loc_A.col_pos],&incy);
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					x[tmp]=x[tmp]+x_r[i];
				}
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					tmp=loc_A.col_set[i];
					x[tmp]=x[tmp]+x_c[i];
				}
				//correct for out of date entries of A
				for (i=0; i<k-first_col; i++)
				{
					tu[i]=tw[i]=std::complex<double> (0.0,0.0);
				}
				transa='n';
				lm=(loc_A.col_num-outj)>0?(loc_A.col_num-outj):0;
				ln=k-first_col;
				beta=std::complex<double> (0.0,0.0);
				lda=NB;
				incx=1;
				incy=1;
				zgemv_(&transa,&ln,&lm,&g,&W_C1[outj][0],&lda,&u_c[outj],&incx,&beta,tw,&incy);
				zgemv_(&transa,&ln,&lm,&g,&U_C1[outj][0],&lda,&u_c[outj],&incx,&beta,tu,&incy);
				for (i=outi; i<loc_A.row_num; i++)
				{
					x_r[i]=std::complex<double> (0.0,0.0);
				}
				transa='t';
				lm=(loc_A.row_num-outi)>0?(loc_A.row_num-outi):0;
				ln=k-first_col;
				alpha=std::complex<double> (1.0,0.0);
				beta=std::complex<double> (1.0,0.0);
				lda=NB;
				incx=1;
				incy=1;
				zgemv_(&transa,&ln,&lm,&alpha,&W_R[outi][0],&lda,tu,&incx,&beta,&x_r[outi],&incy);
				zgemv_(&transa,&ln,&lm,&alpha,&U_R[outi][0],&lda,tw,&incx,&beta,&x_r[outi],&incy);
				for (i=outi; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					x[tmp]=x[tmp]-x_r[i];
				}
				MPI_Allreduce(&x[k+1],&q[k+1],N-k-1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_col);
				MPI_Allreduce(&q[k+1],&p[k+1],N-k-1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_row);
				s=std::complex<double> (0.0,0.0);
				lm=N-k-1;
				//printf("\nincx=%d,incy=%d,lm=%d \n",incx,incy,lm);
				//for(int mm=0;mm<N;mm++)
				//printf("u[%d]=%lf+i%lf,p[%d]=%lf+i%lf ",mm,u[mm].real(),u[mm].imag(),mm,p[mm].real(),p[mm].imag());

				//               s=zdotc_(&lm,&p[k+1],&incx,&u[k+1],&incy);
				for(int yang=0;yang<lm;yang++)
					s+=conj(p[k+1+yang])*u[k+1+yang];           
				//printf("\ns=%lf+i%lf\n",s.real(),s.imag());
				s=std::complex <double> (-s.real()/2.0,-s.imag()/2.0);
				alpha=std::complex<double> (s.real()*g.real()-s.imag()*g.imag(),s.real()*g.imag()+s.imag()*g.real());                 
				//              printf("\nalpha=%lf+i%lf\n",alpha.real(),alpha.imag());
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					w_r[i]=std::complex<double>((p[tmp].real()+(alpha.real()*u_r[i].real()-alpha.imag()*u_r[i].imag())),(p[tmp].imag()+(alpha.real()*u_r[i].imag()+alpha.imag()*u_r[i].real())));
					//                 printf("w_r[%d]=%lf+i%lf",i,w_r[i].real(),w_r[i].imag());         
				}
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					tmp=loc_A.col_set[i];
					w_c[i]=std::complex<double> ((p[tmp].real()+(alpha.real()*u_c[i].real()-alpha.imag()*u_c[i].imag())),(p[tmp].imag()+(alpha.real()*u_c[i].imag()+alpha.imag()*u_c[i].real())));
					//                              printf("w_c[%d]=%lf+i%lf ",i,w_c[i].real(),w_c[i].imag());   
				}
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					U_R[i][k-first_col]=u_r[i];
					W_R[i][k-first_col]=w_r[i];
				}
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					U_C[i][k-first_col]=u_c[i];
					W_C[i][k-first_col]=w_c[i];
					U_C1[i][k-first_col]=conj(u_c[i]);
					W_C1[i][k-first_col]=conj(w_c[i]);
				}
				for(i=loc_A.row_pos;i<loc_A.row_num;i++)
				{
					w_r[i]=conj(w_r[i]);
					u_r[i]=conj(u_r[i]);
				}
				for(i=loc_A.col_pos;i<loc_A.col_num;i++)
				{
					w_c[i]=conj(w_c[i]);
					u_c[i]=conj(u_c[i]);
				}
				//update remainder of panel block of A
				lm=(outi-loc_A.row_pos)>0?(outi-loc_A.row_pos):0;
				ln=loc_A.col_num-loc_A.col_pos;
				alpha=std::complex<double>  (-1.0,0.0);
				lda=loc_A.col_num;
				pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
				incy=1;
				zgerc_(&ln,&lm,&alpha,&u_c[loc_A.col_pos],&incx,&w_r[loc_A.row_pos],&incy,&A[pos],&lda);
				zgerc_(&ln,&lm,&alpha,&w_c[loc_A.col_pos],&incx,&u_r[loc_A.row_pos],&incy,&A[pos],&lda);
				MPI_Barrier(DIAG_HPSEPS_WORLD);

			}

			//rank-2b update of the submatrix A22
			transa='C';
			transb='n';
			lm=loc_A.row_num-loc_A.row_pos;
			ln=loc_A.col_num-loc_A.col_pos;
			in=end_col-first_col+1;
			alpha=std::complex<double> (-1.0,0.0);
			beta=std::complex<double> (1.0,0.0);
			ldb=NB;
			ldc=loc_A.col_num;
			pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
			zgemm_(&transa,&transb,&ln,&lm,&in,&alpha,&W_C[loc_A.col_pos][0],&ldb,
					&U_R[loc_A.row_pos][0],&ldb,&beta,&A[pos],&ldc);
			zgemm_(&transa,&transb,&ln,&lm,&in,&alpha,&U_C[loc_A.col_pos][0],&ldb,
					&W_R[loc_A.row_pos][0],&ldb,&beta,&A[pos],&ldc);
			MPI_Barrier(DIAG_HPSEPS_WORLD);


		}
	//  printf("hello\n");
	if (uplo=='L'||uplo=='l')
		for (bt=0; bt<bn; bt++)
		{

			for (i=0; i<loc_A.row_num; i++)
				for (j=0; j<NB; j++)
				{
					U_R[i][j]=W_R[i][j]=std::complex<double> (0.0,0.0);
				}
			for (i=0; i<loc_A.col_num; i++)
				for (j=0; j<NB; j++)
				{
					U_C[i][j]=W_C[i][j]=std::complex<double> (0.0,0.0);
				}
			// determain algorithm block size m
			first_col=bt*NB;
			end_col=(first_col+NB-1>N-2)?(N-2):(first_col+NB-1);
			m=end_col-first_col+1;
			indxg2p(comm2D,NB,first_col,first_col,&iarow,&iacol);//diagonal element
			if (coord[0]==iarow) outi=loc_A.row_pos+m;
			else outi=loc_A.row_pos;
			if (coord[1]==iacol) outj=loc_A.col_pos+m;
			else outj=loc_A.col_pos;
			for (i=loc_A.row_pos; i<loc_A.row_num; i++)
			{
				tmp1=loc_A.row_set[i];
				for (j=loc_A.col_pos; j<loc_A.col_num; j++)
				{
					tmp2=loc_A.col_set[j];
					if (tmp1>=tmp2)
						L[i*loc_A.col_num+j]=A[i*loc_A.col_num+j];
					if (tmp1>tmp2)
					{
						L1[i*loc_A.col_num+j]=conj(A[i*loc_A.col_num+j]);
					}
				}
			}
			//parallel translation in the panel block
			for (k=first_col; k<=end_col; k++)
			{
				for (i=k; i<N; i++)
				{
					x[i]=p[i]=u[i]=std::complex<double> (0.0,0.0);
				}
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
					x_r[i]=u_r[i]=w_r[i]=std::complex<double> (0.0,0.0);
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
					x_c[i]=u_c[i]=w_c[i]=std::complex<double> (0.0,0.0);
				indxg2p(comm2D,NB,k,k,&iarow,&iacol);
				if (iarow==coord[0]) loc_A.row_pos++;
				//generate u, tau
				if (coord[1]==iacol)
				{
					for (i=loc_A.row_pos; i<loc_A.row_num; i++)
						u_r[i]=A[i*loc_A.col_num+loc_A.col_pos];
					MPI_Bcast(&u_r[loc_A.row_pos],loc_A.row_num-loc_A.row_pos,MPI_DOUBLE_COMPLEX,iacol,comm_row);
				}
				else
					MPI_Bcast(&u_r[loc_A.row_pos],loc_A.row_num-loc_A.row_pos,MPI_DOUBLE_COMPLEX,iacol,comm_row);
				indxg2p(comm2D,NB,k+1,k,&iarow,&iacol);
				indxg2l(k+1,k,NB,dim[0],dim[1],&loc_i,&loc_j);
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					p[tmp]=u_r[i];
				}
				MPI_Allreduce(&p[k+1],&u[k+1],N-k-1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_col);
				alpha=u[k+1];
				lm=N-k-1;
				zlarfg_(&lm,&alpha,&u[k+2],&incx,&g);
				norm[k]=g;
				if (coord[0]==iarow)
				{
					u_r[loc_i]=std::complex<double> (1.0,0.0);
				}
				z[k]=alpha.real();
				u[k+1]=std::complex<double> (1.0,0.0);
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					u_r[i]=u[tmp];
				}
				if (coord[1]==iacol)
				{
					for (i=loc_A.row_pos; i<loc_A.row_num; i++)
						A[i*loc_A.col_num+loc_A.col_pos]=u_r[i];
				}
				if (iacol==coord[1]) loc_A.col_pos++;
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					tmp=loc_A.col_set[i];
					u_c[i]=u[tmp];
				}
				if (k==N-2) {
					fl=1;
					break;
				}
				for (i=k+1; i<N; i++)
					x[i]=p[i]=std::complex<double> (0.0,0.0);

				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					tmp1=loc_A.row_set[i];
					for (j=loc_A.col_pos; j<outj; j++)
					{
						tmp2=loc_A.col_set[j];
						if (tmp1>=tmp2)
							L[i*loc_A.col_num+j]=A[i*loc_A.col_num+j];
						if (tmp1>tmp2)
							L1[i*loc_A.col_num+j]=conj(A[i*loc_A.col_num+j]);
					}
				}
				//local matrix-std::vector production
				transa='t';
				lm=loc_A.row_num-loc_A.row_pos;
				ln=loc_A.col_num-loc_A.col_pos;
				beta=std::complex<double> (0.0,0.0);
				lda=loc_A.col_num;
				pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
				incx=1;
				incy=1;
				zgemv_(&transa,&ln,&lm,&g,&L[pos],&lda,&u_c[loc_A.col_pos],&incx,&beta,
						&x_r[loc_A.row_pos],&incy);
				transa='n';
				zgemv_(&transa,&ln,&lm,&g,&L1[pos],&lda,&u_r[loc_A.row_pos],&incx,&beta,
						&x_c[loc_A.col_pos],&incy);
				//correct for out of date entries of A
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					x[tmp]=x[tmp]+x_r[i];
				}
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					tmp=loc_A.col_set[i];
					x[tmp]=x[tmp]+x_c[i];
				}
				for (i=0; i<k-first_col; i++)
					tu[i]=tw[i]=std::complex<double> (0.0,0.0);
				transa='n';
				lm=(loc_A.col_num-outj)>0?(loc_A.col_num-outj):0;
				ln=k-first_col;
				beta=std::complex<double> (0.0,0.0);
				lda=NB;
				incx=1;
				incy=1;
				zgemv_(&transa,&ln,&lm,&g,&W_C1[outj][0],&lda,&u_c[outj],&incx,&beta,tw,&incy);
				zgemv_(&transa,&ln,&lm,&g,&U_C1[outj][0],&lda,&u_c[outj],&incx,&beta,tu,&incy);
				for (i=outi; i<loc_A.row_num; i++)
					x_r[i]=std::complex<double> (0.0,0.0);
				transa='t';
				lm=(loc_A.row_num-outi)>0?(loc_A.row_num-outi):0;
				ln=k-first_col;
				alpha=std::complex<double> (1.0,0.0);
				beta=std::complex<double> (1.0,0.0);
				lda=NB;
				incx=1;
				incy=1;
				zgemv_(&transa,&ln,&lm,&alpha,&W_R[outi][0],&lda,tu,&incx,&beta,&x_r[outi],&incy);
				zgemv_(&transa,&ln,&lm,&alpha,&U_R[outi][0],&lda,tw,&incx,&beta,&x_r[outi],&incy);
				for (i=outi; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					x[tmp]=x[tmp]-x_r[i];
				}
				MPI_Allreduce(&x[k+1],&q[k+1],N-k-1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_col);
				MPI_Allreduce(&q[k+1],&p[k+1],N-k-1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_row);
				s=std::complex<double> (0.0,0.0);
				lm=N-k-1;

				// mohan modify 2010-10-23
				zdotc_(&s, &lm,&p[k+1],&incx,&u[k+1],&incy);
				//   printf("\ns=%lf+i%lf\n",s.real(),s.imag());
				s=std::complex<double>(-s.real()/2.0,-s.imag()/2.0);
				alpha=std::complex<double>( s.real()*g.real()-s.imag()*g.imag(),s.real()*g.imag()+s.imag()*g.real());

				//update remainder of panel block of A
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					tmp=loc_A.row_set[i];
					w_r[i]=std::complex<double>((p[tmp].real()+(alpha.real()*u_r[i].real()-alpha.imag()*u_r[i].imag())),(p[tmp].imag()+(alpha.real()*u_r[i].imag()+alpha.imag()*u_r[i].real())));
				}
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					tmp=loc_A.col_set[i];
					w_c[i]=std::complex<double>((p[tmp].real()+(alpha.real()*u_c[i].real()-alpha.imag()*u_c[i].imag())),(p[tmp].imag()+(alpha.real()*u_c[i].imag()+alpha.imag()*u_c[i].real())));
				}
				lm=loc_A.row_num-loc_A.row_pos;
				ln=(outj-loc_A.col_pos)>0?(outj-loc_A.col_pos):0;
				alpha=std::complex<double>(-1.0,0.0);
				lda=loc_A.col_num;
				pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
				incy=1;
				for (i=loc_A.row_pos; i<loc_A.row_num; i++)
				{
					U_R[i][k-first_col]=u_r[i];
					W_R[i][k-first_col]=w_r[i];
				}
				for (i=loc_A.col_pos; i<loc_A.col_num; i++)
				{
					U_C[i][k-first_col]=u_c[i];
					W_C[i][k-first_col]=w_c[i];
					U_C1[i][k-first_col]=conj(u_c[i]);
					W_C1[i][k-first_col]=conj(w_c[i]);
				}
				for(i=loc_A.row_pos;i<loc_A.row_num;i++)
				{
					w_r[i]=conj(w_r[i]);
					u_r[i]=conj(u_r[i]);
				}
				for(i=loc_A.col_pos;i<loc_A.col_num;i++)
				{
					w_c[i]=conj(w_c[i]);
					u_c[i]=conj(u_c[i]);
				}
				zgerc_(&ln,&lm,&alpha,&u_c[loc_A.col_pos],&incx,&w_r[loc_A.row_pos],&incy,&A[pos],&lda);
				zgerc_(&ln,&lm,&alpha,&w_c[loc_A.col_pos],&incx,&u_r[loc_A.row_pos],&incy,&A[pos],&lda);
				MPI_Barrier(DIAG_HPSEPS_WORLD);
			}
			// rank-2b update of the submatrix A22
			transa='C';
			transb='n';
			lm=loc_A.row_num-loc_A.row_pos;
			ln=loc_A.col_num-loc_A.col_pos;
			in=end_col-first_col+1;
			alpha=  std::complex<double> (-1.0,0.0);
			beta= std::complex<double> (1.0,0.0);
			ldb=NB;
			ldc=loc_A.col_num;
			pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;

			zgemm_(&transa,&transb,&ln,&lm,&in,&alpha,&W_C[loc_A.col_pos][0],&ldb,
					&U_R[loc_A.row_pos][0],&ldb,&beta,&A[pos],&ldc);
			zgemm_(&transa,&transb,&ln,&lm,&in,&alpha,&U_C[loc_A.col_pos][0],&ldb,
					&W_R[loc_A.row_pos][0],&ldb,&beta,&A[pos],&ldc);
			MPI_Barrier(DIAG_HPSEPS_WORLD);
		}

	//    printf("hello");

	for (i=0; i<N; i++)
		dg[i]=0.0;
	for (i=0; i<loc_A.row_num; i++)
		for (j=0; j<loc_A.col_num; j++)
		{
			tmp1=loc_A.row_set[i];
			tmp2=loc_A.col_set[j];
			if (tmp1==tmp2)
				dg[tmp1]=A[i*loc_A.col_num+j].real();
		}


	MPI_Allreduce(dg,diag,N,MPI_DOUBLE,MPI_SUM,comm_row);
	MPI_Allreduce(diag,dg,N,MPI_DOUBLE,MPI_SUM,comm_col);
	//  for(int haha=0;haha<N;haha++)
	// {
	//     printf("diag[%d]=%lf,dg[%d]=%lf,off_diag[%d]=%lf,z[%d]=%lf\n ",haha,diag[haha],haha,dg[haha],haha,off_diag[haha],haha,z[haha]);
	// }  
	BlasConnector::copy(N,dg,indx,diag,indy);
	BlasConnector::copy(N,z,indx,off_diag,indy);
	//  for(int haha=0;haha<N;haha++)
	//  {
	//     printf("diag[%d]=%lf,dg[%d]=%lf,off_diag[%d]=%lf,z[%d]=%lf\n ",haha,diag[haha],haha,dg[haha],haha,off_diag[haha],haha,z[haha]);
	// }  
	off_diag[N-1]=0.0;
	MPI_Comm_free(&comm_row);
	MPI_Comm_free(&comm_col);

	delete[] L;
	delete[] L1;
}

