#include "pdsytrd.h"

void pdsytrd(MPI_Comm comm2D, LocalMatrix loc_A,int N,int NB,
             double *a,double *diag,double *off_diag,double  *norm,char uplo)
{
    /*
     * PSEPS routine (version 1.0) --
     * Computer Network Information Center, CAS.
     * December 15, 2005
     *
     * Purpose
     * ========
     * pzhetrd reduces a Hermitian matrix A to Hermitian  tridiagonal form T by
     * an unitary similarity transformation inparallel: Q'*A*Q = T.
     * Arguments
     * =========
     * Notes: Each global data object is described by an associated struct variable
     * loc_A. This struct stores the information required to establish the mapping
     * between an object element and its corresponding process
     * uplo    (global input) char
     *         = 'U':  Upper triangles of sub( A ) and sub( B ) are stored;
     *         = 'L':  Lower triangles of sub( A ) and sub( B ) are stored.
     * comm2D  (input) MPI_Comm
     * MPI 2D  communicator
     * N       (global input) INTEGER
     *         The number of columns and rows to be operated on matrices A��N >= 0.
     * NB      (input) INTEGER
     *         blocked size of 2D blocked cyclic matrix
     * A       (local input/local output) double precision pointer,
     *          pointer into the local memory to an array of dimension
     *         (loc_A.row_num, loc_A.col_num).
     *         On entry, this array contains the local pieces of the
     *         N-by-N Hermitian distributed matrix sub( A ). If UPLO = 'U',
     *         the leading N-by-N upper triangular part of sub( A ) contains
     *         the upper triangular part of the matrix.  If UPLO = 'L', the
     *         leading N-by-N lower triangular part of sub( A ) contains
     *         the lower triangular part of the matrix.
     * loc_A   (local input) struct Loc_A
     *         This struct avaible stores the information required to establish
     *         the mapping between an object element and its corresponding
     *         process and memory location.
     * diag    (local output) double precisoin array,
     *         The diagonal elements of the tridiagonal matrix T:
     * off_diag (local output) DOUBLE PRECISION array,
     *         The off-diagonal elements of the tridiagonal matrix T
     * norm    (local output) double precision array,
     *         This array contains the scalar factors TAU of  the elementary
     *         reflectors.
     */
    double  W_R[loc_A.row_num][NB],U_R[loc_A.row_num][NB],
    W_C[loc_A.col_num][NB],U_C[loc_A.col_num][NB];
    double  w_r[loc_A.row_num],w_c[loc_A.col_num],u_r[loc_A.row_num],u_c[loc_A.col_num],
    x_r[loc_A.row_num],x_c[loc_A.col_num];
    double  *L,*L1;
    double  beta,alpha;
    double  u[N],tw[N],tu[N],x[N],p[N],w[N],q[N],s,g;
    int     bt;
    double z[N],dg[N];
    double t1,t2,off;
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


    MPI_Cart_get(comm2D,2,dim,period,coord);
    mpi_sub_col(comm2D,&comm_col);
    mpi_sub_row(comm2D,&comm_row);

    L=(double*)malloc(loc_A.col_num*loc_A.row_num*sizeof(double));
    L1=(double*)malloc(loc_A.col_num*loc_A.row_num*sizeof(double));

    bn=N/NB;
    if (bn*NB<N) bn++;

    for (i=loc_A.row_pos; i<loc_A.row_num; i++)
        for (j=loc_A.col_pos; j<loc_A.col_num; j++)
            L1[i*loc_A.col_num+j]=L[i*loc_A.col_num+j]=0.0;

    if (uplo=='U'||uplo=='u')
        for (bt=0; bt<bn; bt++)
        {
            for (i=0; i<loc_A.row_num; i++)
                for (j=0; j<NB; j++)
                    U_R[i][j]=W_R[i][j]=0.0;

            for (i=0; i<loc_A.col_num; i++)
                for (j=0; j<NB; j++)
                    U_C[i][j]=W_C[i][j]=0.0;
            /*determain algorithm block size m*/
            first_col=bt*NB;
            end_col=(first_col+NB-1>N-2)?(N-2):(first_col+NB-1);
            m=end_col-first_col+1;
            indxg2p(comm2D,NB,first_col,first_col,&iarow,&iacol);/*diagonal element*/
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
                        L[i*loc_A.col_num+j]=a[i*loc_A.col_num+j];
                    if (tmp1<tmp2)
                        L1[i*loc_A.col_num+j]=a[i*loc_A.col_num+j];
                }
            }
            /*parallel translation in the panel block*/
            for (k=first_col; k<=end_col; k++)
            {
                for (i=k; i<N; i++)
                    x[i]=p[i]=u[i]=0.0;
                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                    x_r[i]=u_r[i]=w_r[i]=0.0;
                for (i=loc_A.col_pos; i<loc_A.col_num; i++)
                    x_c[i]=u_c[i]=w_c[i]=0.0;
                indxg2p(comm2D,NB,k,k,&iarow,&iacol);
                if (iacol==coord[1]) loc_A.col_pos++;
                /*generate u,tau*/
                if (iarow==coord[0])
                {
                    for (i=loc_A.col_pos; i<loc_A.col_num; i++)
                        u_c[i]=a[loc_A.row_pos*loc_A.col_num+i];
                    MPI_Bcast(&u_c[loc_A.col_pos],loc_A.col_num-loc_A.col_pos,MPI_DOUBLE,iarow,comm_col);
                }
                else
                    MPI_Bcast(&u_c[loc_A.col_pos],loc_A.col_num-loc_A.col_pos,MPI_DOUBLE,iarow,comm_col);
                indxg2p(comm2D,NB,k,k+1,&iarow,&iacol);
                indxg2l(k,k+1,NB,dim[0],dim[1],&loc_i,&loc_j);
                for (i=loc_A.col_pos; i<loc_A.col_num; i++)
                {
                    tmp=loc_A.col_set[i];
                    p[tmp]=u_c[i];
                }
                MPI_Allreduce(&p[k+1],&u[k+1],N-k-1,MPI_DOUBLE,MPI_SUM,comm_row);
                alpha=u[k+1];
                lm=N-k-1;
                dlarfg_(&lm,&alpha,&u[k+2],&incx,&g);
                norm[k]=g;
                if (coord[1]==iacol)
                    u_c[loc_j]=1.0;
                z[k]=alpha;
                u[k+1]=1.0;
                for (i=loc_A.col_pos; i<loc_A.col_num; i++)
                {
                    tmp=loc_A.col_set[i];
                    u_c[i]=u[tmp];
                }
                if (iarow==coord[0])
                    for (i=loc_A.col_pos; i<loc_A.col_num; i++)
                        a[loc_A.row_pos*loc_A.col_num+i]=u_c[i];
                if (iarow==coord[0]) loc_A.row_pos++;
                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                {
                    tmp=loc_A.row_set[i];
                    u_r[i]=u[tmp];
                }
                if (k==N-2) {
                    fl=1;
                    break;
                }
                for (i=k+1; i<N; i++)
                    x[i]=p[i]=0.0;
                for (i=loc_A.row_pos; i<outi; i++)
                {
                    tmp1=loc_A.row_set[i];
                    for (j=loc_A.col_pos; j<loc_A.col_num; j++)
                    {
                        tmp2=loc_A.col_set[j];
                        if (tmp1<=tmp2)
                            L[i*loc_A.col_num+j]=a[i*loc_A.col_num+j];
                        if (tmp1<tmp2)
                            L1[i*loc_A.col_num+j]=a[i*loc_A.col_num+j];
                    }
                }
                /*local matrix-std::vector production Au*/
                transa='t';
                lm=loc_A.row_num-loc_A.row_pos;
                ln=loc_A.col_num-loc_A.col_pos;
                beta=1.0;
                lda=loc_A.col_num;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
                incx=1;
                incy=1;
                dgemv_(&transa,&ln,&lm,&g,&L[pos],&lda,&u_c[loc_A.col_pos],&incx,&beta,
                       &x_r[loc_A.row_pos],&incy);
                transa='n';
                dgemv_(&transa,&ln,&lm,&g,&L1[pos],&lda,&u_r[loc_A.row_pos],&incx,&beta,
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
                /*correct for out of date entries of A*/
                for (i=0; i<k-first_col; i++)
                    tu[i]=tw[i]=0.0;
                transa='n';
                lm=(loc_A.col_num-outj)>0?(loc_A.col_num-outj):0;
                ln=k-first_col;
                beta=0.0;
                lda=NB;
                incx=1;
                incy=1;
                dgemv_(&transa,&ln,&lm,&g,&W_C[outj][0],&lda,&u_c[outj],&incx,&beta,tw,&incy);
                dgemv_(&transa,&ln,&lm,&g,&U_C[outj][0],&lda,&u_c[outj],&incx,&beta,tu,&incy);
                for (i=outi; i<loc_A.row_num; i++)
                    x_r[i]=0.0;
                transa='t';
                lm=(loc_A.row_num-outi)>0?(loc_A.row_num-outi):0;
                ln=k-first_col;
                alpha=1.0;
                beta=1.0;
                lda=NB;
                incx=1;
                incy=1;
                dgemv_(&transa,&ln,&lm,&alpha,&W_R[outi][0],&lda,tu,&incx,&beta,&x_r[outi],&incy);
                dgemv_(&transa,&ln,&lm,&alpha,&U_R[outi][0],&lda,tw,&incx,&beta,&x_r[outi],&incy);
                for (i=outi; i<loc_A.row_num; i++)
                {
                    tmp=loc_A.row_set[i];
                    x[tmp]=x[tmp]-x_r[i];
                }
                MPI_Allreduce(&x[k+1],&q[k+1],N-k-1,MPI_DOUBLE,MPI_SUM,comm_col);
                MPI_Allreduce(&q[k+1],&p[k+1],N-k-1,MPI_DOUBLE,MPI_SUM,comm_row);
                s=0.0;
                lm=N-k-1;
                for (i=k+1; i<N; i++)
                    s+=p[i]*u[i];
                s=-s/2.0;
                alpha=s*g;
                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                {
                    tmp=loc_A.row_set[i];
                    w_r[i]=p[tmp]+(alpha*u_r[i]);
                }
                for (i=loc_A.col_pos; i<loc_A.col_num; i++)
                {
                    tmp=loc_A.col_set[i];
                    w_c[i]=p[tmp]+(alpha*u_c[i]);
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
                }
                /*update remainder of panel block of A*/
                lm=(outi-loc_A.row_pos)>0?(outi-loc_A.row_pos):0;
                ln=loc_A.col_num-loc_A.col_pos;
                alpha=-1.0;
                lda=loc_A.col_num;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
                incy=1;
                dger_(&ln,&lm,&alpha,&u_c[loc_A.col_pos],&incx,&w_r[loc_A.row_pos],&incy,&a[pos],&lda);
                dger_(&ln,&lm,&alpha,&w_c[loc_A.col_pos],&incx,&u_r[loc_A.row_pos],&incy,&a[pos],&lda);
                MPI_Barrier(DIAG_HPSEPS_WORLD);
            }
            /*rank-2b update of the submatrix A22*/
            transa='t';
            transb='n';
            lm=loc_A.row_num-loc_A.row_pos;
            ln=loc_A.col_num-loc_A.col_pos;
            in=end_col-first_col+1;
            alpha=-1.0;
            beta=1.0;
            ldb=NB;
            ldc=loc_A.col_num;
            pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
            dgemm_(&transa,&transb,&ln,&lm,&in,&alpha,&W_C[loc_A.col_pos][0],&ldb,
                   &U_R[loc_A.row_pos][0],&ldb,&beta,&a[pos],&ldc);
            dgemm_(&transa,&transb,&ln,&lm,&in,&alpha,&U_C[loc_A.col_pos][0],&ldb,
                   &W_R[loc_A.row_pos][0],&ldb,&beta,&a[pos],&ldc);
            MPI_Barrier(DIAG_HPSEPS_WORLD);
        }
    if (uplo=='L'||uplo=='l')
        for (bt=0; bt<bn; bt++)
        {
            for (i=0; i<loc_A.row_num; i++)
                for (j=0; j<NB; j++)
                    U_R[i][j]=W_R[i][j]=0.0;
            for (i=0; i<loc_A.col_num; i++)
                for (j=0; j<NB; j++)
                    U_C[i][j]=W_C[i][j]=0.0;
            /* determain algorithm block size m*/
            first_col=bt*NB;
            end_col=(first_col+NB-1>N-2)?(N-2):(first_col+NB-1);
            m=end_col-first_col+1;
            indxg2p(comm2D,NB,first_col,first_col,&iarow,&iacol);/*diagonal element*/
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
                        L[i*loc_A.col_num+j]=a[i*loc_A.col_num+j];
                    if (tmp1>tmp2)
                        L1[i*loc_A.col_num+j]=a[i*loc_A.col_num+j];
                }
            }
            /*parallel translation in the panel block*/
            for (k=first_col; k<=end_col; k++)
            {
                for (i=k; i<N; i++)
                    x[i]=p[i]=u[i]=0.0;
                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                    x_r[i]=u_r[i]=w_r[i]=0.0;
                for (i=loc_A.col_pos; i<loc_A.col_num; i++)
                    x_c[i]=u_c[i]=w_c[i]=0.0;
                indxg2p(comm2D,NB,k,k,&iarow,&iacol);
                if (iarow==coord[0]) loc_A.row_pos++;
                /*generate u, tau*/
                if (coord[1]==iacol)
                {
                    for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                        u_r[i]=a[i*loc_A.col_num+loc_A.col_pos];
                    MPI_Bcast(&u_r[loc_A.row_pos],loc_A.row_num-loc_A.row_pos,MPI_DOUBLE,iacol,comm_row);
                }
                else
                    MPI_Bcast(&u_r[loc_A.row_pos],loc_A.row_num-loc_A.row_pos,MPI_DOUBLE,iacol,comm_row);
                indxg2p(comm2D,NB,k+1,k,&iarow,&iacol);
                indxg2l(k+1,k,NB,dim[0],dim[1],&loc_i,&loc_j);
                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                {
                    tmp=loc_A.row_set[i];
                    p[tmp]=u_r[i];
                }
                MPI_Allreduce(&p[k+1],&u[k+1],N-k-1,MPI_DOUBLE,MPI_SUM,comm_col);
                alpha=u[k+1];
                lm=N-k-1;
                dlarfg_(&lm,&alpha,&u[k+2],&incx,&g);
                norm[k]=g;
                if (coord[0]==iarow)
                {
                    u_r[loc_i]=1.0;
                }
                z[k]=alpha;
                u[k+1]=1.0;
                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                {
                    tmp=loc_A.row_set[i];
                    u_r[i]=u[tmp];
                }
                if (coord[1]==iacol)
                    for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                        a[i*loc_A.col_num+loc_A.col_pos]=u_r[i];
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
                    x[i]=p[i]=0.0;

                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                {
                    tmp1=loc_A.row_set[i];
                    for (j=loc_A.col_pos; j<outj; j++)
                    {
                        tmp2=loc_A.col_set[j];
                        if (tmp1>=tmp2)
                            L[i*loc_A.col_num+j]=a[i*loc_A.col_num+j];
                        if (tmp1>tmp2)
                            L1[i*loc_A.col_num+j]=a[i*loc_A.col_num+j];
                    }
                }
                /*local matrix-std::vector production*/
                transa='t';
                lm=loc_A.row_num-loc_A.row_pos;
                ln=loc_A.col_num-loc_A.col_pos;
                beta=0.0;
                lda=loc_A.col_num;
                pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;
                incx=1;
                incy=1;
                dgemv_(&transa,&ln,&lm,&g,&L[pos],&lda,&u_c[loc_A.col_pos],&incx,&beta,
                       &x_r[loc_A.row_pos],&incy);
                transa='n';
                dgemv_(&transa,&ln,&lm,&g,&L1[pos],&lda,&u_r[loc_A.row_pos],&incx,&beta,
                       &x_c[loc_A.col_pos],&incy);
                /*correct for out of date entries of A*/
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
                    tu[i]=tw[i]=0.0;
                transa='n';
                lm=(loc_A.col_num-outj)>0?(loc_A.col_num-outj):0;
                ln=k-first_col;
                beta=0.0;
                lda=NB;
                incx=1;
                incy=1;
                dgemv_(&transa,&ln,&lm,&g,&W_C[outj][0],&lda,&u_c[outj],&incx,&beta,tw,&incy);
                dgemv_(&transa,&ln,&lm,&g,&U_C[outj][0],&lda,&u_c[outj],&incx,&beta,tu,&incy);
                for (i=outi; i<loc_A.row_num; i++)
                    x_r[i]=0.0;
                transa='t';
                lm=(loc_A.row_num-outi)>0?(loc_A.row_num-outi):0;
                ln=k-first_col;
                alpha=1.0;
                beta=1.0;
                lda=NB;
                incx=1;
                incy=1;
                dgemv_(&transa,&ln,&lm,&alpha,&W_R[outi][0],&lda,tu,&incx,&beta,&x_r[outi],&incy);
                dgemv_(&transa,&ln,&lm,&alpha,&U_R[outi][0],&lda,tw,&incx,&beta,&x_r[outi],&incy);
                for (i=outi; i<loc_A.row_num; i++)
                {
                    tmp=loc_A.row_set[i];
                    x[tmp]=x[tmp]-x_r[i];
                }
                MPI_Allreduce(&x[k+1],&q[k+1],N-k-1,MPI_DOUBLE,MPI_SUM,comm_col);
                MPI_Allreduce(&q[k+1],&p[k+1],N-k-1,MPI_DOUBLE,MPI_SUM,comm_row);
                s=0.0;
                lm=N-k-1;
                for (i=k+1; i<N; i++)
                    s+=p[i]*u[i];
                s=-s/2.0;
                alpha=s*g;
                /*update remainder of panel block of A*/
                for (i=loc_A.row_pos; i<loc_A.row_num; i++)
                {
                    tmp=loc_A.row_set[i];
                    w_r[i]=(p[tmp]+(alpha*u_r[i]));
                }
                for (i=loc_A.col_pos; i<loc_A.col_num; i++)
                {
                    tmp=loc_A.col_set[i];
                    w_c[i]=(p[tmp]+(alpha*u_c[i]));
                }
                lm=loc_A.row_num-loc_A.row_pos;
                ln=(outj-loc_A.col_pos)>0?(outj-loc_A.col_pos):0;
                alpha=-1.0;
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
                }
                dger_(&ln,&lm,&alpha,&u_c[loc_A.col_pos],&incx,&w_r[loc_A.row_pos],&incy,&a[pos],&lda);
                dger_(&ln,&lm,&alpha,&w_c[loc_A.col_pos],&incx,&u_r[loc_A.row_pos],&incy,&a[pos],&lda);
                MPI_Barrier(DIAG_HPSEPS_WORLD);
            }
            /*rank-2b update of the submatrix A22*/
            transa='t';
            transb='n';
            lm=loc_A.row_num-loc_A.row_pos;
            ln=loc_A.col_num-loc_A.col_pos;
            in=end_col-first_col+1;
            alpha=-1.0;
            beta=1.0;
            ldb=NB;
            ldc=loc_A.col_num;
            pos=loc_A.row_pos*loc_A.col_num+loc_A.col_pos;

            dgemm_(&transa,&transb,&ln,&lm,&in,&alpha,&W_C[loc_A.col_pos][0],&ldb,
                   &U_R[loc_A.row_pos][0],&ldb,&beta,&a[pos],&ldc);
            dgemm_(&transa,&transb,&ln,&lm,&in,&alpha,&U_C[loc_A.col_pos][0],&ldb,
                   &W_R[loc_A.row_pos][0],&ldb,&beta,&a[pos],&ldc);
            MPI_Barrier(DIAG_HPSEPS_WORLD);
        }
    for (i=0; i<N; i++)
        dg[i]=0.0;
    for (i=0; i<loc_A.row_num; i++)
        for (j=0; j<loc_A.col_num; j++)
        {
            tmp1=loc_A.row_set[i];
            tmp2=loc_A.col_set[j];
            if (tmp1==tmp2)
                dg[tmp1]=a[i*loc_A.col_num+j];
        }
    MPI_Allreduce(dg,diag,N,MPI_DOUBLE,MPI_SUM,comm_row);
    MPI_Allreduce(diag,dg,N,MPI_DOUBLE,MPI_SUM,comm_col);
    BlasConnector::copy(N,dg,indx,diag,indy);
    BlasConnector::copy(N,z,indx,off_diag,indy);
    off_diag[N-1]=0.0;
    MPI_Comm_free(&comm_row);
    MPI_Comm_free(&comm_col);
    free(L);
    free(L1);
}









