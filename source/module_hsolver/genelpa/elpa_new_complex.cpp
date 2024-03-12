#include <complex>
#include <map>
#include <regex>
#include <fstream>
#include <cfloat>
#include <cstring>
#include <mpi.h>

#include "elpa_new.h"
#include "elpa_solver.h"

#include "my_math.hpp"
#include "utils.h"

extern std::map<int, elpa_t> NEW_ELPA_HANDLE_POOL;

int ELPA_Solver::eigenvector(std::complex<double>* A, double* EigenValue, std::complex<double>* EigenVector)
{
    int info;
    int allinfo;
    double t;

    if((loglevel>0 && myid==0) || loglevel>1)
    {
        t=-1;
        timer(myid, "elpa_eigenvectors_dc", "1", t);
    }
    elpa_eigenvectors(NEW_ELPA_HANDLE_POOL[handle_id], A, EigenValue, EigenVector, &info);
    if((loglevel>0 && myid==0) || loglevel>1)
    {
        timer(myid, "elpa_eigenvectors_dc", "1", t);
    }
    MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
    allinfo = info;
    return allinfo;
}

int ELPA_Solver::generalized_eigenvector(std::complex<double>* A, std::complex<double>* B, int& DecomposedState,
                                         double* EigenValue, std::complex<double>* EigenVector)
{
    int info, allinfo;
    double t;

    if((loglevel>0 && myid==0) || loglevel>1)
    {
        t=-1;
        timer(myid, "decomposeRightMatrix", "1", t);
    }
    if(DecomposedState==0) // B is not decomposed
        allinfo=decomposeRightMatrix(B, EigenValue, EigenVector, DecomposedState);
    else
        allinfo=0;

    if((loglevel>0 && myid==0) || loglevel>1)
    {
        timer(myid, "decomposeRightMatrix", "1", t);
    }
    if(allinfo != 0)
        return allinfo;

    // transform A to A~
    if((loglevel>0 && myid==0) || loglevel>1)
    {
        t=-1;
        timer(myid, "transform A to A~", "2", t);
    }
    if(DecomposedState == 1 || DecomposedState == 2)
    {
        // calculate A*U^-1, put to work
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "A*U^-1", "2.1a", t);
        }
        Cpzgemm('C', 'N', nFull, 1.0, A, B, 0.0, zwork.data(), desc);
        if(loglevel>1)
        {
            timer(myid, "A*U^-1", "2.1a", t);
        }

        // calculate U^-C^(A*U^-1), put to a
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "U^-T*(A*U^-1)", "2.2a", t);
        }
        Cpzgemm('C', 'N', nFull, 1.0, B, zwork.data(), 0.0, A, desc);
        if(loglevel>1)
        {
            timer(myid, "U^-T*(A*U^-1)", "2.2a", t);
        }
    }
    else
    {
        // calculate b*a^C and put to work
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "B*A^T", "2.1b", t);
        }
        Cpzgemm('N', 'C', nFull, 1.0, B, A, 0.0, zwork.data(), desc);
        if(loglevel>1)
        {
            timer(myid, "B*A^T", "2.1b", t);
        }
        // calculate b*work^C and put to a -- original A*x=v*B*x was transform to a*x'=v*x'
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "B*(B*A^T)^T", "2.2b", t);
        }
        Cpzgemm('N', 'C', nFull, 1.0, B, zwork.data(), 0.0, A, desc);
        if(loglevel>1)
        {
            timer(myid, "B*(B*A^T)^T", "2.2b", t);
        }
    }
    if((loglevel>0 && myid==0) || loglevel>1)
    {
        timer(myid, "transform A to A~", "2", t);
    }

    // calculate the eigenvalue and eigenvector of A~
    if((loglevel>0 && myid==0) || loglevel>1)
    {
        t=-1;
        timer(myid, "elpa_eigenvectors", "3", t);
    }
    if(loglevel>2) saveMatrix("A_tilde.dat", nFull, A, desc, cblacs_ctxt);
    info=eigenvector(A, EigenValue, EigenVector);
    if((loglevel>0 && myid==0) || loglevel>1)
    {
        timer(myid, "elpa_eigenvectors", "3", t);
    }
#ifdef __MPI
    MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
#else
    allinfo = info;
#endif
    if(loglevel>2) saveMatrix("EigenVector_tilde.dat", nFull, EigenVector, desc, cblacs_ctxt);

    if((loglevel>0 && myid==0) || loglevel>1)
    {
        t=-1;
        timer(myid, "composeEigenVector", "4", t);
    }
    // transform eigenvector c~ to original eigenvector c
    allinfo=composeEigenVector(DecomposedState, B, EigenVector);
    if((loglevel>0 && myid==0) || loglevel>1)
    {
        timer(myid, "composeEigenVector", "4", t);
    }
    return allinfo;
}

int ELPA_Solver::decomposeRightMatrix(std::complex<double>* B, double* EigenValue, std::complex<double>* EigenVector, int& DecomposedState)
{
    int info=0;
    int allinfo=0;
    double t;

    // first try cholesky decomposing
    if(nFull<CHOLESKY_CRITICAL_SIZE) // use pdpotrf for small matrix
    {
        DecomposedState=1;
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "pzpotrf_", "1", t);
        }
        Cpzpotrf('U', nFull, B, desc);
        if(loglevel>1)
        {
            timer(myid, "pzpotrf_", "1", t);
        }
        MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
        allinfo = info;
        if(allinfo != 0) //if pdpotrf fail, try elpa_cholesky_real
        {
            DecomposedState=2;
            if(loglevel>1)
            {
                t=-1;
                timer(myid, "elpa_cholesky_dc", "2", t);
            }
            elpa_cholesky(NEW_ELPA_HANDLE_POOL[handle_id], B, &info);
            if(loglevel>1)
            {
                timer(myid, "elpa_cholesky_dc", "2", t);
            }
            MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
            allinfo = info;
        }
    } else
    {
        DecomposedState=2;
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "elpa_cholesky_dc", "1", t);
        }
        elpa_cholesky(NEW_ELPA_HANDLE_POOL[handle_id], B, &info);
        if(loglevel>1)
        {
            timer(myid, "elpa_cholesky_dc", "1", t);
        }
        MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
        allinfo = info;
        if(allinfo != 0)
        {
            DecomposedState=1;
            if(loglevel>1)
            {
                t=-1;
                timer(myid, "pzpotrf_", "2", t);
            }
            Cpzpotrf('U', nFull, B, desc);
            if(loglevel>1)
            {
                timer(myid, "pzpotrf_", "2", t);
            }
            MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
            allinfo = info;
        }
    }

    if(allinfo==0) // calculate U^{-1}
    {
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "clear low triangle", "1", t);
        }
        for(int j=0; j<nacols; ++j)
        {
            int jGlobal=globalIndex(j, nblk, npcols, mypcol);
            for(int i=0; i<narows; ++i)
            {
                int iGlobal=globalIndex(i, nblk, nprows, myprow);
                if(iGlobal>jGlobal) B[i+j*narows]=0;
            }
        }
        if(loglevel>1)
        {
            timer(myid, "clear low triangle", "1", t);
        }
        if(loglevel>2) saveMatrix("U.dat", nFull, B, desc, cblacs_ctxt);
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "invert U", "1", t);
        }
        elpa_invert_triangular(NEW_ELPA_HANDLE_POOL[handle_id], B, &info);
        if(loglevel>1)
        {
            timer(myid, "invert U", "1", t);
        }
        if(loglevel>2) saveMatrix("U_inv.dat", nFull, B, desc, cblacs_ctxt);
    } else {
        // if cholesky decomposing failed, try diagonalize
        // calculate B^{-1/2}_{i,j}=\sum_k q_{i,k}*ev_k^{-1/2}*q_{j,k} and put to b, which will be b^-1/2
        DecomposedState=3;
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "calculate eigenvalue and eigenvector of B", "1", t);
        }
        //elpa_eigenvectors_all_host_arrays_dc(NEW_ELPA_HANDLE_POOL[handle_id], b,
        //                     EigenValue, q, &info);
        info=eigenvector(B, EigenValue, EigenVector);
        if(loglevel>1)
        {
            timer(myid, "calculate eigenvalue and eigenvector of B", "1", t);
        }
        MPI_Allreduce(&info, &allinfo, 1, MPI_INT, MPI_MAX, comm);
        allinfo = info;
        // calculate q*ev and put to work
        for(int i=0; i<nacols; ++i)
        {
            int eidx=globalIndex(i, nblk, npcols, mypcol);
            //double ev_sqrt=1.0/sqrt(ev[eidx]);
            double ev_sqrt=EigenValue[eidx]>DBL_MIN?1.0/sqrt(EigenValue[eidx]):0;
            for(int j=0; j<narows; ++j)
                zwork[i*lda+j]=EigenVector[i*lda+j]*ev_sqrt;
        }

        // calculate qevq=qev*q^T, put to b, which is B^{-1/2}
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "qevq=qev*q^T", "2", t);
        }
        Cpzgemm('N', 'C', nFull, 1.0, zwork.data(), EigenVector, 0.0, B, desc);
        if(loglevel>1)
        {
            timer(myid, "qevq=qev*q^T", "2", t);
        }
    }
    return allinfo;
}

int ELPA_Solver::composeEigenVector(int DecomposedState, std::complex<double>* B, std::complex<double>* EigenVector)
{
    double t;
    if(DecomposedState==1 || DecomposedState==2)
    {
        // transform the eigenvectors to original general equation, let U^-1*q, and put to q
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "Cpztrmm", "1", t);
        }
        Cpztrmm('L', 'U', 'N', 'N', nFull, nev, 1.0, B, EigenVector, desc);
        if(loglevel>1)
        {
            timer(myid, "Cpztrmm", "1", t);
        }
    } else {
        // transform the eigenvectors to original general equation, let b^C*q, and put to q
        if(loglevel>1)
        {
            t=-1;
            timer(myid, "Cpzgemm", "1", t);
        }
        Cpzgemm('C', 'N', nFull, nev, nFull, 1.0, B, zwork.data(), 0.0, EigenVector, desc);
        if(loglevel>1)
        {
            timer(myid, "Cpzgemm", "1", t);
        }
    }
    return 0;
}

// calculate the error
// $ \ket{ \delta \psi_i } = H\ket{\psi_i} $
// $ \delta_i = \braket{ \delta \psi_i | \delta \psi_i } $
//
// V: eigenvector matrix
// D: Diagonal matrix of eigenvalue
// maxError: maximum absolute value of error
// meanError: mean absolute value of error
void ELPA_Solver::verify(std::complex<double>* A, double* EigenValue, std::complex<double>* EigenVector,
                         double &maxError, double &meanError)
{
    std::complex<double>* V=EigenVector;
    const int naloc=narows*nacols;
    std::complex<double>* D=new std::complex<double>[naloc];
    std::complex<double>* R=zwork.data();

    for(int i=0; i<naloc; ++i)
        D[i]=0;

    for(int i=0; i<nFull; ++i)
    {
        int localRow, localCol;
        int localProcRow, localProcCol;

        localRow=localIndex(i, nblk, nprows, localProcRow);
        if(myprow==localProcRow)
        {
            localCol=localIndex(i, nblk, npcols, localProcCol);
            if(mypcol==localProcCol)
            {
                int idx = localRow + localCol*narows;
                D[idx]=EigenValue[i];
            }
        }
    }

    // R=V*D
    Cpzhemm('R', 'U', nFull, 1.0, D, V, 0.0, R, desc);
    if(loglevel>2) saveMatrix("VD.dat", nFull, R, desc, cblacs_ctxt);
    // R=A*V-V*D=A*V-R
    Cpzhemm('L', 'U', nFull, 1.0, A, V, -1.0, R, desc);
    if(loglevel>2) saveMatrix("AV-VD.dat", nFull, R, desc, cblacs_ctxt);
    // calculate the maximum and mean value of sum_i{R(:,i)*R(:,i)}
    double sumError=0;
    maxError=0;
    for(int i=1; i<=nev; ++i)
    {
        std::complex<double> E;
        Cpzdotc(nFull, E, R, 1, i, 1,
                         R, 1, i, 1, desc);
        double abs_E=std::abs(E);
        sumError+=abs_E;
        maxError=std::max(maxError, abs_E);
    }
    meanError=sumError/nFull;
    delete[] D;
}

// calculate the error
// $ \ket{ \delta \psi_i } = (H - \epsilon_i S)\ket{\psi_i} $
// $ \delta_i = \braket{ \delta \psi_i | \delta \psi_i } $
//
// V: eigenvector matrix
// D: Diagonal matrix of eigenvalue
// maxError: maximum absolute value of error
// meanError: mean absolute value of error
void ELPA_Solver::verify(std::complex<double>* A, std::complex<double>* B,
                        double* EigenValue, std::complex<double>* EigenVector,
                        double &maxError, double &meanError)
{
    std::complex<double>* V=EigenVector;
    const int naloc=narows*nacols;
    std::complex<double>* D=new std::complex<double>[naloc];
    std::complex<double>* R=new std::complex<double>[naloc];

    for(int i=0; i<naloc; ++i)
        D[i]=0;

    for(int i=0; i<nFull; ++i)
    {
        int localRow, localCol;
        int localProcRow, localProcCol;

        localRow=localIndex(i, nblk, nprows, localProcRow);
        if(myprow==localProcRow)
        {
            localCol=localIndex(i, nblk, npcols, localProcCol);
            if(mypcol==localProcCol)
            {
                int idx = localRow + localCol*narows;
                D[idx]=EigenValue[i];
            }
        }
    }

    // zwork=B*V
    Cpzhemm('L', 'U', nFull, 1.0, B, V, 0.0, zwork.data(), desc);
    if(loglevel>2) saveMatrix("BV.dat", nFull, zwork.data(), desc, cblacs_ctxt);
    // R=B*V*D=zwork*D
    Cpzhemm('R', 'U', nFull, 1.0, D, zwork.data(), 0.0, R, desc);
    if(loglevel>2) saveMatrix("BVD.dat", nFull, R, desc, cblacs_ctxt);
    // R=A*V-B*V*D=A*V-R
    Cpzhemm('L', 'U', nFull, 1.0, A, V, -1.0, R, desc);
    if(loglevel>2) saveMatrix("AV-BVD.dat", nFull, R, desc, cblacs_ctxt);
    // calculate the maximum and mean value of sum_i{R(:,i)*R(:,i)}
    double sumError=0;
    maxError=0;
    for(int i=1; i<=nev; ++i)
    {
        std::complex<double> E;
        Cpzdotc(nFull, E, R, 1, i, 1,
                         R, 1, i, 1, desc);
        double abs_E=std::abs(E);
        sumError+=abs_E;
        maxError=std::max(maxError, abs_E);
    }
    meanError=sumError/nFull;

    delete[] D;
    delete[] R;
}
