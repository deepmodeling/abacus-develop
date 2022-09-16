#pragma once
#include "mpi.h"

#include <complex>
#include <fstream>
#include <vector>

class ELPA_Solver
{
  public:
    ELPA_Solver(const bool isReal,
                const MPI_Comm comm,
                const int nev,
                const int narows,
                const int nacols,
                const int* desc);
    ELPA_Solver(const bool isReal,
                const MPI_Comm comm,
                const int nev,
                const int narows,
                const int nacols,
                const int* desc,
                const int* otherParameter);

    int eigenvector(double* A, double* EigenValue, double* EigenVector);
    int generalized_eigenvector(double* A, double* B, int& DecomposedState, double* EigenValue, double* EigenVector);
    int eigenvector(std::complex<double>* A, double* EigenValue, std::complex<double>* EigenVector);
    int generalized_eigenvector(std::complex<double>* A,
                                std::complex<double>* B,
                                int& DecomposedState,
                                double* EigenValue,
                                std::complex<double>* EigenVector);
    void setLoglevel(int loglevel);
    void setKernel(bool isReal, int Kernel);
    void setQR(int useQR);
    void outputParameters();
    void verify(double* A, double* EigenValue, double* EigenVector, double& maxRemain, double& meanRemain);
    void verify(double* A, double* B, double* EigenValue, double* EigenVector, double& maxRemain, double& meanRemain);
    void verify(std::complex<double>* A,
                double* EigenValue,
                std::complex<double>* EigenVector,
                double& maxError,
                double& meanError);
    void verify(std::complex<double>* A,
                std::complex<double>* B,
                double* EigenValue,
                std::complex<double>* EigenVector,
                double& maxError,
                double& meanError);
    void exit();

  private:
    const int CHOLESKY_CRITICAL_SIZE = 1000;
    bool isReal;
    MPI_Comm comm;
    int nFull;
    int nev;
    int narows;
    int nacols;
    int desc[9];
    int method;
    int kernel_id;
    int cblacs_ctxt;
    int nblk;
    int lda;
    std::vector<double> dwork;
    std::vector<std::complex<double>> zwork;
    int myid;
    int nprows;
    int npcols;
    int myprow;
    int mypcol;
    int useQR;
    int wantDebug;
    int loglevel;
    std::ofstream logfile;
    // for legacy interface
    int comm_f;
    int mpi_comm_rows;
    int mpi_comm_cols;
    // for new elpa handle
    int handle_id;

    // toolbox
    int read_cpuflag();
    int read_real_kernel();
    int read_complex_kernel();
    int allocate_work();
    int decomposeRightMatrix(double* B, double* EigenValue, double* EigenVector, int& DecomposedState);
    int decomposeRightMatrix(std::complex<double>* B,
                             double* EigenValue,
                             std::complex<double>* EigenVector,
                             int& DecomposedState);
    int composeEigenVector(int DecomposedState, double* B, double* EigenVector);
    int composeEigenVector(int DecomposedState, std::complex<double>* B, std::complex<double>* EigenVector);
    // debug tool
    void timer(int myid, const char function[], const char step[], double& t0);
};
