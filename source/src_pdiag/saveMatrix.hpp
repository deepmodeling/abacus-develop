#ifdef __MPI
#include "mpi.h"
#endif
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <complex>

inline void saveMatrix(char* filePrefix, char* LAYOUT, int narows, int nacols, double* a)
{
    using namespace std;
    char fileName[80];
    int myid;
    ofstream matrixFile;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    sprintf(fileName, "%s_%3.3d.dat", filePrefix, myid);
    matrixFile.open(fileName);
    matrixFile.flags(ios_base::scientific);
    matrixFile.precision(17);
    matrixFile.width(24);
    if(*LAYOUT=='R'||*LAYOUT=='r')
    {
        for(int i=0; i<narows; ++i)
        {
            for(int j=0; j<nacols; ++j)
            {
                matrixFile<<a[i*nacols+j]<<" ";
            }
            matrixFile<<endl;
        }
    }
    else
    {
        for(int i=0; i<narows; ++i)
        {
            for(int j=0; j<nacols; ++j)
            {
                matrixFile<<a[i+j*narows]<<" ";
            }
            matrixFile<<endl;
        }
    }
    matrixFile.close();
}

inline void saveMatrix(char* filePrefix, char* LAYOUT, int narows, int nacols, complex<double>* a)
{
    using namespace std;
    char fileName[80];
    int myid;
    ofstream matrixFile;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    sprintf(fileName, "%s_%3.3d.dat", filePrefix, myid);
    matrixFile.open(fileName);
    matrixFile.flags(ios_base::scientific);
    matrixFile.precision(17);
    matrixFile.width(24);
    if(*LAYOUT=='R'||*LAYOUT=='r')
    {
        for(int i=0; i<narows; ++i)
        {
            for(int j=0; j<nacols; ++j)
            {
                matrixFile<<a[i*nacols+j]<<" ";
            }
            matrixFile<<endl;
        }
    }
    else
    {
        for(int i=0; i<narows; ++i)
        {
            for(int j=0; j<nacols; ++j)
            {
                matrixFile<<a[i+j*narows]<<" ";
            }
            matrixFile<<endl;
        }
    }
    matrixFile.close();
}

