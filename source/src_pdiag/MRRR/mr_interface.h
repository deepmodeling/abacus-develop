#ifdef __MPI
#include "mpi.h"
#endif
extern int pdstemr_mpi(MPI_Comm comm1D, char *jobz, char *range, int *n,
        double * d__, double *e, double *vl, double *vu, int *il, int *iu,
        int *m, double *w, double *z__, int *ldz, int *nzc, int *isuppz,
        bool *tryrac, double *work, int *lwork, int *iwork, int *liwork,
        int *info);
