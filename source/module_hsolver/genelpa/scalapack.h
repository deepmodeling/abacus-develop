#pragma once
// scalapack
int numroc_(const int *N, const int *NB, const int *IPROC, const int *ISRCPROC, const int *NPROCS);
void descinit_(int *DESC, const int *M, const int *N, const int *MB, const int *NB, const int *IRSRC, const int *ICSRC, const int *ICTXT, const int *LLD, int *INFO);
void pdpotrf_(const char *UPLO, const int *N, double *A, const int *IA, const int *JA, const int *DESCA, int *INFO);
void pzpotrf_(const char *UPLO, const int *N, double _Complex *A, const int *IA, const int *JA, const int *DESCA, int *INFO);
void pdsyev_(const char *JOBZ, const char *UPLO, int *N, double *A, int *IA, int *JA, int *DESCA, 
             double *W, double *Z, int *IZ, int *JZ, int *DESCZ, double *WORK, int *LWORK, int *INFO);                          
void pdgemr2d_(int *M, int *N, double *A, int *IA, int *JA, int *DESCA, 
			   double *B, int *IB, int *JB, int *DESCB, int *ICTXT);			   
void pzgemr2d_(int *M, int *N, double _Complex *A, int *IA, int *JA, int *DESCA, 
			   double _Complex *B, int *IB, int *JB, int *DESCB, int *ICTXT);
