#pragma once
void pddot_(int* n, double* dot, double* x, int* ix, int* jx, int* descx, int* incx,
								 double* y, int* iy, int* jy, int* descy, int* incy);
								 
void pzdotc_(int* n, double _Complex* dot, double _Complex* x, int* ix, int* jx, int* descx, int* incx,
								 double _Complex* y, int* iy, int* jy, int* descy, int* incy);
void pdsymv_(char* uplo, int* n, 
			 double* alpha, double* a, int* ia, int* ja, int* desca,
						    double* x, int* ix, int* jx, int* descx, int* incx,
             double* beta,  double* y, int* iy, int* jy, int* descy, int* incy);
void pdtran_(int* m , int* n ,
             double* alpha , double* a , int* ia , int* ja , int* desca ,
             double* beta ,  double* c , int* ic , int* jc , int* descc );

void pdgemm_(char* transa , char* transb , int* m , int* n , int* k ,
             double* alpha , double* a , int* ia , int* ja , int* desca ,
                             double* b , int* ib , int* jb , int* descb ,
             double* beta ,  double* c , int* ic , int* jc , int* descc );
void pzgemm_(char* transa , char* transb , int* m , int* n , int* k ,
             double _Complex* alpha , double _Complex* a , int* ia , int* ja , int* desca ,
									  double _Complex* b , int* ib , int* jb , int* descb ,
             double _Complex* beta ,  double _Complex* c , int* ic , int* jc , int* descc );
void pdsymm_(char* side , char* uplo , int* m , int* n ,
             double* alpha , double* a , int* ia , int* ja , int* desca ,
                             double* b , int* ib , int* jb , int* descb ,
             double* beta ,  double* c , int* ic , int* jc , int* descc );
void pzsymm_(char* side , char* uplo , int* m , int* n ,
             double _Complex* alpha , double _Complex* a , int* ia , int* ja , int* desca ,
									  double _Complex* b , int* ib , int* jb , int* descb ,
             double _Complex* beta ,  double _Complex* c , int* ic , int* jc , int* descc );
void pzhemm_(char* side , char* uplo , int* m , int* n ,
             double _Complex* alpha , double _Complex* a , int* ia , int* ja , int* desca ,
									  double _Complex* b , int* ib , int* jb , int* descb ,
             double _Complex* beta ,  double _Complex* c , int* ic , int* jc , int* descc );
void pdtrmm_(char* side , char* uplo , char* transa , char* diag , int* m , int* n ,
             double* alpha , double* a , int* ia , int* ja , int* desca ,
                             double* b , int* ib , int* jb , int* descb );
void pztrmm_(char* side , char* uplo , char* transa , char* diag , int* m , int* n ,
             double _Complex* alpha ,  double _Complex* a , int* ia , int* ja , int* desca ,
									   double _Complex* b , int* ib , int* jb , int* descb );
