# ifndef DIAGO_CG_CUH
# define DIAGO_CG_CUH

#ifdef __CUDA
#include "cufft.h"
#include "cublas_v2.h"
#endif

template<class T, class T2>
class Diago_CG_CUDA
{
public:

    Diago_CG_CUDA();
    ~Diago_CG_CUDA();

    static int moved;

    // Calculates the inner product of two float vectors, dim = 2*len.
    float ddot_real(
        const int & dim,
        const float2* psi_L,
        const float2* psi_R,
        const bool reduce = true) ;
    
    // Calculates the inner product of two double vectors, dim = 2*len.
    double ddot_real(
        const int & dim,
        const double2* psi_L,
        const double2* psi_R,
        const bool reduce = true) ;
    
    // Calculates the inner product of two float vectors on device, dim = 2*len.
    void ddot_real_on_device(
        const int &dim,
        const float2* psi_L,
        const float2* psi_R,
        float *res);
    
    // Calculates the inner product of two double vectors on device, dim = 2*len.
    void ddot_real_on_device(
        const int &dim,
        const double2* psi_L,
        const double2* psi_R,
        double *res);
    
    // Calculates the inner product of two complex float vectors, dim = len.
    float2 ddot(
        const int & dim,
        const float2* psi_L,
        const float2* psi_R ) ;
    
    // Calculates the inner product of two complex double vectors, dim = len.
    double2 ddot(
        const int & dim,
        const double2* psi_L,
        const double2* psi_R ) ;

    // Calculates the inner product of complex float vectors psi[m] and psik.
    float2 ddot(
        const int & dim,
        const float2 *psi, // matrix
        const int & m,
        float2 *psik ) ;

    // Calculates the inner product of complex double vectors psi[m] and psik.
    double2 ddot(
        const int & dim,
        const double2 *psi, // matrix
        const int & m,
        double2 *psik ) ;

    // Calculates the inner product of complex float vectors psi_L[m] and psi_R[n].
    float2 ddot(
        const int & dim,
        const float2 *psi_L, // matrix
        const int & m,
        const float2 *psi_R, // matrix
        const int & n) ;
    
    // Calculates the inner product of complex double vectors psi_L[m] and psi_R[n].
    double2 ddot(
        const int & dim,
        const double2 *psi_L, // matrix
        const int & m,
        const double2 *psi_R, // matrix
        const int & n) ;
    
    // Solving eigenvalues using the Conjugate Gradient Method. 
    void diag(
        T2 *phi, // matrix
        T *e,
        T2 *vkb_c,
        const int &dim,
        const int &dmx,
        const int &n_band,
        const T *precondition,
        const T &eps,
        const int &maxter,
        const bool &reorder,
        int &notconv,
        double &avg_iter);

    // Schmidt Orthogonalization (FP32)
    void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &m,
        const float2 *psi, // matrix
        float2 *spsi,
        float2 *psi_m
    );

    // Schmidt Orthogonalization (FP64)
    void schmit_orth(
        const int &dim,
        const int &dmx,
        const int &m,
        const double2 *psi, // matrix
        double2 *spsi,
        double2 *psi_m
    );

private:

    int test_cg;
#ifdef __CUDA
    cublasHandle_t diag_handle; // cublas handle
    // cublasHandle_t diag_handle_device; // pointer_mode = device
    T* tmp_3value;
#endif

#ifdef __ROCM
    hipblasHandle_t diag_handle; // hipblas handle
#endif
    // Calculate the gardient during each CG Loop.
    void calculate_gradient(
        const T* precondition,
        const int dim,
        const T2 *hpsi,
        const T2 *spsi,
        T2 *g,
        T2 *pspsi);
    
    // Schmidt Orthogonalization of the gradient(FP32).
    void orthogonal_gradient(
        const int &dim,
        const int &dmx,
        float2 *g,
        float2 *sg,
        float2 *lagrange,
        const float2 *eigenfunction, // matrix
        const int m);
    
    // Schmidt Orthogonalization of the gradient(FP64).
    void orthogonal_gradient(
        const int &dim,
        const int &dmx,
        double2 *g,
        double2 *sg,
        double2 *lagrange,
        const double2 *eigenfunction, // matrix
        const int m);

    // calculate parameter alpha and beta during each CG Loop.
    void calculate_gamma_cg(
        const int iter,
        const int dim,
        const T *precondition,
        const T2 *g,
        const T2 *sg,
        T2 *psg,
        T2 *cg,
        T &gg_last,
        const T &cg0,
        const T &theta,
        const T2 *psi_m);

    // Update psi in each CG Loop.
    bool update_psi(
        const int dim,
        T &cg_norm,
        T &theta,
        T2 *hcg,
        const T2 *cg,
        T2 *scg,
        T2 *psi_m ,
        T &eigenvalue,
        const T &threshold,
        T2 *hpsi,
        T2 *spsi,
        T2 *vkb_c);

};
# endif
