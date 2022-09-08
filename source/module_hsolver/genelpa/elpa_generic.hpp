// `elpa_generic.h` replacement for version 2021.11.002
#pragma once
static inline void elpa_set(elpa_t e, const char *name, int value, int *error)
{
	elpa_set_integer(e, name, value, error);
}

static inline void elpa_set(elpa_t e, const char *name, double value, int *error)
{
	elpa_set_double(e, name, value, error);
}

static inline void elpa_get(elpa_t e, const char *name, int *value, int *error)
{
	elpa_get_integer(e, name, value, error);
}

static inline void elpa_get(elpa_t e, const char *name, double *value, int *error)
{
	elpa_get_double(e, name, value, error);
}

static inline void elpa_eigenvectors(elpa_t handle, double *a, double *ev, double *q, int *error)
{
	elpa_eigenvectors_all_host_arrays_d(handle, a, ev, q, error);
}

static inline void elpa_eigenvectors(elpa_t handle, float *a, float *ev, float *q, int *error)
{
	elpa_eigenvectors_all_host_arrays_f(handle, a, ev, q, error);
}

static inline void elpa_eigenvectors(elpa_t handle, double complex *a, double *ev, double complex *q, int *error)
{
	elpa_eigenvectors_all_host_arrays_dc(handle, a, ev, q, error);
}

static inline void elpa_eigenvectors(elpa_t handle, float complex *a, float *ev, float complex *q, int *error)
{
	elpa_eigenvectors_all_host_arrays_fc(handle, a, ev, q, error);
}

static inline void elpa_eigenvectors_double(elpa_t handle, double *a, double *ev, double *q, int *error)
{
	elpa_eigenvectors_device_pointer_d(handle, a, ev, q, error);
}

static inline void elpa_eigenvectors_float(elpa_t handle, float *a, float *ev, float *q, int *error)
{
	elpa_eigenvectors_device_pointer_f(handle, a, ev, q, error);
}

static inline void elpa_eigenvectors_double_complex(elpa_t handle, double complex *a, double *ev, double complex *q, int *error)
{
	elpa_eigenvectors_device_pointer_dc(handle, a, ev, q, error);
}

static inline void elpa_eigenvectors_float_complex(elpa_t handle, float complex *a, float *ev, float complex *q, int *error)
{
	elpa_eigenvectors_device_pointer_fc(handle, a, ev, q, error);
}

static inline void elpa_skew_eigenvectors(elpa_t handle, double *a, double *ev, double *q, int *error)
{
	elpa_eigenvectors_all_host_arrays_d(handle, a, ev, q, error);
}

static inline void elpa_skew_eigenvectors(elpa_t handle, float *a, float *ev, float *q, int *error)
{
	elpa_eigenvectors_all_host_arrays_f(handle, a, ev, q, error);
}

static inline void elpa_skew_eigenvectors_double(elpa_t handle, double *a, double *ev, double *q, int *error)
{
	elpa_eigenvectors_device_pointer_d(handle, a, ev, q, error);
}

static inline void elpa_skew_eigenvectors_float(elpa_t handle, float *a, float *ev, float *q, int *error)
{
	elpa_eigenvectors_device_pointer_f(handle, a, ev, q, error);
}

static inline void elpa_generalized_eigenvectors(elpa_t handle, double *a, double *b, double *ev, double *q,
	int is_already_decomposed, int *error)
{
	elpa_generalized_eigenvectors_d(handle, a, b, ev, q, is_already_decomposed, error);
}

static inline void elpa_generalized_eigenvectors(elpa_t handle, float *a, float *b, float *ev, float *q,
 	int is_already_decomposed, int *error)
{
	elpa_generalized_eigenvectors_f(handle, a, b, ev, q, is_already_decomposed, error);
}

static inline void elpa_generalized_eigenvectors(elpa_t handle, double complex *a, double complex *b, 
	double *ev, double complex *q, int is_already_decomposed, int *error)
{
	elpa_generalized_eigenvectors_dc(handle, a, b, ev, q, is_already_decomposed, error);
}

static inline void elpa_generalized_eigenvectors(elpa_t handle, float complex *a, float complex *b, 
	float *ev, float complex *q, int is_already_decomposed, int *error)
{
	elpa_generalized_eigenvectors_fc(handle, a, b, ev, q, is_already_decomposed, error);
}

static inline void elpa_eigenvalues(elpa_t handle, double *a, double *ev, int *error)
{
	elpa_eigenvalues_all_host_arrays_d(handle, a, ev, error);
}

static inline void elpa_eigenvalues(elpa_t handle, float *a, float *ev, int *error)
{
	elpa_eigenvalues_all_host_arrays_f(handle, a, ev, error);
}

static inline void elpa_eigenvalues(elpa_t handle, double complex *a, double *ev, int *error)
{
	elpa_eigenvalues_all_host_arrays_dc(handle, a, ev, error);
}

static inline void elpa_eigenvalues(elpa_t handle, float complex *a, float *ev, int *error)
{
	elpa_eigenvalues_all_host_arrays_fc(handle, a, ev, error);
}

static inline void elpa_eigenvalues_double(elpa_t handle, double *a, double *ev, int *error)
{
	elpa_eigenvalues_device_pointer_d(handle, a, ev, error);
}

static inline void elpa_eigenvalues_float(elpa_t handle, float *a, float *ev, int *error)
{
	elpa_eigenvalues_device_pointer_f(handle, a, ev, error);
}

static inline void elpa_eigenvalues_double_complex(elpa_t handle, double complex *a, double *ev, int *error)
{
	elpa_eigenvalues_device_pointer_dc(handle, a, ev, error);
}

static inline void elpa_eigenvalues_float_complex(elpa_t handle, float complex *a, float *ev, int *error)
{
	elpa_eigenvalues_device_pointer_fc(handle, a, ev, error);
}

static inline void elpa_skew_eigenvalues(elpa_t handle, double *a, double *ev, int *error)
{
	elpa_eigenvalues_all_host_arrays_d(handle, a, ev, error);
}

static inline void elpa_skew_eigenvalues(elpa_t handle, float *a, float *ev, int *error)
{
	elpa_eigenvalues_all_host_arrays_f(handle, a, ev, error);
}

static inline void elpa_skew_eigenvalues_double(elpa_t handle, double *a, double *ev, int *error)
{
	elpa_eigenvalues_device_pointer_d(handle, a, ev, error);
}

static inline void elpa_skew_eigenvalues_float(elpa_t handle, float *a, float *ev, int *error)
{
	elpa_eigenvalues_device_pointer_f(handle, a, ev, error);
}

static inline void elpa_cholesky(elpa_t handle, double *a, int *error)
{
	elpa_cholesky_d(handle, a, error);
}

static inline void elpa_cholesky(elpa_t handle, float *a, int *error)
{
	elpa_cholesky_f(handle, a, error);
}

static inline void elpa_cholesky(elpa_t handle, double complex *a, int *error)
{
	elpa_cholesky_dc(handle, a, error);
}

static inline void elpa_cholesky(elpa_t handle, float complex *a, int *error)
{
	elpa_cholesky_fc(handle, a, error);
}

static inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, 
	double *a, double *b, int nrows_b, int ncols_b, double *c, int nrows_c, int ncols_c, int *error)
{
	elpa_hermitian_multiply_d(handle, uplo_a, uplo_c, ncb,
					a, b, nrows_b, ncols_b, c, nrows_c, ncols_c, error);
}

static inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, 
	float *a, float *b, int nrows_b, int ncols_b, float *c, int nrows_c, int ncols_c, int *error)
{
	elpa_hermitian_multiply_df(handle, uplo_a, uplo_c, ncb,
					a, b, nrows_b, ncols_b, c, nrows_c, ncols_c, error);
}

static inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, 
	double complex *a, double complex *b, int nrows_b, int ncols_b, double complex *c, int nrows_c, int ncols_c, int *error)
{
	elpa_hermitian_multiply_dc(handle, uplo_a, uplo_c, ncb,
					a, b, nrows_b, ncols_b, c, nrows_c, ncols_c, error);
}

static inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, 
	float complex *a, float complex *b, int nrows_b, int ncols_b, float complex *c, int nrows_c, int ncols_c, int *error)
{
	elpa_hermitian_multiply_fc(handle, uplo_a, uplo_c, ncb,
					a, b, nrows_b, ncols_b, c, nrows_c, ncols_c, error);
}

static inline void elpa_invert_triangular(elpa_t handle, double *a, int *error)
{
	elpa_invert_trm_d(handle, a, error);
}

static inline void elpa_invert_triangular(elpa_t handle, float *a, int *error)
{
	elpa_invert_trm_f(handle, a, error);
}

static inline void elpa_invert_triangular(elpa_t handle, double complex *a, int *error)
{
	elpa_invert_trm_dc(handle, a, error);
}

static inline void elpa_invert_triangular(elpa_t handle, float complex *a, int *error)
{
	elpa_invert_trm_fc(handle, a, error);
}
