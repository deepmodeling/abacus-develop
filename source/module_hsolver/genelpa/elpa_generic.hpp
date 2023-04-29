#pragma once
#include "elpa_new.h"
#include <complex>
/*! \brief generic C method for elpa_set
 *
 *  \details
 *  \param  handle  handle of the ELPA object for which a key/value pair should be set
 *  \param  name    the name of the key
 *  \param  value   integer/double value to be set for the key
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
inline void elpa_set(elpa_t handle, const char *name, int value, int *error)
{
    elpa_set_integer(handle, name, value, error);
}
inline void elpa_set(elpa_t handle, const char *name, double value, int *error)
{
    elpa_set_double(handle, name, value, error);
}

/*! \brief generic C method for elpa_get
 *
 *  \details
 *  \param  handle  handle of the ELPA object for which a key/value pair should be queried
 *  \param  name    the name of the key
 *  \param  value   integer/double value to be queried
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
inline void elpa_get(elpa_t handle, const char *name, int *value, int *error)
{
    elpa_get_integer(handle, name, value, error);
}
inline void elpa_get(elpa_t handle, const char *name, double *value, int *error)
{
    elpa_get_double(handle, name, value, error);
}

/*! \brief generic C method for elpa_eigenvectors
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a
 *  \param  ev      on return: float/double pointer to eigenvalues
 *  \param  q       on return: float/double float complex/double complex pointer to eigenvectors
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
#if ELPA_API_VERSION <= 20210502 // ELPA 2021.05.002 and earlier versions
inline void elpa_eigenvectors(const elpa_t handle, double *a, double *ev, double *q, int *error)
{
    elpa_eigenvectors_d(handle, a, ev, q, error);
}

inline void elpa_eigenvectors(const elpa_t handle, float  *a, float  *ev, float  *q, int *error)
{
    elpa_eigenvectors_f(handle, a, ev, q, error);
}

inline void elpa_eigenvectors(const elpa_t handle, std::complex<double> *a, double *ev, std::complex<double> *q, int *error)
{
    elpa_eigenvectors_dc(handle, reinterpret_cast<double _Complex*>(a), ev, reinterpret_cast<double _Complex*>(q), error);
}

inline void elpa_eigenvectors(const elpa_t handle, std::complex<float>  *a, float  *ev, std::complex<float>  *q, int *error)
{
    elpa_eigenvectors_fc(handle, reinterpret_cast<float _Complex*>(a), ev, reinterpret_cast<float _Complex*>(q), error);
}
#elif ELPA_API_VERSION < 20220501   // ELPA version between 2021.11.001 and 2022.05.001
inline void elpa_eigenvectors(const elpa_t handle, double *a, double *ev, double *q, int *error)
{
    elpa_eigenvectors_all_host_arrays_d(handle, a, ev, q, error);
}

inline void elpa_eigenvectors(const elpa_t handle, float *a, float *ev, float *q, int *error)
{
    elpa_eigenvectors_all_host_arrays_f(handle, a, ev, q, error);
}

inline void elpa_eigenvectors(const elpa_t handle, std::complex<double> *a, double *ev, std::complex<double> *q, int *error)
{
    elpa_eigenvectors_all_host_arrays_dc(handle, reinterpret_cast<double _Complex*>(a),
                                         ev, reinterpret_cast<double _Complex*>(q), error);
}

inline void elpa_eigenvectors(const elpa_t handle, std::complex<float>  *a, float  *ev, std::complex<float>  *q, int *error)
{
    elpa_eigenvectors_all_host_arrays_fc(handle, reinterpret_cast<float _Complex*>(a),
                                         ev, reinterpret_cast<float _Complex*>(q), error);
}
#else // ELPA version 2022.05.001, ELPA has its own c++ interface from version 2022.11.001
inline void elpa_eigenvectors(const elpa_t handle, double *a, double *ev, double *q, int *error)
{
    elpa_eigenvectors_a_h_a_d(handle, a, ev, q, error);
}

inline void elpa_eigenvectors(const elpa_t handle, float *a, float *ev, float *q, int *error)
{
    elpa_eigenvectors_a_h_a_f(handle, a, ev, q, error);
}

inline void elpa_eigenvectors(const elpa_t handle, std::complex<double> *a, double *ev, std::complex<double> *q, int *error)
{
    elpa_eigenvectors_a_h_a_dc(handle, reinterpret_cast<double _Complex*>(a),
                               ev, reinterpret_cast<double _Complex*>(q), error);
}

inline void elpa_eigenvectors(const elpa_t handle, std::complex<float>  *a, float  *ev, std::complex<float>  *q, int *error)
{
    elpa_eigenvectors_a_h_a_fc(handle, reinterpret_cast<float _Complex*>(a),
                               ev, reinterpret_cast<float _Complex*>(q), error);
}
#endif

/*! \brief generic C method for elpa_skew_eigenvectors
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a
 *  \param  ev      on return: float/double pointer to eigenvalues
 *  \param  q       on return: float/double float complex/double complex pointer to eigenvectors
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
#if ELPA_API_VERSION <= 20210502 // ELPA 2021.05.002 and earlier versions
inline void elpa_skew_eigenvectors(const elpa_t handle, double *a, double *ev, double *q, int *error)
{
    elpa_eigenvectors_d(handle, a, ev, q, error);
}

inline void elpa_skew_eigenvectors(const elpa_t handle, float  *a, float  *ev, float  *q, int *error)
{
    elpa_eigenvectors_f(handle, a, ev, q, error);
}
#elif ELPA_API_VERSION < 20220501   // ELPA version between 2021.11.001 and 2022.05.001
inline void elpa_skew_eigenvectors(const elpa_t handle, double *a, double *ev, double *q, int *error)
{
    elpa_eigenvectors_all_host_arrays_d(handle, a, ev, q, error);
}

inline void elpa_skew_eigenvectors(const elpa_t handle, float  *a, float  *ev, float  *q, int *error)
{
    elpa_eigenvectors_all_host_arrays_f(handle, a, ev, q, error);
}
#else // ELPA version 2022.05.001, ELPA has its own c++ interface from version 2022.11.001
inline void elpa_skew_eigenvectors(const elpa_t handle, double *a, double *ev, double *q, int *error)
{
    elpa_skew_eigenvectors_a_h_a_d(handle, a, ev, q, error);
}

inline void elpa_skew_eigenvectors(const elpa_t handle, float  *a, float  *ev, float  *q, int *error)
{
    elpa_skew_eigenvectors_a_h_a_f(handle, a, ev, q, error);
}
#endif



/*! \brief generic C method for elpa_generalized_eigenvectors
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a
 *  \param  b       float/double float complex/double complex pointer to matrix b
 *  \param  ev      on return: float/double pointer to eigenvalues
 *  \param  q       on return: float/double float complex/double complex pointer to eigenvectors
 *  \param  is_already_decomposed   set to 1, if b already decomposed by previous call to elpa_generalized
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
inline void elpa_generalized_eigenvectors(elpa_t handle, double *a, double *b, double *ev, double *q, int is_already_decomposed, int *error)
{
    elpa_generalized_eigenvectors_d(handle, a, b, ev, q, is_already_decomposed, error);
}

inline void elpa_generalized_eigenvectors(elpa_t handle, float  *a, float  *b, float  *ev, float  *q, int is_already_decomposed, int *error)
{
    elpa_generalized_eigenvectors_f(handle, a, b, ev, q, is_already_decomposed, error);
}

inline void elpa_generalized_eigenvectors(elpa_t handle, std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, int is_already_decomposed, int *error)
{
    elpa_generalized_eigenvectors_dc(handle, reinterpret_cast<double _Complex*>(a), reinterpret_cast<double _Complex*>(b),
                                     ev, reinterpret_cast<double _Complex*>(q), is_already_decomposed, error);
}

inline void elpa_generalized_eigenvectors(elpa_t handle, std::complex<float>  *a, std::complex<float>  *b, float  *ev, std::complex<float>  *q, int is_already_decomposed, int *error)
{
    elpa_generalized_eigenvectors_fc(handle, reinterpret_cast<float _Complex*>(a), reinterpret_cast<float _Complex*>(b),
                                     ev, reinterpret_cast<float _Complex*>(q), is_already_decomposed, error);
}

/*! \brief generic C method for elpa_eigenvalues
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a
 *  \param  ev      on return: float/double pointer to eigenvalues
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
#if ELPA_API_VERSION <= 20210502 // ELPA 2021.05.002 and earlier versions
inline void elpa_eigenvalues(elpa_t handle, double *a, double *ev, int *error)
{
    elpa_eigenvalues_d(handle, a, ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, float  *a, float  *ev, int *error)
{
    elpa_eigenvalues_f(handle, a, ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, std::complex<double> *a, double *ev, int *error)
{
    elpa_eigenvalues_dc(handle, reinterpret_cast<double _Complex*>(a), ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, std::complex<float>  *a, float  *ev, int *error)
{
    elpa_eigenvalues_fc (handle, reinterpret_cast<float _Complex*>(a), ev, error);
}
#elif ELPA_API_VERSION < 20220501   // ELPA version between 2021.11.001 and 2022.05.001
inline void elpa_eigenvalues(elpa_t handle, double *a, double *ev, int *error)
{
    elpa_eigenvalues_all_host_arrays_d(handle, a, ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, float  *a, float  *ev, int *error)
{
    elpa_eigenvalues_all_host_arrays_f(handle, a, ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, std::complex<double> *a, double *ev, int *error)
{
    elpa_eigenvalues_all_host_arrays_dc(handle, reinterpret_cast<double _Complex*>(a), ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, std::complex<float>  *a, float  *ev, int *error)
{
    elpa_eigenvalues_all_host_arrays_fc(handle, reinterpret_cast<float _Complex*>(a), ev, error);
}
#else // ELPA version 2022.05.001, ELPA has its own c++ interface from version 2022.11.001
inline void elpa_eigenvalues(elpa_t handle, double *a, double *ev, int *error)
{
    elpa_eigenvalues_a_h_a_d(handle, a, ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, float  *a, float  *ev, int *error)
{
    elpa_eigenvalues_a_h_a_f(handle, a, ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, std::complex<double> *a, double *ev, int *error)
{
    elpa_eigenvalues_a_h_a_dc(handle, reinterpret_cast<double _Complex*>(a), ev, error);
}
inline void elpa_eigenvalues(elpa_t handle, std::complex<float>  *a, float  *ev, int *error)
{
    elpa_eigenvalues_a_h_a_fc(handle, reinterpret_cast<float _Complex*>(a), ev, error);
}
#endif

/*! \brief generic C method for elpa_skew_eigenvalues
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a
 *  \param  ev      on return: float/double pointer to eigenvalues
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
#if ELPA_API_VERSION <= 20210502 // ELPA 2021.05.002 and earlier versions
inline void elpa_skew_eigenvalues(elpa_t handle, double *a, double *ev, int *error)
{
    elpa_eigenvalues_d(handle, a, ev, error);
}
inline void elpa_skew_eigenvalues(elpa_t handle, float  *a, float  *ev, int *error)
{
    elpa_eigenvalues_f(handle, a, ev, error);
}
#elif ELPA_API_VERSION < 20220501   // ELPA version between 2021.11.001 and 2022.05.001
inline void elpa_skew_eigenvalues(elpa_t handle, double *a, double *ev, int *error)
{
    elpa_eigenvalues_all_host_arrays_d(handle, a, ev, error);
}
inline void elpa_skew_eigenvalues(elpa_t handle, float  *a, float  *ev, int *error)
{
    elpa_eigenvalues_all_host_arrays_f(handle, a, ev, error);
}
#else // ELPA version 2022.05.001, ELPA has its own c++ interface from version 2022.11.001
inline void elpa_skew_eigenvalues(elpa_t handle, double *a, double *ev, int *error)
{
    elpa_eigenvalues_a_h_a_d(handle, a, ev, error);
}
inline void elpa_skew_eigenvalues(elpa_t handle, float  *a, float  *ev, int *error)
{
    elpa_eigenvalues_a_h_a_f(handle, a, ev, error);
}
#endif

/*! \brief generic C method for elpa_cholesky
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a, for which
 *                  the cholesky factorizaion will be computed
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */

#if ELPA_API_VERSION < 20220501   // ELPA version before 2022.05.001
inline void elpa_cholesky(elpa_t handle, double *a, int *error)
{
    elpa_cholesky_d(handle, a, error);
}
inline void elpa_cholesky(elpa_t handle, float  *a, int *error)
{
    elpa_cholesky_f(handle, a, error);
}
inline void elpa_cholesky(elpa_t handle, std::complex<double> *a, int *error)
{
    elpa_cholesky_dc(handle, reinterpret_cast<double _Complex*>(a), error);
}
inline void elpa_cholesky(elpa_t handle, std::complex<float>  *a, int *error)
{
    elpa_cholesky_fc(handle, reinterpret_cast<float _Complex*>(a), error);
}
#else
inline void elpa_cholesky(elpa_t handle, double *a, int *error)
{
    elpa_cholesky_a_h_a_d(handle, a, error);
}
inline void elpa_cholesky(elpa_t handle, float  *a, int *error)
{
    elpa_cholesky_a_h_a_f(handle, a, error);
}
inline void elpa_cholesky(elpa_t handle, std::complex<double> *a, int *error)
{
    elpa_cholesky_a_h_a_dc(handle, reinterpret_cast<double _Complex*>(a), error);
}
inline void elpa_cholesky(elpa_t handle, std::complex<float>  *a, int *error)
{
    elpa_cholesky_a_h_a_fc(handle, reinterpret_cast<float _Complex*>(a), error);
}
#endif

/*! \brief generic C method for elpa_hermitian_multiply
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  uplo_a  descriptor for matrix a
 *  \param  uplo_c  descriptor for matrix c
 *  \param  ncb     int
 *  \param  a       float/double float complex/double complex pointer to matrix a
 *  \param  b       float/double float complex/double complex pointer to matrix b
 *  \param  nrows_b number of rows for matrix b
 *  \param  ncols_b number of cols for matrix b
 *  \param  c       float/double float complex/double complex pointer to matrix c
 *  \param  nrows_c number of rows for matrix c
 *  \param  ncols_c number of cols for matrix c
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
#if ELPA_API_VERSION < 20220501   // ELPA version before 2022.05.001
inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, double *a, double *b, int nrows_b, int ncols_b, double *c, int nrows_c, int ncols_c, int *error)
{
    elpa_hermitian_multiply_d(handle, uplo_a, uplo_c, ncb, a, b, nrows_b, ncols_b, c, nrows_c, ncols_c, error);
}
inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, float  *a, float  *b, int nrows_b, int ncols_b, float  *c, int nrows_c, int ncols_c, int *error)
{
    elpa_hermitian_multiply_df(handle, uplo_a, uplo_c, ncb, a, b, nrows_b, ncols_b, c, nrows_c, ncols_c, error);
}
inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, std::complex<double> *a, std::complex<double> *b, int nrows_b, int ncols_b, std::complex<double> *c, int nrows_c, int ncols_c, int *error)
{
    elpa_hermitian_multiply_dc(handle, uplo_a, uplo_c, ncb, reinterpret_cast<double _Complex*>(a),
                               reinterpret_cast<double _Complex*>(b), nrows_b, ncols_b,
                               reinterpret_cast<double _Complex*>(c), nrows_c, ncols_c, error);
}
inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, std::complex<float>  *a, std::complex<float>  *b, int nrows_b, int ncols_b, std::complex<float>  *c, int nrows_c, int ncols_c, int *error)
{
    elpa_hermitian_multiply_fc(handle, uplo_a, uplo_c, ncb, reinterpret_cast<float _Complex*>(a),
                               reinterpret_cast<float _Complex*>(b), nrows_b, ncols_b,
                               reinterpret_cast<float _Complex*>(c), nrows_c, ncols_c, error);
}
#else
inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, double *a, double *b, int nrows_b, int ncols_b, double *c, int nrows_c, int ncols_c, int *error)
{
    elpa_hermitian_multiply_a_h_a_d(handle, uplo_a, uplo_c, ncb, a, b, nrows_b, ncols_b, c, nrows_c, ncols_c, error);
}
inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, float  *a, float  *b, int nrows_b, int ncols_b, float  *c, int nrows_c, int ncols_c, int *error)
{
    elpa_hermitian_multiply_a_h_a_f(handle, uplo_a, uplo_c, ncb, a, b, nrows_b, ncols_b, c, nrows_c, ncols_c, error);
}
inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, std::complex<double> *a, std::complex<double> *b, int nrows_b, int ncols_b, std::complex<double> *c, int nrows_c, int ncols_c, int *error)
{
    elpa_hermitian_multiply_a_h_a_dc(handle, uplo_a, uplo_c, ncb, reinterpret_cast<double _Complex*>(a),
                               reinterpret_cast<double _Complex*>(b), nrows_b, ncols_b,
                               reinterpret_cast<double _Complex*>(c), nrows_c, ncols_c, error);
}
inline void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, std::complex<float>  *a, std::complex<float>  *b, int nrows_b, int ncols_b, std::complex<float>  *c, int nrows_c, int ncols_c, int *error)
{
    elpa_hermitian_multiply_a_h_a_fc(handle, uplo_a, uplo_c, ncb, reinterpret_cast<float _Complex*>(a),
                               reinterpret_cast<float _Complex*>(b), nrows_b, ncols_b,
                               reinterpret_cast<float _Complex*>(c), nrows_c, ncols_c, error);
}
#endif

/*! \brief generic C method for elpa_invert_triangular
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a, which
 *                  should be inverted
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
#if ELPA_API_VERSION < 20220501   // ELPA version before 2022.05.001
inline void elpa_invert_triangular(elpa_t handle, double *a, int *error)
{
    elpa_invert_trm_d(handle, a, error);
}
inline void elpa_invert_triangular(elpa_t handle, float  *a, int *error)
{
    elpa_invert_trm_f(handle, a, error);
}
inline void elpa_invert_triangular(elpa_t handle, std::complex<double> *a, int *error)
{
    elpa_invert_trm_dc(handle, reinterpret_cast<double _Complex*>(a), error);
}
inline void elpa_invert_triangular(elpa_t handle, std::complex<float>  *a, int *error)
{
    elpa_invert_trm_fc(handle, reinterpret_cast<float _Complex*>(a), error);
}
#else
inline void elpa_invert_triangular(elpa_t handle, double *a, int *error)
{
    elpa_invert_trm_a_h_a_d(handle, a, error);
}
inline void elpa_invert_triangular(elpa_t handle, float  *a, int *error)
{
    elpa_invert_trm_a_h_a_f(handle, a, error);
}
inline void elpa_invert_triangular(elpa_t handle, std::complex<double> *a, int *error)
{
    elpa_invert_trm_a_h_a_dc(handle, reinterpret_cast<double _Complex*>(a), error);
}
inline void elpa_invert_triangular(elpa_t handle, std::complex<float>  *a, int *error)
{
    elpa_invert_trm_a_h_a_fc(handle, reinterpret_cast<float _Complex*>(a), error);
}
#endif
