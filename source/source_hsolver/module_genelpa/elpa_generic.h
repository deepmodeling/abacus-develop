#ifndef ABACUS_SOURCE_HSOLVER_MODULE_GENELPA_ELPA_GENERIC_H
#define ABACUS_SOURCE_HSOLVER_MODULE_GENELPA_ELPA_GENERIC_H

#include "elpa_new.h"

#include <complex>

// C++ overloads that mirror the C11 _Generic dispatch in `elpa/elpa_generic.h`.
// Only used with ELPA versions older than 2022.11.001, which ship no C++ interface;
// see elpa_new.h.

/*! \brief generic C method for elpa_set
 *
 *  \details
 *  \param  handle  handle of the ELPA object for which a key/value pair should be set
 *  \param  name    the name of the key
 *  \param  value   integer/double value to be set for the key
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
void elpa_set(elpa_t handle, const char *name, int value, int *error);
void elpa_set(elpa_t handle, const char *name, double value, int *error);

/*! \brief generic C method for elpa_get
 *
 *  \details
 *  \param  handle  handle of the ELPA object for which a key/value pair should be queried
 *  \param  name    the name of the key
 *  \param  value   integer/double value to be queried
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
void elpa_get(elpa_t handle, const char *name, int *value, int *error);
void elpa_get(elpa_t handle, const char *name, double *value, int *error);

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
void elpa_eigenvectors(const elpa_t handle, double *a, double *ev, double *q, int *error);
void elpa_eigenvectors(const elpa_t handle, float  *a, float  *ev, float  *q, int *error);
void elpa_eigenvectors(const elpa_t handle, std::complex<double> *a, double *ev, std::complex<double> *q, int *error);
void elpa_eigenvectors(const elpa_t handle, std::complex<float>  *a, float  *ev, std::complex<float>  *q, int *error);

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
void elpa_skew_eigenvectors(const elpa_t handle, double *a, double *ev, double *q, int *error);
void elpa_skew_eigenvectors(const elpa_t handle, float  *a, float  *ev, float  *q, int *error);

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
void elpa_generalized_eigenvectors(elpa_t handle, double *a, double *b, double *ev, double *q, int is_already_decomposed, int *error);
void elpa_generalized_eigenvectors(elpa_t handle, float  *a, float  *b, float  *ev, float  *q, int is_already_decomposed, int *error);
void elpa_generalized_eigenvectors(elpa_t handle, std::complex<double> *a, std::complex<double> *b, double *ev, std::complex<double> *q, int is_already_decomposed, int *error);
void elpa_generalized_eigenvectors(elpa_t handle, std::complex<float>  *a, std::complex<float>  *b, float  *ev, std::complex<float>  *q, int is_already_decomposed, int *error);

/*! \brief generic C method for elpa_eigenvalues
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a
 *  \param  ev      on return: float/double pointer to eigenvalues
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
void elpa_eigenvalues(elpa_t handle, double *a, double *ev, int *error);
void elpa_eigenvalues(elpa_t handle, float  *a, float  *ev, int *error);
void elpa_eigenvalues(elpa_t handle, std::complex<double> *a, double *ev, int *error);
void elpa_eigenvalues(elpa_t handle, std::complex<float>  *a, float  *ev, int *error);

/*! \brief generic C method for elpa_skew_eigenvalues
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a
 *  \param  ev      on return: float/double pointer to eigenvalues
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
void elpa_skew_eigenvalues(elpa_t handle, double *a, double *ev, int *error);
void elpa_skew_eigenvalues(elpa_t handle, float  *a, float  *ev, int *error);

/*! \brief generic C method for elpa_cholesky
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a, for which
 *                  the cholesky factorizaion will be computed
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */

void elpa_cholesky(elpa_t handle, double *a, int *error);
void elpa_cholesky(elpa_t handle, float  *a, int *error);
void elpa_cholesky(elpa_t handle, std::complex<double> *a, int *error);
void elpa_cholesky(elpa_t handle, std::complex<float>  *a, int *error);

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
void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, double *a, double *b, int nrows_b, int ncols_b, double *c, int nrows_c, int ncols_c, int *error);
void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, float  *a, float  *b, int nrows_b, int ncols_b, float  *c, int nrows_c, int ncols_c, int *error);
void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, std::complex<double> *a, std::complex<double> *b, int nrows_b, int ncols_b, std::complex<double> *c, int nrows_c, int ncols_c, int *error);
void elpa_hermitian_multiply(elpa_t handle, char uplo_a, char uplo_c, int ncb, std::complex<float>  *a, std::complex<float>  *b, int nrows_b, int ncols_b, std::complex<float>  *c, int nrows_c, int ncols_c, int *error);

/*! \brief generic C method for elpa_invert_triangular
 *
 *  \details
 *  \param  handle  handle of the ELPA object, which defines the problem
 *  \param  a       float/double float complex/double complex pointer to matrix a, which
 *                  should be inverted
 *  \param  error   on return the error code, which can be queried with elpa_strerr()
 *  \result void
 */
void elpa_invert_triangular(elpa_t handle, double *a, int *error);
void elpa_invert_triangular(elpa_t handle, float  *a, int *error);
void elpa_invert_triangular(elpa_t handle, std::complex<double> *a, int *error);
void elpa_invert_triangular(elpa_t handle, std::complex<float>  *a, int *error);

#endif // ABACUS_SOURCE_HSOLVER_MODULE_GENELPA_ELPA_GENERIC_H
