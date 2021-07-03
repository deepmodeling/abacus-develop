 #include<complex>
 /*! \brief C old, deprecated interface, will be deleted. Use "elpa_get_communicators"
 *
 * \param mpi_comm_word    MPI global communicator (in)
 * \param my_prow          Row coordinate of the calling process in the process grid (in)
 * \param my_pcol          Column coordinate of the calling process in the process grid (in)
 * \param mpi_comm_rows    Communicator for communicating within rows of processes (out)
 * \result int             integer error value of mpi_comm_split function
 */
 int get_elpa_row_col_comms(MPI_Comm mpi_comm_world, int my_prow, int my_pcol, int *mpi_comm_rows, int *mpi_comm_cols);
 /*! \brief C old, deprecated interface, will be deleted. Use "elpa_get_communicators"
 *
 * \param mpi_comm_word    MPI global communicator (in)
 * \param my_prow          Row coordinate of the calling process in the process grid (in)
 * \param my_pcol          Column coordinate of the calling process in the process grid (in)
 * \param mpi_comm_rows    Communicator for communicating within rows of processes (out)
 * \result int             integer error value of mpi_comm_split function
 */
 int get_elpa_communicators(MPI_Comm mpi_comm_world, int my_prow, int my_pcol, int *mpi_comm_rows, int *mpi_comm_cols);
 /*! \brief C interface to create ELPA communicators
 *
 * \param mpi_comm_word    MPI global communicator (in)
 * \param my_prow          Row coordinate of the calling process in the process grid (in)
 * \param my_pcol          Column coordinate of the calling process in the process grid (in)
 * \param mpi_comm_rows    Communicator for communicating within rows of processes (out)
 * \result int             integer error value of mpi_comm_split function
 */
 int elpa_get_communicators(MPI_Comm mpi_comm_world, int my_prow, int my_pcol, int *mpi_comm_rows, int *mpi_comm_cols);
  /*! \brief C interface to solve the real eigenvalue problem with 1-stage solver
  *
 *  \param  na                   Order of matrix a
 *  \param  nev                  Number of eigenvalues needed.
 *                               The smallest nev eigenvalues/eigenvectors are calculated.
 *  \param  a                    Distributed matrix for which eigenvalues are to be computed.
 *                               Distribution is like in Scalapack.
 *                               The full matrix must be set (not only one half like in scalapack).
 *  \param lda                   Leading dimension of a
 *  \param ev(na)                On output: eigenvalues of a, every processor gets the complete set
 *  \param q                     On output: Eigenvectors of a
 *                               Distribution is like in Scalapack.
 *                               Must be always dimensioned to the full size (corresponding to (na,na))
 *                               even if only a part of the eigenvalues is needed.
 *  \param ldq                   Leading dimension of q
 *  \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
 *  \param matrixCols           distributed number of matrix columns
 *  \param mpi_comm_rows        MPI-Communicator for rows
 *  \param mpi_comm_cols        MPI-Communicator for columns
 *
 *  \result                     int: 1 if error occured, otherwise 0
*/
 int elpa_solve_evp_real_1stage(int na, int nev, double *a, int lda, double *ev, double *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols);
 /*! \brief C interface to solve the complex eigenvalue problem with 1-stage solver
 *
 *  \param  na                   Order of matrix a
 *  \param  nev                  Number of eigenvalues needed.
 *                               The smallest nev eigenvalues/eigenvectors are calculated.
 *  \param  a                    Distributed matrix for which eigenvalues are to be computed.
 *                               Distribution is like in Scalapack.
 *                               The full matrix must be set (not only one half like in scalapack).
 *  \param lda                   Leading dimension of a
 *  \param ev(na)                On output: eigenvalues of a, every processor gets the complete set
 *  \param q                     On output: Eigenvectors of a
 *                               Distribution is like in Scalapack.
 *                               Must be always dimensioned to the full size (corresponding to (na,na))
 *                               even if only a part of the eigenvalues is needed.
 *  \param ldq                   Leading dimension of q
 *  \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
 *  \param matrixCols           distributed number of matrix columns
 *  \param mpi_comm_rows        MPI-Communicator for rows
 *  \param mpi_comm_cols        MPI-Communicator for columns
 *
 *  \result                     int: 1 if error occured, otherwise 0
 */
 int elpa_solve_evp_complex_1stage(int na, int nev, std::complex<double> *a, int lda, double *ev, std::complex<double> *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols);
 /*! \brief C interface to solve the real eigenvalue problem with 2-stage solver
 *
 *  \param  na                        Order of matrix a
 *  \param  nev                       Number of eigenvalues needed.
 *                                    The smallest nev eigenvalues/eigenvectors are calculated.
 *  \param  a                         Distributed matrix for which eigenvalues are to be computed.
 *                                    Distribution is like in Scalapack.
 *                                    The full matrix must be set (not only one half like in scalapack).
 *  \param lda                        Leading dimension of a
 *  \param ev(na)                     On output: eigenvalues of a, every processor gets the complete set
 *  \param q                          On output: Eigenvectors of a
 *                                    Distribution is like in Scalapack.
 *                                    Must be always dimensioned to the full size (corresponding to (na,na))
 *                                    even if only a part of the eigenvalues is needed.
 *  \param ldq                        Leading dimension of q
 *  \param nblk                       blocksize of cyclic distribution, must be the same in both directions!
 *  \param matrixCols                 distributed number of matrix columns
 *  \param mpi_comm_rows              MPI-Communicator for rows
 *  \param mpi_comm_cols              MPI-Communicator for columns
 *  \param mpi_coll_all               MPI communicator for the total processor set
 *  \param THIS_REAL_ELPA_KERNEL_API  specify used ELPA2 kernel via API
 *  \param use_qr                     use QR decomposition 1 = yes, 0 = no
 *
 *  \result                     int: 1 if error occured, otherwise 0
 */
 int elpa_solve_evp_real_2stage(int na, int nev, double *a, int lda, double *ev, double *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, MPI_Comm mpi_comm_all, int THIS_REAL_ELPA_KERNEL_API, int useQR);
 /*! \brief C interface to solve the complex eigenvalue problem with 2-stage solver
 *
 *  \param  na                        Order of matrix a
 *  \param  nev                       Number of eigenvalues needed.
 *                                    The smallest nev eigenvalues/eigenvectors are calculated.
 *  \param  a                         Distributed matrix for which eigenvalues are to be computed.
 *                                    Distribution is like in Scalapack.
 *                                    The full matrix must be set (not only one half like in scalapack).
 *  \param lda                        Leading dimension of a
 *  \param ev(na)                     On output: eigenvalues of a, every processor gets the complete set
 *  \param q                          On output: Eigenvectors of a
 *                                    Distribution is like in Scalapack.
 *                                    Must be always dimensioned to the full size (corresponding to (na,na))
 *                                    even if only a part of the eigenvalues is needed.
 *  \param ldq                        Leading dimension of q
 *  \param nblk                       blocksize of cyclic distribution, must be the same in both directions!
 *  \param matrixCols                 distributed number of matrix columns
 *  \param mpi_comm_rows              MPI-Communicator for rows
 *  \param mpi_comm_cols              MPI-Communicator for columns
 *  \param mpi_coll_all               MPI communicator for the total processor set
 *  \param THIS_COMPLEX_ELPA_KERNEL_API  specify used ELPA2 kernel via API
 *
 *  \result                     int: 1 if error occured, otherwise 0
 */
 int elpa_solve_evp_complex_2stage(int na, int nev, std::complex<double> *a, int lda, double *ev, std::complex<double> *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, MPI_Comm mpi_comm_all, int THIS_COMPLEX_ELPA_KERNEL_API);
 /*! \brief C interface to driver function "elpa_solve_evp_real"
 *
 *  \param  na                        Order of matrix a
 *  \param  nev                       Number of eigenvalues needed.
 *                                    The smallest nev eigenvalues/eigenvectors are calculated.
 *  \param  a                         Distributed matrix for which eigenvalues are to be computed.
 *                                    Distribution is like in Scalapack.
 *                                    The full matrix must be set (not only one half like in scalapack).
 *  \param lda                        Leading dimension of a
 *  \param ev(na)                     On output: eigenvalues of a, every processor gets the complete set
 *  \param q                          On output: Eigenvectors of a
 *                                    Distribution is like in Scalapack.
 *                                    Must be always dimensioned to the full size (corresponding to (na,na))
 *                                    even if only a part of the eigenvalues is needed.
 *  \param ldq                        Leading dimension of q
 *  \param nblk                       blocksize of cyclic distribution, must be the same in both directions!
 *  \param matrixCols                 distributed number of matrix columns
 *  \param mpi_comm_rows              MPI-Communicator for rows
 *  \param mpi_comm_cols              MPI-Communicator for columns
 *  \param mpi_coll_all               MPI communicator for the total processor set
 *  \param THIS_REAL_ELPA_KERNEL_API  specify used ELPA2 kernel via API
 *  \param use_qr                     use QR decomposition 1 = yes, 0 = no
 *  \param method                      choose whether to use ELPA 1stage or 2stage solver
 *                                     possible values: "1stage" => use ELPA 1stage solver
 *                                                      "2stage" => use ELPA 2stage solver
 *                                                       "auto"   => (at the moment) use ELPA 2stage solver
 *
 *  \result                     int: 1 if error occured, otherwise 0
 */
 int elpa_solve_evp_real(int na, int nev, double *a, int lda, double *ev, double *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, MPI_Comm mpi_comm_all, int THIS_REAL_ELPA_KERNEL_API, int useQR, char *method);
 /*! \brief C interface to driver function "elpa_solve_evp_complex"
 *
 *  \param  na                           Order of matrix a
 *  \param  nev                          Number of eigenvalues needed.
 *                                       The smallest nev eigenvalues/eigenvectors are calculated.
 *  \param  a                            Distributed matrix for which eigenvalues are to be computed.
 *                                       Distribution is like in Scalapack.
 *                                       The full matrix must be set (not only one half like in scalapack).
 *  \param lda                           Leading dimension of a
 *  \param ev(na)                        On output: eigenvalues of a, every processor gets the complete set
 *  \param q                             On output: Eigenvectors of a
 *                                       Distribution is like in Scalapack.
 *                                       Must be always dimensioned to the full size (corresponding to (na,na))
 *                                       even if only a part of the eigenvalues is needed.
 *  \param ldq                           Leading dimension of q
 *  \param nblk                          blocksize of cyclic distribution, must be the same in both directions!
 *  \param matrixCols                    distributed number of matrix columns
 *  \param mpi_comm_rows                 MPI-Communicator for rows
 *  \param mpi_comm_cols                 MPI-Communicator for columns
 *  \param mpi_coll_all                  MPI communicator for the total processor set
 *  \param THIS_COMPLEX_ELPA_KERNEL_API  specify used ELPA2 kernel via API
 *  \param method                        choose whether to use ELPA 1stage or 2stage solver
 *                                       possible values: "1stage" => use ELPA 1stage solver
 *                                                        "2stage" => use ELPA 2stage solver
 *                                                         "auto"   => (at the moment) use ELPA 2stage solver
 *
 *  \result                     int: 1 if error occured, otherwise 0
 */
 int elpa_solve_evp_complex(int na, int nev, std::complex<double> *a, int lda, double *ev, std::complex<double> *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, MPI_Comm mpi_comm_all, int THIS_COMPLEX_ELPA_KERNEL_API, char *method);
 /*! \brief  C interface to solve tridiagonal eigensystem with divide and conquer method
 *\details

 *\param na                    Matrix dimension
 *\param nev                   number of eigenvalues/vectors to be computed
 *\param d                     array d(na) on input diagonal elements of tridiagonal matrix, on
 *                             output the eigenvalues in ascending order
 *\param e                     array e(na) on input subdiagonal elements of matrix, on exit destroyed
 *\param q                     on exit : matrix q(ldq,matrixCols) contains the eigenvectors
 *\param ldq                   leading dimension of matrix q
 *\param nblk                  blocksize of cyclic distribution, must be the same in both directions!
 *\param matrixCols            columns of matrix q
 *\param mpi_comm_rows         MPI communicator for rows
 *\param mpi_comm_cols         MPI communicator for columns
 *\param wantDebug             give more debug information if 1, else 0
 *\result success              int 1 on success, else 0
 */
 int elpa_solve_tridi(int na, int nev, double *d, double *e, double *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
 /*! \brief  C interface for elpa_mult_at_b_real: Performs C : = A**T * B
 *        where   A is a square matrix (na,na) which is optionally upper or lower triangular
 *                B is a (na,ncb) matrix
 *                C is a (na,ncb) matrix where optionally only the upper or lower
 *                  triangle may be computed
 *\details
 *\param  uplo_a               'U' if A is upper triangular
 *                             'L' if A is lower triangular
 *                             anything else if A is a full matrix
 *                             Please note: This pertains to the original A (as set in the calling program)
 *                                          whereas the transpose of A is used for calculations
 *                             If uplo_a is 'U' or 'L', the other triangle is not used at all,
 *                             i.e. it may contain arbitrary numbers
 *\param uplo_c                'U' if only the upper diagonal part of C is needed
 *                             'L' if only the upper diagonal part of C is needed
 *                             anything else if the full matrix C is needed
 *                             Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
 *                                           written to a certain extent, i.e. one shouldn't rely on the content there!
 *\param na                    Number of rows/columns of A, number of rows of B and C
 *\param ncb                   Number of columns  of B and C
 *\param a                     matrix a
 *\param lda                   leading dimension of matrix a
 *\param ldaCols               columns of matrix a
 *\param b                     matrix b
 *\param ldb                   leading dimension of matrix b
 *\param ldbCols               columns of matrix b
 *\param nblk                  blocksize of cyclic distribution, must be the same in both directions!
 *\param  mpi_comm_rows        MPI communicator for rows
 *\param  mpi_comm_cols        MPI communicator for columns
 *\param c                     matrix c
 *\param ldc                   leading dimension of matrix c
 *\param ldcCols               columns of matrix c
 *\result success              int report success (1) or failure (0)
 */
 int elpa_mult_at_b_real(char uplo_a, char uplo_c, int na, int ncb, double *a, int lda, int ldaCols, double *b, int ldb, int ldbCols, int nlbk, int mpi_comm_rows, int mpi_comm_cols, double *c, int ldc, int ldcCols);
 /*! \brief C interface for elpa_mult_ah_b_complex: Performs C : = A**H * B
 *        where   A is a square matrix (na,na) which is optionally upper or lower triangular
 *                B is a (na,ncb) matrix
 *                C is a (na,ncb) matrix where optionally only the upper or lower
 *                  triangle may be computed
 *\details

 *\param  uplo_a               'U' if A is upper triangular
 *                             'L' if A is lower triangular
 *                             anything else if A is a full matrix
 *                             Please note: This pertains to the original A (as set in the calling program)
 *                                          whereas the transpose of A is used for calculations
 *                             If uplo_a is 'U' or 'L', the other triangle is not used at all,
 *                             i.e. it may contain arbitrary numbers
 *\param uplo_c                'U' if only the upper diagonal part of C is needed
 *                             'L' if only the upper diagonal part of C is needed
 *                             anything else if the full matrix C is needed
 *                             Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
 *                                           written to a certain extent, i.e. one shouldn't rely on the content there!
 *\param na                    Number of rows/columns of A, number of rows of B and C
 *\param ncb                   Number of columns  of B and C
 *\param a                     matrix a
 *\param lda                   leading dimension of matrix a
 *\param ldaCols               columns of matrix a
 *\param b                     matrix b
 *\param ldb                   leading dimension of matrix b
 *\param ldbCols               columns of matrix b
 *\param nblk                  blocksize of cyclic distribution, must be the same in both directions!
 *\param  mpi_comm_rows        MPI communicator for rows
 *\param  mpi_comm_cols        MPI communicator for columns
 *\param c                     matrix c
 *\param ldc                   leading dimension of matrix c
 *\param ldcCols               columns of matrix c
 *\result success              int reports success (1) or failure (0)
 */
 int elpa_mult_ah_b_complex(char uplo_a, char uplo_c, int na, int ncb, std::complex<double> *a, int lda, std::complex<double> *b, int ldb, int nblk, int mpi_comm_rows, int mpi_comm_cols, std::complex<double> *c, int ldc);
 /*! \brief  C interface to elpa_invert_trm_real: Inverts a upper triangular matrix
 *\details
 *\param  na                   Order of matrix
 *\param  a(lda,matrixCols)    Distributed matrix which should be inverted
 *                             Distribution is like in Scalapack.
 *                             Only upper triangle is needs to be set.
 *                             The lower triangle is not referenced.
 *\param  lda                  Leading dimension of a
 *\param                       matrixCols  local columns of matrix a
 *\param  nblk                 blocksize of cyclic distribution, must be the same in both directions!
 *\param  mpi_comm_rows        MPI communicator for rows
 *\param  mpi_comm_cols        MPI communicator for columns
 *\param wantDebug             int more debug information on failure if 1, else 0
 *\result succes               int reports success (1) or failure (0)
 */
 int elpa_invert_trm_real(int na, double *a, int lda, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
 /*! \brief  C interface to elpa_invert_trm_complex: Inverts a complex upper triangular matrix
 *\details
 *\param  na                   Order of matrix
 *\param  a(lda,matrixCols)    Distributed matrix which should be inverted
 *                             Distribution is like in Scalapack.
 *                             Only upper triangle is needs to be set.
 *                             The lower triangle is not referenced.
 *\param  lda                  Leading dimension of a
 *\param                       matrixCols  local columns of matrix a
 *\param  nblk                 blocksize of cyclic distribution, must be the same in both directions!
 *\param  mpi_comm_rows        MPI communicator for rows
 *\param  mpi_comm_cols        MPI communicator for columns
 *\param wantDebug             int more debug information on failure if 1, else 0
 *\result succes               int reports success (1) or failure (0)
 */
 int elpa_invert_trm_complex(int na, std::complex<double> *a, int lda, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
 /*! \brief  elpa_cholesky_real: Cholesky factorization of a real symmetric matrix
 *\details

 *\param  na                   Order of matrix
 *\param  a(lda,matrixCols)    Distributed matrix which should be factorized.
 *                             Distribution is like in Scalapack.
 *                             Only upper triangle is needs to be set.
 *                             On return, the upper triangle contains the Cholesky factor
 *                             and the lower triangle is set to 0.
 *\param  lda                  Leading dimension of a
 *\param  matrixCols           local columns of matrix a
 *\param  nblk                 blocksize of cyclic distribution, must be the same in both directions!
 *\param  mpi_comm_rows        MPI communicator for rows
 *\param  mpi_comm_cols        MPI communicator for columns
 *\param wantDebug             int more debug information on failure if 1, else 0
 *\result succes               int reports success (1) or failure (0)
 */
 int elpa_cholesky_real(int na, double *a, int lda, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
 /*! \brief  C interface elpa_cholesky_complex: Cholesky factorization of a complex hermitian matrix
 *\details
 *\param  na                   Order of matrix
 *\param  a(lda,matrixCols)    Distributed matrix which should be factorized.
 *                             Distribution is like in Scalapack.
 *                             Only upper triangle is needs to be set.
 *                             On return, the upper triangle contains the Cholesky factor
 *                             and the lower triangle is set to 0.
 *\param  lda                  Leading dimension of a
 *\param                       matrixCols  local columns of matrix a
 *\param  nblk                 blocksize of cyclic distribution, must be the same in both directions!
 *\param  mpi_comm_rows        MPI communicator for rows
 *\param  mpi_comm_cols        MPI communicator for columns
 *\param wantDebug             int more debug information on failure, if 1, else 0
 *\result succes               int reports success (1) or failure (0)
 */
 int elpa_cholesky_complex(int na, std::complex<double> *a, int lda, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
