#include "diago_iter_assist.h"

#include "module_base/blas_connector.h"
#include "module_base/complexmatrix.h"
#include "module_base/constants.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"

namespace hsolver
{

double DiagoIterAssist::avg_iter = 0.0;
int DiagoIterAssist::PW_DIAG_NMAX = 30;
double DiagoIterAssist::PW_DIAG_THR = 1.0e-2;
bool DiagoIterAssist::need_subspace = false;

//----------------------------------------------------------------------
// Hamiltonian diagonalization in the subspace spanned
// by nstart states psi (atomic or random wavefunctions).
// Produces on output n_band eigenvectors (n_band <= nstart) in evc.
//----------------------------------------------------------------------
void DiagoIterAssist::diagH_subspace(hamilt::Hamilt* pHamilt,
                                     const psi::Psi<std::complex<double>> &psi,
                                     psi::Psi<std::complex<double>> &evc,
                                     double *en,
                                     int n_band)
{
    ModuleBase::TITLE("DiagoIterAssist", "diagH_subspace");
    ModuleBase::timer::tick("DiagoIterAssist", "diagH_subspace");

    // two case:
    // 1. pw base: nstart = n_band, psi(nbands * npwx)
    // 2. lcao_in_pw base: nstart >= n_band, psi(NLOCAL * npwx)
    const int nstart = psi.get_nbands();
    if (n_band == 0)
        n_band = nstart;
    assert(n_band <= nstart);

    ModuleBase::ComplexMatrix hc(nstart, nstart);
    ModuleBase::ComplexMatrix sc(nstart, nstart);
    ModuleBase::ComplexMatrix hvec(nstart, n_band);

    const int dmin = psi.get_current_nbas();
    const int dmax = psi.get_nbasis();

    // qianrui improve this part 2021-3-14
    std::complex<double> *aux = new std::complex<double>[dmax * nstart];
    const std::complex<double> *paux = aux;
    const std::complex<double> *ppsi = psi.get_pointer();

    pHamilt->hPsi(ppsi, aux, size_t(dmax * nstart));

    char trans1 = 'C';
    char trans2 = 'N';
    zgemm_(&trans1,
           &trans2,
           &nstart,
           &nstart,
           &dmin,
           &ModuleBase::ONE,
           ppsi,
           &dmax,
           paux,
           &dmax,
           &ModuleBase::ZERO,
           hc.c,
           &nstart);
    hc = transpose(hc, false);

    zgemm_(&trans1,
           &trans2,
           &nstart,
           &nstart,
           &dmin,
           &ModuleBase::ONE,
           ppsi,
           &dmax,
           ppsi,
           &dmax,
           &ModuleBase::ZERO,
           sc.c,
           &nstart);
    sc = transpose(sc, false);

    delete[] aux;

    if (GlobalV::NPROC_IN_POOL > 1)
    {
        Parallel_Reduce::reduce_complex_double_pool(hc.c, nstart * nstart);
        Parallel_Reduce::reduce_complex_double_pool(sc.c, nstart * nstart);
    }

    // after generation of H and S matrix, diag them
    DiagoIterAssist::diagH_LAPACK(nstart, n_band, hc, sc, nstart, en, hvec);

    //=======================
    // diagonize the H-matrix
    //=======================
    if ((GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw") && GlobalV::CALCULATION == "nscf")
    {
        GlobalV::ofs_running << " Not do zgemm to get evc." << std::endl;
    }
    else if ((GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw")
             && (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md"
                 || GlobalV::CALCULATION == "relax")) // pengfei 2014-10-13
    {
        // because psi and evc are different here,
        // I think if psi and evc are the same,
        // there may be problems, mohan 2011-01-01
        char transa = 'N';
        char transb = 'T';
        zgemm_(&transa,
               &transb,
               &dmax, // m: row of A,C
               &n_band, // n: col of B,C
               &nstart, // k: col of A, row of B
               &ModuleBase::ONE, // alpha
               ppsi, // A
               &dmax, // LDA: if(N) max(1,m) if(T) max(1,k)
               hvec.c, // B
               &n_band, // LDB: if(N) max(1,k) if(T) max(1,n)
               &ModuleBase::ZERO, // belta
               evc.get_pointer(), // C
               &dmax); // LDC: if(N) max(1, m)
    }
    else
    {
        // As the evc and psi may refer to the same matrix, we first
        // create a temporary matrix to store the result. (by wangjp)
        // qianrui improve this part 2021-3-13
        char transa = 'N';
        char transb = 'T';
        ModuleBase::ComplexMatrix evctmp(n_band, dmin, false);
        zgemm_(&transa,
               &transb,
               &dmin,
               &n_band,
               &nstart,
               &ModuleBase::ONE,
               ppsi,
               &dmax,
               hvec.c,
               &n_band,
               &ModuleBase::ZERO,
               evctmp.c,
               &dmin);
        for (int ib = 0; ib < n_band; ib++)
        {
            for (int ig = 0; ig < dmin; ig++)
            {
                evc(ib, ig) = evctmp(ib, ig);
            }
        }
    }

    ModuleBase::timer::tick("DiagoIterAssist", "diagH_subspace");
    return;
}

void DiagoIterAssist::diagH_subspace_init(hamilt::Hamilt* pHamilt,
                                     const ModuleBase::ComplexMatrix &psi,
                                     psi::Psi<std::complex<double>> &evc,
                                     double *en)
{
    ModuleBase::TITLE("DiagoIterAssist", "diagH_subspace_init");
    ModuleBase::timer::tick("DiagoIterAssist", "diagH_subspace");

    // two case:
    // 1. pw base: nstart = n_band, psi(nbands * npwx)
    // 2. lcao_in_pw base: nstart >= n_band, psi(NLOCAL * npwx)
    const int nstart = psi.nr;
    const int n_band = evc.get_nbands();

    ModuleBase::ComplexMatrix hc(nstart, nstart);
    ModuleBase::ComplexMatrix sc(nstart, nstart);
    ModuleBase::ComplexMatrix hvec(nstart, n_band);

    const int dmin = evc.get_current_nbas();
    const int dmax = evc.get_nbasis();

    // qianrui improve this part 2021-3-14
    std::complex<double> *aux = new std::complex<double>[dmax * nstart];
    const std::complex<double> *paux = aux;
    const std::complex<double> *ppsi = psi.c;

    pHamilt->hPsi(ppsi, aux, size_t(dmax * nstart));

    char trans1 = 'C';
    char trans2 = 'N';
    zgemm_(&trans1,
           &trans2,
           &nstart,
           &nstart,
           &dmin,
           &ModuleBase::ONE,
           ppsi,
           &dmax,
           paux,
           &dmax,
           &ModuleBase::ZERO,
           hc.c,
           &nstart);
    hc = transpose(hc, false);

    zgemm_(&trans1,
           &trans2,
           &nstart,
           &nstart,
           &dmin,
           &ModuleBase::ONE,
           ppsi,
           &dmax,
           ppsi,
           &dmax,
           &ModuleBase::ZERO,
           sc.c,
           &nstart);
    sc = transpose(sc, false);

    delete[] aux;

    if (GlobalV::NPROC_IN_POOL > 1)
    {
        Parallel_Reduce::reduce_complex_double_pool(hc.c, nstart * nstart);
        Parallel_Reduce::reduce_complex_double_pool(sc.c, nstart * nstart);
    }

    // after generation of H and S matrix, diag them
    DiagoIterAssist::diagH_LAPACK(nstart, n_band, hc, sc, nstart, en, hvec);

    //=======================
    // diagonize the H-matrix
    //=======================
    if ((GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw") && GlobalV::CALCULATION == "nscf")
    {
        GlobalV::ofs_running << " Not do zgemm to get evc." << std::endl;
    }
    else if ((GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw")
             && (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md"
                 || GlobalV::CALCULATION == "relax")) // pengfei 2014-10-13
    {
        // because psi and evc are different here,
        // I think if psi and evc are the same,
        // there may be problems, mohan 2011-01-01
        char transa = 'N';
        char transb = 'T';
        zgemm_(&transa,
               &transb,
               &dmax, // m: row of A,C
               &n_band, // n: col of B,C
               &nstart, // k: col of A, row of B
               &ModuleBase::ONE, // alpha
               ppsi, // A
               &dmax, // LDA: if(N) max(1,m) if(T) max(1,k)
               hvec.c, // B
               &n_band, // LDB: if(N) max(1,k) if(T) max(1,n)
               &ModuleBase::ZERO, // belta
               evc.get_pointer(), // C
               &dmax); // LDC: if(N) max(1, m)
    }
    else
    {
        // As the evc and psi may refer to the same matrix, we first
        // create a temporary matrix to store the result. (by wangjp)
        // qianrui improve this part 2021-3-13
        char transa = 'N';
        char transb = 'T';
        ModuleBase::ComplexMatrix evctmp(n_band, dmin, false);
        zgemm_(&transa,
               &transb,
               &dmin,
               &n_band,
               &nstart,
               &ModuleBase::ONE,
               ppsi,
               &dmax,
               hvec.c,
               &n_band,
               &ModuleBase::ZERO,
               evctmp.c,
               &dmin);
        for (int ib = 0; ib < n_band; ib++)
        {
            for (int ig = 0; ig < dmin; ig++)
            {
                evc(ib, ig) = evctmp(ib, ig);
            }
        }
    }

    ModuleBase::timer::tick("DiagoIterAssist", "diagH_subspace");
    return;
}

void DiagoIterAssist::diagH_LAPACK(const int nstart,
                                   const int nbands,
                                   const ModuleBase::ComplexMatrix &hc,
                                   const ModuleBase::ComplexMatrix &sc,
                                   const int ldh, // nstart
                                   double *e,
                                   ModuleBase::ComplexMatrix &hvec)
{
    ModuleBase::TITLE("DiagoIterAssist", "diagH_LAPACK");
    ModuleBase::timer::tick("DiagoIterAssist", "diagH_LAPACK");

    int lwork = 0;

    ModuleBase::ComplexMatrix sdum(nstart, ldh);
    ModuleBase::ComplexMatrix hdum;

    sdum = sc;

    const bool all_eigenvalues = (nstart == nbands);

    // workspace query
    int nb = LapackConnector::ilaenv(1, "ZHETRD", "U", nstart, -1, -1, -1);

    if (nb < 1)
    {
        nb = std::max(1, nstart);
    }

    if (nb == 1 || nb >= nstart)
    {
        lwork = 2 * nstart; // mohan modify 2009-08-02
    }
    else
    {
        lwork = (nb + 1) * nstart;
    }

    std::complex<double> *work = new std::complex<double>[lwork];
    ModuleBase::GlobalFunc::ZEROS(work, lwork);

    //=====================================================================
    // input s and (see below) h are copied so that they are not destroyed
    //=====================================================================

    int info = 0;
    int rwork_dim;
    if (all_eigenvalues)
    {
        rwork_dim = 3 * nstart - 2;
    }
    else
    {
        rwork_dim = 7 * nstart;
    }

    double *rwork = new double[rwork_dim];
    ModuleBase::GlobalFunc::ZEROS(rwork, rwork_dim);

    if (all_eigenvalues)
    {
        //===========================
        // calculate all eigenvalues
        //===========================
        hvec = hc;
        LapackConnector::zhegv(1, 'V', 'U', nstart, hvec, ldh, sdum, ldh, e, work, lwork, rwork, info);
    }
    else
    {
        //=====================================
        // calculate only m lowest eigenvalues
        //=====================================
        int *iwork = new int[5 * nstart];
        int *ifail = new int[nstart];

        ModuleBase::GlobalFunc::ZEROS(rwork, 7 * nstart);
        ModuleBase::GlobalFunc::ZEROS(iwork, 5 * nstart);
        ModuleBase::GlobalFunc::ZEROS(ifail, nstart);

        hdum.create(nstart, ldh);
        hdum = hc;

        //=============================
        // Number of calculated bands
        //=============================
        int mm = nbands;

        LapackConnector::zhegvx(1, // INTEGER
                                'V', // CHARACTER*1
                                'I', // CHARACTER*1
                                'U', // CHARACTER*1
                                nstart, // INTEGER
                                hdum, // COMPLEX*16 array
                                ldh, // INTEGER
                                sdum, // COMPLEX*16 array
                                ldh, // INTEGER
                                0.0, // DOUBLE PRECISION
                                0.0, // DOUBLE PRECISION
                                1, // INTEGER
                                nbands, // INTEGER
                                0.0, // DOUBLE PRECISION
                                mm, // INTEGER
                                e, // DOUBLE PRECISION array
                                hvec, // COMPLEX*16 array
                                ldh, // INTEGER
                                work, // DOUBLE array, dimension (MAX(1,LWORK))
                                lwork, // INTEGER
                                rwork, // DOUBLE PRECISION array, dimension (7*N)
                                iwork, // INTEGER array, dimension (5*N)
                                ifail, // INTEGER array, dimension (N)
                                info // INTEGER
        );

        delete[] iwork;
        delete[] ifail;
    }
    delete[] rwork;
    delete[] work;

    ModuleBase::timer::tick("DiagoIterAssist", "diagH_LAPACK");
    return;
}

bool DiagoIterAssist::test_exit_cond(const int &ntry, const int &notconv)
{
    //================================================================
    // If this logical function is true, need to do diagH_subspace
    // and cg again.
    //================================================================

    bool scf = true;
    if (GlobalV::CALCULATION == "nscf")
        scf = false;

    // If ntry <=5, try to do it better, if ntry > 5, exit.
    const bool f1 = (ntry <= 5);

    // In non-self consistent calculation, do until totally converged.
    const bool f2 = ((!scf && (notconv > 0)));

    // if self consistent calculation, if not converged > 5,
    // using diagH_subspace and cg method again. ntry++
    const bool f3 = ((scf && (notconv > 5)));
    return (f1 && (f2 || f3));
}

} // namespace hsolver