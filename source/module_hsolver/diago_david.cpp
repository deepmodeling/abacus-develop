#include "diago_david.h"

#include "diago_iter_assist.h"
#include "module_base/constants.h"
#include "module_base/blas_connector.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "src_parallel/parallel_common.h"

namespace hsolver
{
typedef hamilt::Operator<std::complex<double>>::hpsi_info hp_info;

int DiagoDavid::PW_DIAG_NDIM = 4;

DiagoDavid::DiagoDavid(const double *precondition_in)
{
    this->precondition = precondition_in;

    test_david = 2;
    // 1: check which function is called and which step is executed
    // 2: check the eigenvalues of the result of each iteration
    // 3: check the eigenvalues and errors of the last result
    // default: no check
}

void DiagoDavid::diag_mock(hamilt::Hamilt* phm_in, psi::Psi<std::complex<double>> &psi, double *eigenvalue_in)
{
    if (test_david == 1)
        ModuleBase::TITLE("DiagoDavid", "diag_mock");
    ModuleBase::timer::tick("DiagoDavid", "diag_mock");

    const int dim = psi.get_current_nbas();
    const int nband = psi.get_nbands();

    assert(DiagoDavid::PW_DIAG_NDIM > 1);
    assert(DiagoDavid::PW_DIAG_NDIM * nband < dim * GlobalV::NPROC_IN_POOL);
    // qianrui change it 2021-7-25.
    // In strictly speaking, it shoule be PW_DIAG_NDIM*nband < npw sum of all pools. We roughly estimate it here.
    // However, in most cases, total number of plane waves should be much larger than nband*PW_DIAG_NDIM

    int nbase_x = DiagoDavid::PW_DIAG_NDIM * nband; // maximum dimension of the reduced basis set

    psi::Psi<std::complex<double>> basis(1, nbase_x, dim, &(psi.get_ngk(0))); // the reduced basis set
    ModuleBase::ComplexMatrix hp(nbase_x, dim); // the product of H and psi in the reduced basis set
    ModuleBase::ComplexMatrix sp(nbase_x, dim); // the Product of S and psi in the reduced basis set

    ModuleBase::ComplexMatrix hc(nbase_x, nbase_x); // Hamiltonian on the reduced basis
    ModuleBase::ComplexMatrix sc(nbase_x, nbase_x); // Overlap on the reduced basis
    ModuleBase::ComplexMatrix vc(nbase_x, nbase_x); // Eigenvectors of hc
    std::vector<double> eigenvalue(nbase_x); // the lowest N eigenvalues of hc

    std::vector<bool> convflag(nband, false); // convflag[m] = true if the m th band is convergent
    std::vector<int> unconv(nband); // unconv[m] store the number of the m th unconvergent band

    int nbase = 0; // the dimension of the reduced basis set
    this->notconv = nband; // the number of the unconvergent bands
    for (int m = 0; m < nband; m++)
        unconv[m] = m;

    ModuleBase::timer::tick("DiagoDavid", "first");
    // orthogonalise the initial trial psi(0~nband-1)
    //plan for SchmitOrth
    ModuleBase::ComplexMatrix lagrange_matrix(nband, nband);
    std::vector<int> pre_matrix_mm_m(nband, 0);
    std::vector<int> pre_matrix_mv_m(nband, 1);
    this->planSchmitOrth(nband, pre_matrix_mm_m.data(), pre_matrix_mv_m.data());
    for( int m = 0; m < nband; m++)
    {
        phm_in->sPsi(&psi(m, 0), &sp(m, 0),  (size_t)dim);
    }
    //begin SchmitOrth
    for (int m = 0; m < nband; m++)
    {
        // psi_m = psi(m)
        ModuleBase::GlobalFunc::COPYARRAY(&psi(m, 0), &basis(m, 0), dim);
        /*for (int ig = 0; ig < dim; ig++)
        {
            psi_m[ig] = psi(m, ig);
        }*/

        //phm_in->sPsi(psi_m.data(), spsi.data(),  (size_t)dim);
        this->SchmitOrth(
            dim,
            nband,
            m,
            basis,
            sp,
            &lagrange_matrix(m, 0),
            pre_matrix_mm_m[m],
            pre_matrix_mv_m[m]
        );
        phm_in->sPsi(&basis(m, 0), &sp(m, 0),  (size_t)dim);

        // basis(m) = psi_m, hp(m) = H |psi_m>, sp(m) = S |psi_m>
        //ModuleBase::GlobalFunc::COPYARRAY(psi_m.data(), &basis(m, 0), dim);
        //ModuleBase::GlobalFunc::COPYARRAY(spsi.data(), &sp(m, 0), dim);
        /*std::complex<double>* sp_p = &sp(m, 0);
        std::complex<double>* basis_p = &basis(m, 0);
        for (int ig = 0; ig < dim; ig++)
        {
            basis_p[ig] = psi_m[ig];
            //hp(m, ig) = hpsi[ig];
            sp_p[ig] = spsi[ig];
        }*/
    }
    //end of SchmitOrth and calculate H|psi>
    hp_info dav_hpsi_in(&basis, psi::Range(1, 0, 0, nband-1), &hp(0, 0));
    phm_in->ops->hPsi(dav_hpsi_in);

    hc.zero_out();
    sc.zero_out();

    this->cal_elem(dim, nbase, this->notconv, basis, hp, sp, hc, sc);

    this->diag_zhegvx(nbase, nband, hc, sc, nbase_x, eigenvalue.data(), vc);

    for (int m = 0; m < nband; m++)
    {
        eigenvalue_in[m] = eigenvalue[m];
    }

    ModuleBase::timer::tick("DiagoDavid", "first");

    int dav_iter = 0;
    do
    {
        dav_iter++;

        this->cal_grad(phm_in,
                       dim,
                       nbase,
                       this->notconv,
                       basis,
                       hp,
                       sp,
                       vc,
                       unconv.data(),
                       eigenvalue.data());

        this->cal_elem(dim, nbase, this->notconv, basis, hp, sp, hc, sc);

        this->diag_zhegvx(nbase, nband, hc, sc, nbase_x, eigenvalue.data(), vc);

        // check convergence and update eigenvalues
        ModuleBase::timer::tick("DiagoDavid", "check_update");

        this->notconv = 0;
        for (int m = 0; m < nband; m++)
        {
            convflag[m] = (abs(eigenvalue[m] - eigenvalue_in[m]) < DiagoIterAssist::PW_DIAG_THR);

            if (!convflag[m])
            {
                unconv[this->notconv] = m;
                this->notconv++;
            }

            eigenvalue_in[m] = eigenvalue[m];
        }

        ModuleBase::timer::tick("DiagoDavid", "check_update");
        if (!this->notconv || (nbase + this->notconv > nbase_x) || (dav_iter == DiagoIterAssist::PW_DIAG_NMAX))
        {
            ModuleBase::timer::tick("DiagoDavid", "last");

            // updata eigenvectors of Hamiltonian
            ModuleBase::GlobalFunc::ZEROS(psi.get_pointer(), psi.get_nbands() * psi.get_nbasis());
            char transa = 'N';
            char transb = 'T';
            zgemm_(&transa,
                    &transb,
                    &dim, // m: row of A,C
                    &nband, // n: col of B,C
                    &nbase, // k: col of A, row of B
                    &ModuleBase::ONE, // alpha
                    basis.get_pointer(), // A
                    &basis.get_nbasis(), // LDA: if(N) max(1,m) if(T) max(1,k)
                    vc.c, // B
                    &nbase_x, // LDB: if(N) max(1,k) if(T) max(1,n)
                    &ModuleBase::ZERO, // belta
                    psi.get_pointer(), // C
                    &psi.get_nbasis()); // LDC: if(N) max(1, m)
            /*for (int m = 0; m < nband; m++)
            {
                for (int j = 0; j < nbase; j++)
                {
                    for (int ig = 0; ig < dim; ig++)
                    {
                        psi(m, ig) += vc(j, m) * basis(j, ig);
                    }
                }
            }*/

            if (!this->notconv || (dav_iter == DiagoIterAssist::PW_DIAG_NMAX))
            {
                // overall convergence or last iteration: exit the iteration

                ModuleBase::timer::tick("DiagoDavid", "last");
                break;
            }
            else
            {
                // if the dimension of the reduced basis set is becoming too large,
                // then replace the first N (=nband) basis vectors with the current
                // estimate of the eigenvectors and set the basis dimension to N;

                this->refresh(dim, nband, nbase, eigenvalue_in, psi, basis, hp, sp, hc, sc, vc);
                ModuleBase::timer::tick("DiagoDavid", "last");
            }

        } // end of if

    } while (1);

    DiagoIterAssist::avg_iter += static_cast<double>(dav_iter);

    ModuleBase::timer::tick("DiagoDavid", "diag_mock");
    return;
}

void DiagoDavid::cal_grad(hamilt::Hamilt* phm_in,
                          const int &npw,
                          const int &nbase, // current dimension of the reduced basis
                          const int &notconv,
                          psi::Psi<std::complex<double>> &basis,
                          ModuleBase::ComplexMatrix &hp,
                          ModuleBase::ComplexMatrix &sp,
                          const ModuleBase::ComplexMatrix &vc,
                          const int *unconv,
                          const double *eigenvalue)
{
    if (test_david == 1)
        ModuleBase::TITLE("DiagoDavid", "cal_grad");
    if(notconv == 0) return;
    ModuleBase::timer::tick("DiagoDavid", "cal_grad");

    //use template pointer for accelerate
    std::complex<double> *hpsi;
    std::complex<double> *spsi;
    std::complex<double> *ppsi;

    // expand the reduced basis set with the new basis vectors P|R(psi)>...
    // in which psi are the last eigenvectors
    // we define |R(psi)> as (H-ES)*|Psi>, E = <psi|H|psi>/<psi|S|psi>
    ModuleBase::ComplexMatrix vc_ev_vector(notconv, nbase);
    for (int m = 0; m < notconv; m++)
    {
        for(int i = 0; i < nbase; i++)
        {
            vc_ev_vector(m, i) = vc(i, unconv[m]);
        }
    }
    ppsi = &basis(nbase, 0);
    int inc = 1;
    char trans = 'N';
    char transb = 'N';
    zgemm_(&trans,
            &transb,
            &npw, // m: row of A,C
            &notconv, // n: col of B,C
            &nbase, // k: col of A, row of B
            &ModuleBase::ONE, // alpha
            hp.c, // A
            &hp.nc, // LDA: if(N) max(1,m) if(T) max(1,k)
            vc_ev_vector.c, // B
            &vc_ev_vector.nc, // LDB: if(N) max(1,k) if(T) max(1,n)
            &ModuleBase::ZERO, // belta
            ppsi, // C
            &basis.get_nbasis()); // LDC: if(N) max(1, m)
    /*zgemv_(&trans,
        &npw,
        &nbase,
        &ModuleBase::ONE,
        hp.c,
        &hp.nc,
        vc_ev_vector.data(),
        &inc,
        &ModuleBase::ZERO,
        respsi,
        &inc);*/
    for (int m = 0; m < notconv; m++)
    {
        for(int i = 0; i < nbase; i++)
        {
            vc_ev_vector(m, i) *= -1 * eigenvalue[unconv[m]];
        }
    }
    zgemm_(&trans,
            &transb,
            &npw, // m: row of A,C
            &notconv, // n: col of B,C
            &nbase, // k: col of A, row of B
            &ModuleBase::ONE, // alpha
            sp.c, // A
            &sp.nc, // LDA: if(N) max(1,m) if(T) max(1,k)
            vc_ev_vector.c, // B
            &vc_ev_vector.nc, // LDB: if(N) max(1,k) if(T) max(1,n)
            &ModuleBase::ONE, // belta
            ppsi, // C
            &basis.get_nbasis()); // LDC: if(N) max(1, m)
        /*zgemv_(&trans,
            &npw,
            &nbase,
            &ModuleBase::ONE,
            sp.c,
            &sp.nc,
            vc_ev_vector.data(),
            &inc,
            &ModuleBase::ONE,
            respsi,
            &inc);*/

        /*ModuleBase::GlobalFunc::ZEROS(respsi, npw);
        for (int i = 0; i < nbase; i++)
        {
            hpsi = &(hp(i, 0));
            spsi = &(sp(i, 0));
            auto vc_value = vc(i, unconv[m]);
            auto ev_value = eigenvalue[unconv[m]];
            for (int ig = 0; ig < npw; ig++)
            {
                respsi[ig] += vc_value * (hpsi[ig] - ev_value * spsi[ig]);
            }
        }*/

    for (int m = 0; m < notconv; m++)
    {
        ppsi = &basis(nbase + m, 0);
        for (int ig = 0; ig < npw; ig++)
        {
            ppsi[ig] /= this->precondition[ig];
        }
    }

    //there is a nbase to nbase + notconv band orthogonalise
    //plan for SchmitOrth
    ModuleBase::ComplexMatrix lagrange_matrix(notconv, nbase + notconv);
    std::vector<int> pre_matrix_mm_m(notconv, 0);
    std::vector<int> pre_matrix_mv_m(notconv, 1);
    this->planSchmitOrth(notconv, pre_matrix_mm_m.data(), pre_matrix_mv_m.data());
    for( int m = 0; m < notconv; m++)
    {
        phm_in->sPsi(&basis(nbase + m, 0), &sp(nbase + m, 0),  (size_t)npw);
    }
    //first nbase bands psi* dot notconv bands spsi to prepare lagrange_matrix
    trans = 'C';
    transb = 'N';
    //calculate the square matrix for future lagranges
    zgemm_(&trans,
            &transb,
            &nbase, // m: row of A,C
            &notconv, // n: col of B,C
            &npw, // k: col of A, row of B
            &ModuleBase::ONE, // alpha
            &basis(0, 0), // A
            &basis.get_nbasis(), // LDA: if(N) max(1,m) if(T) max(1,k)
            &sp(nbase, 0), // B
            &sp.nc, // LDB: if(N) max(1,k) if(T) max(1,n)
            &ModuleBase::ZERO, // belta
            &lagrange_matrix(0, 0), // C
            &lagrange_matrix.nc); // LDC: if(N) max(1, m)

    for (int m = 0; m < notconv; m++)
    {
        ppsi = &basis(nbase + m, 0);
        spsi = &sp(nbase + m, 0);

        this->SchmitOrth(
            npw,
            nbase + notconv,
            nbase + m,
            basis,
            sp,
            &lagrange_matrix(m, 0),
            pre_matrix_mm_m[m],
            pre_matrix_mv_m[m]
        );
        phm_in->sPsi(ppsi, spsi, (size_t)npw);

    }
    //calculate H|psi> for not convergence bands
    hp_info dav_hpsi_in(&basis, psi::Range(1, 0, nbase, nbase + notconv-1), &hp(nbase, 0));
    phm_in->ops->hPsi(dav_hpsi_in);
    
    ModuleBase::timer::tick("DiagoDavid", "cal_grad");
    return;
}

void DiagoDavid::cal_elem(const int &npw,
                          int &nbase, // current dimension of the reduced basis
                          const int &notconv, // number of newly added basis vectors
                          const psi::Psi<std::complex<double>> &basis,
                          const ModuleBase::ComplexMatrix &hp,
                          const ModuleBase::ComplexMatrix &sp,
                          ModuleBase::ComplexMatrix &hc,
                          ModuleBase::ComplexMatrix &sc)
{
    if (test_david == 1)
        ModuleBase::TITLE("DiagoDavid", "cal_elem");

    if(notconv == 0) return;
    ModuleBase::timer::tick("DiagoDavid", "cal_elem");

    // updat the reduced Hamiltonian
    int offset_h = nbase * hc.nr;
    int offset_s = nbase * sc.nr;
    //	ModuleBase::GlobalFunc::ZEROS( hc.c+offset_h, notconv*hc.nr );
    //	ModuleBase::GlobalFunc::ZEROS( sc.c+offset_s, notconv*sc.nr );

    char trans1 = 'C';
    char trans2 = 'N';
    const int nb_notc = (nbase + notconv);
    hc = transpose(hc, false);
    zgemm_(&trans1,
           &trans2,
           &notconv,
           &nb_notc,
           &npw,
           &ModuleBase::ONE,
           &basis(nbase, 0),
           &basis.get_nbasis(),
           hp.c,
           &hp.nc,
           &ModuleBase::ONE,
           hc.c + nbase,
           &hc.nr);
    hc = transpose(hc, false);

    sc = transpose(sc, false);
    zgemm_(&trans1,
           &trans2,
           &notconv,
           &nb_notc,
           &npw,
           &ModuleBase::ONE,
           &basis(nbase, 0),
           &basis.get_nbasis(),
           sp.c,
           &sp.nc,
           &ModuleBase::ONE,
           sc.c + nbase,
           &sc.nr);
    sc = transpose(sc, false);

//    Parallel_Reduce::reduce_complex_double_pool(hc.c + offset_h, notconv * hc.nr);
//    Parallel_Reduce::reduce_complex_double_pool(sc.c + offset_s, notconv * sc.nr);
#ifdef __MPI
	std::complex<double> *swap = new std::complex<double>[notconv * hc.nc];
    ModuleBase::GlobalFunc::COPYARRAY(hc.c+offset_h, swap, notconv * hc.nc);
    MPI_Reduce(swap,hc.c+offset_h,notconv * hc.nc,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD);
    ModuleBase::GlobalFunc::COPYARRAY(sc.c+offset_h, swap, notconv * hc.nc);
    MPI_Reduce(swap,sc.c+offset_h,notconv * hc.nc,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD);
	delete[] swap;
#endif
    /*
        for( int i = nbase; i < nbase+notconv; i++ )
        {
            for( int j = 0; j <i; j++ )
            {
                hc(j,i) = conj( hc(i,j) );
                sc(j,i) = conj( sc(i,j) );
            }
        }
    */
    nbase += notconv;
    ModuleBase::timer::tick("DiagoDavid", "cal_elem");
    return;
}

//==============================================================================
// optimize diag_zhegvx().

// 09-05-09 wangjp
// fixed a bug in diag_zhegvx().
// modify the dimension of h and s as (n,n) and copy the leading N*N
// part of hc & sc into h & s

// 09-05-10 wangjp
// As the complexmatrixs will be copied again in the subroutine ZHEGVX(...  ),
// i.e ZHEGVX(...) will not destroy the input complexmatrixs,
// we needn't creat another two complexmatrixs in diag_zhegvx().
//==============================================================================
void DiagoDavid::diag_zhegvx(const int &n,
                             const int &m,
                             const ModuleBase::ComplexMatrix &hc,
                             const ModuleBase::ComplexMatrix &sc,
                             const int &ldh,
                             double *eigenvalue,
                             ModuleBase::ComplexMatrix &vc)
{
    //	ModuleBase::TITLE("DiagoDavid","diag_zhegvx");
    ModuleBase::timer::tick("DiagoDavid", "diag_zhegvx");
if(GlobalV::MY_RANK == 0)
{
    assert(ldh >= max(1, n));
    int lwork;
    int info = 0;
    int mm = m;
    std::string name1 = "ZHETRD";
    std::string name2 = "L";

    int nb = LapackConnector::ilaenv(1, name1.c_str(), name2.c_str(), n, -1, -1, -1);
    if (nb < 1)
    {
        nb = max(1, n);
    }
    if (nb == 1 || nb >= n)
    {
        lwork = 2 * n; // qianrui fix a bug 2021-7-25 : lwork should be at least max(1,2*n)
    }
    else
    {
        lwork = (nb + 1) * n;
    }
    std::complex<double> *work = new std::complex<double>[2 * lwork];
    assert(work != 0);
    double *rwork = new double[7 * n];
    assert(rwork != 0);
    int *iwork = new int[5 * n];
    assert(iwork != 0);
    int *ifail = new int[n];
    assert(ifail != 0);
    ModuleBase::GlobalFunc::ZEROS(work, lwork); // qianrui change it, only first lwork numbers are used in zhegvx
    ModuleBase::GlobalFunc::ZEROS(rwork, 7 * n);
    ModuleBase::GlobalFunc::ZEROS(iwork, 5 * n);
    ModuleBase::GlobalFunc::ZEROS(ifail, n);

    // ModuleBase::WARNING_QUIT("divid","open zhegvx!");

    LapackConnector::zhegvx(1,
                            'V',
                            'I',
                            'L',
                            n,
                            hc,
                            n,
                            sc,
                            n,
                            0.0,
                            0.0,
                            1,
                            m,
                            0.0,
                            mm,
                            eigenvalue,
                            vc,
                            n,
                            work,
                            lwork,
                            rwork,
                            iwork,
                            ifail,
                            info);

    delete[] work;
    delete[] rwork;
    delete[] iwork;
    delete[] ifail;
}

#ifdef __MPI
    for (int i=0;i<n;i++)
    {
        Parallel_Common::bcast_complex_double(&vc(i, 0), m);
    }
    Parallel_Common::bcast_double(eigenvalue,m);
#endif

    ModuleBase::timer::tick("DiagoDavid", "diag_zhegvx");
    return;
}

void DiagoDavid::refresh(const int &npw,
                         const int &nband,
                         int &nbase,
                         const double *eigenvalue_in,
                         const psi::Psi<std::complex<double>> &psi,
                         psi::Psi<std::complex<double>> &basis,
                         ModuleBase::ComplexMatrix &hp,
                         ModuleBase::ComplexMatrix &sp,
                         ModuleBase::ComplexMatrix &hc,
                         ModuleBase::ComplexMatrix &sc,
                         ModuleBase::ComplexMatrix &vc)
{
    if (test_david == 1)
        ModuleBase::TITLE("DiagoDavid", "refresh");
    ModuleBase::timer::tick("DiagoDavid", "refresh");

    // update hp,sp
    basis.zero_out();
    char transa = 'N';
    char transb = 'T';
    zgemm_(&transa,
            &transb,
            &npw, // m: row of A,C
            &nband, // n: col of B,C
            &nbase, // k: col of A, row of B
            &ModuleBase::ONE, // alpha
            hp.c, // A
            &hp.nc, // LDA: if(N) max(1,m) if(T) max(1,k)
            vc.c, // B
            &vc.nc, // LDB: if(N) max(1,k) if(T) max(1,n)
            &ModuleBase::ZERO, // belta
            basis.get_pointer(), // C
            &basis.get_nbasis()); // LDC: if(N) max(1, m)

    zgemm_(&transa,
            &transb,
            &npw, // m: row of A,C
            &nband, // n: col of B,C
            &nbase, // k: col of A, row of B
            &ModuleBase::ONE, // alpha
            sp.c, // A
            &sp.nc, // LDA: if(N) max(1,m) if(T) max(1,k)
            vc.c, // B
            &vc.nc, // LDB: if(N) max(1,k) if(T) max(1,n)
            &ModuleBase::ZERO, // belta
            &basis(nband, 0), // C
            &basis.get_nbasis()); // LDC: if(N) max(1, m)

    /*for (int m = 0; m < nband; m++)
    {
        for (int j = 0; j < nbase; j++)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                basis(m, ig) += vc(j, m) * hp(j, ig);
                basis(m + nband, ig) += vc(j, m) * sp(j, ig);
            }
        }
    }*/

    ModuleBase::GlobalFunc::COPYARRAY(&basis(0, 0), &hp(0, 0), npw * nband);
    ModuleBase::GlobalFunc::COPYARRAY(&basis(nband, 0), &sp(0, 0), npw * nband);
    /*for (int m = 0; m < nband; m++)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            hp(m, ig) = basis(m, ig);
            sp(m, ig) = basis(m + nband, ig);
        }
    }*/

    // update basis
    basis.zero_out();
    for (int m = 0; m < nband; m++)
    {
        ModuleBase::GlobalFunc::COPYARRAY(&psi(m, 0), &basis(m, 0), npw);
        /*for (int ig = 0; ig < npw; ig++)
            basis(m, ig) = psi(m, ig);*/
    }

    // updata the reduced Hamiltonian
    nbase = nband;
    hc.zero_out();
    sc.zero_out();
    for (int i = 0; i < nbase; i++)
    {
        hc(i, i) = eigenvalue_in[i];
        sc(i, i) = ModuleBase::ONE;
        vc(i, i) = ModuleBase::ONE;
    }

    ModuleBase::timer::tick("DiagoDavid", "refresh");
    return;
}

void DiagoDavid::cal_err(const int &npw,
                         const int &nband,
                         const int &nbase,
                         const ModuleBase::ComplexMatrix &vc,
                         const ModuleBase::ComplexMatrix &hp,
                         const psi::Psi<std::complex<double>> &basis,
                         const double *eigenvalue_in,
                         std::complex<double> *respsi)
{
    ModuleBase::timer::tick("DiagoDavid", "cal_err");
    double *err = new double[nband];
    assert(err != 0);

    for (int m = 0; m < nband; m++)
    {
        ModuleBase::GlobalFunc::ZEROS(respsi, npw);
        for (int j = 0; j < nbase; j++)
        {
            for (int ig = 0; ig < npw; ig++)
            {
                respsi[ig] += vc(j, m) * (hp(j, ig) - eigenvalue_in[m] * basis(j, ig));
            }
        }

        err[m] = ModuleBase::GlobalFunc::ddot_real(npw, respsi, respsi);
        err[m] = sqrt(err[m]);
    }

    std::cout << "i       eigenvalues       err (||*||)   " << std::endl;

    for (int i = 0; i < nband; i++)
    {
        std::cout << i << "\t\t" << eigenvalue_in[i] << "\t\t" << err[i] << std::endl;
    }

    delete[] err;
    ModuleBase::timer::tick("DiagoDavid", "cal_err");
    return;
}

void DiagoDavid::SchmitOrth(const int &npw,
                            const int n_band,
                            const int m,
                            psi::Psi<std::complex<double>>& psi,
                            const ModuleBase::ComplexMatrix& spsi,
                            std::complex<double>* lagrange_m,
                            const int mm_size,
                            const int mv_size)
{
    //	if(test_david == 1) ModuleBase::TITLE("DiagoDavid","SchmitOrth");
    ModuleBase::timer::tick("DiagoDavid", "SchmitOrth");

    // orthogonalize starting eigenfunction to those already calculated
    // psi_m orthogonalize to psi(0) ~ psi(m-1)
    // Attention, the orthogonalize here read as
    // psi(m) -> psi(m) - \sum_{i < m} \langle psi(i)|S|psi(m) \rangle psi(i)
    // so the orthogonalize is performed about S.

    assert(psi.get_nbands() >= n_band);
    assert(m >= 0);
    assert(m < n_band);

    std::complex<double>* psi_m = &psi(m, 0);

    //std::complex<double> *lagrange = new std::complex<double>[m + 1];
    //ModuleBase::GlobalFunc::ZEROS(lagrange, m + 1);

    int inc = 1;
    char trans = 'C';
    char transb = 'N';
    //calculate the square matrix for future lagranges
    if(mm_size != 0)
    {
        zgemm_(&trans,
            &transb,
            &mm_size, // m: row of A,C
            &mm_size, // n: col of B,C
            &npw, // k: col of A, row of B
            &ModuleBase::ONE, // alpha
            &psi(m-mv_size+1-mm_size, 0), // A
            &psi.get_nbasis(), // LDA: if(N) max(1,m) if(T) max(1,k)
            &spsi(m, 0), // B
            &spsi.nc, // LDB: if(N) max(1,k) if(T) max(1,n)
            &ModuleBase::ZERO, // belta
            &lagrange_m[m-mv_size+1-mm_size], // C
            &n_band); // LDC: if(N) max(1, m)
    }
    //calculate other lagranges for this band
    zgemv_(&trans,
           &npw,
           &mv_size,
           &ModuleBase::ONE,
           &psi(m-mv_size+1, 0),
           &psi.get_nbasis(),
           &spsi(m, 0),
           &inc,
           &ModuleBase::ZERO,
           &lagrange_m[m-mv_size+1],
           &inc);
    /*for (int j = 0; j < m; j++)
    {
        const std::complex<double>* psi_p = &(psi(j, 0));
        zdotc_(&lagrange[j], &npw, psi_p, &one, spsi, &one);
        for (int ig = 0; ig < npw; ig++)
        {
            lagrange[j] += conj(psi(j, ig)) * spsi[ig];
        }
        //	lagrange[j] = Diago_CG::ddot( npw, psi, j, spsi );
    }*/
    //zdotc_(&lagrange[m], &npw, psi_m, &inc, spsi, &inc);
    /*for (int ig = 0; ig < npw; ig++)
    {
        lagrange[m] += conj(psi_m[ig]) * spsi[ig];
    }*/
    //	lagrange[m] = Diago_CG::ddot( npw, psi_m, spsi );

    Parallel_Reduce::reduce_complex_double_pool(lagrange_m, m + 1);

    //	out.printr1_d("lagrange", lagrange, m+1 );

    double psi_norm = lagrange_m[m].real();
    assert(psi_norm > 0.0);
    //	std::cout << "m = " << m << std::endl;

    for (int j = 0; j < m; j++)
    {
        const std::complex<double> alpha = std::complex<double>(-1, 0) * lagrange_m[j];
        zaxpy_(&npw, &alpha, &psi(j,0), &inc, psi_m, &inc);
        /*for (int ig = 0; ig < npw; ig++)
        {
            psi_m[ig] -= lagrange[j] * psi(j, ig);
        }*/
        psi_norm -= (conj(lagrange_m[j]) * lagrange_m[j]).real();
    }

    assert(psi_norm > 0.0);

    psi_norm = sqrt(psi_norm);

    if (psi_norm < 1.0e-12)
    {
        std::cout << "DiagoDavid::SchmitOrth:aborted for psi_norm <1.0e-12" << std::endl;
        std::cout << "n_band = " << n_band << std::endl;
        std::cout << "m = " << m << std::endl;
        exit(0);
    }
    else
    {
        for (int i = 0; i < npw; i++)
        {
            psi_m[i] /= psi_norm;
        }
    }

    //delete[] lagrange;
    ModuleBase::timer::tick("DiagoDavid", "SchmitOrth");
    return;
}

void DiagoDavid::planSchmitOrth(
    const int nband,
    int* pre_matrix_mm_m,
    int* pre_matrix_mv_m
)
{
    if(nband<=0)return;
    ModuleBase::GlobalFunc::ZEROS(pre_matrix_mm_m, nband);
    ModuleBase::GlobalFunc::ZEROS(pre_matrix_mv_m, nband);
    int last_matrix_size = nband;
    int matrix_size = int(nband / 2);
    int divide_times = 0;
    std::vector<int> divide_points(nband);
    int res_nband = nband - matrix_size;
    while(matrix_size>1)
    {
        int index = nband - matrix_size;
        if(divide_times == 0)
        {
            divide_points[0] = index;
            pre_matrix_mm_m[index] = matrix_size;
            if(res_nband == matrix_size) pre_matrix_mv_m[index] = 1;
            else pre_matrix_mv_m[index] = 2;
            divide_times = 1;
        }
        else
        {
        for(int i=divide_times-1; i>=0; i--)
        {
            divide_points[i*2] = divide_points[i] - matrix_size;
            divide_points[i*2+1] = divide_points[i*2] + last_matrix_size;
            pre_matrix_mm_m[ divide_points[i*2] ] = matrix_size;
            pre_matrix_mm_m[ divide_points[i*2+1]] = matrix_size;
            if(res_nband == matrix_size)
            {
                pre_matrix_mv_m[divide_points[i*2]] = 1;
                pre_matrix_mv_m[divide_points[i*2+1]] = 1;
            }
            else
            {
                pre_matrix_mv_m[divide_points[i*2]] = 2;
                pre_matrix_mv_m[divide_points[i*2+1]] = 2;
            }
        }
            divide_times *= 2;
        }
        last_matrix_size = matrix_size;
        matrix_size = int(res_nband / 2);
        res_nband -= matrix_size;
    }
    //fill the pre_matrix_mv_m array
    pre_matrix_mv_m[0] = 1;
    for(int m = 1; m < nband; m++)
    {
        if(pre_matrix_mv_m[m] == 0)
        {
            pre_matrix_mv_m[m] = pre_matrix_mv_m[m-1]+1;
        }
    }
}

void DiagoDavid::diag(hamilt::Hamilt *phm_in, psi::Psi<std::complex<double>> &psi, double *eigenvalue_in)
{
    /// record the times of trying iterative diagonalization
    int ntry = 0;
    this->notconv = 0;
    do
    {
        this->diag_mock(phm_in, psi, eigenvalue_in);
        ++ntry;
    } while (DiagoIterAssist::test_exit_cond(ntry, this->notconv));

    if (notconv > max(5, psi.get_nbands() / 4))
    {
        std::cout << "\n notconv = " << this->notconv;
        std::cout << "\n DiagoDavid::diag', too many bands are not converged! \n";
    }
    return;
}

} // namespace hsolver
