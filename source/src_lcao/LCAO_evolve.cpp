#include "LCAO_evolve.h"

#include "../input.h"
#include "../module_base/scalapack_connector.h"
#include "../src_io/write_HS.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_pw/global.h"
#include "module_hamilt/hamilt_lcao.h"

#include <complex>
// fuxiang add 2016-10-28

Evolve_LCAO_Matrix::~Evolve_LCAO_Matrix()
{
}

inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

void Evolve_LCAO_Matrix::evolve_complex_matrix(const int& ik,
                                               hamilt::Hamilt* phami,
                                               Local_Orbital_wfc& lowf,
                                               psi::Psi<std::complex<double>>* psi,
                                               psi::Psi<std::complex<double>>* psi_laststep,
                                               double* ekb) const
{
    ModuleBase::TITLE("Evolve_LCAO_Matrix", "evolve_complex_matrix");
    time_t time_start = time(NULL);
    GlobalV::ofs_running << " Start Time : " << ctime(&time_start);

    if (INPUT.tddft == 1)
    {
#ifdef __MPI
        this->using_ScaLAPACK_complex(ik,
                                      phami,
                                      lowf.wfc_k_grid[ik],
                                      psi[0].get_pointer(),
                                      psi_laststep[0].get_pointer(),
                                      lowf,
                                      ekb);
        // this->using_ScaLAPACK_complex(ik, lowf.wfc_k_grid[ik], lowf.wfc_k[ik], lowf.wfc_k_laststep[ik], lowf, ekb);
#else
        this->using_LAPACK_complex(ik,
                                   phami,
                                   lowf.wfc_k_grid[ik],
                                   psi[0].get_pointer(),
                                   psi_laststep[0].get_pointer(),
                                   lowf,
                                   ekb);
#endif
    }
    else
    {
        ModuleBase::WARNING_QUIT("Evolve_LCAO_Matrix::evolve_complex_matrix", "only tddft==1 cando evolve");
    }

    time_t time_end = time(NULL);
    ModuleBase::GlobalFunc::OUT_TIME("evolve(std::complex)", time_start, time_end);

    return;
}

void Evolve_LCAO_Matrix::using_LAPACK_complex(const int& ik,
                                              hamilt::Hamilt* phami,
                                              std::complex<double>** wfc_k_grid,
                                              std::complex<double>* wfc_k,
                                              std::complex<double>* wfc_k_laststep,
                                              Local_Orbital_wfc& lowf,
                                              double* ekb) const
{
    ModuleBase::TITLE("Evolve_LCAO_Matrix", "using_LAPACK_complex");

    //	Calculate the U operator

    ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    ModuleBase::ComplexMatrix Stmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
    phami->matrix(h_mat, s_mat);
    for (int i = 0; i < GlobalV::NLOCAL; i++)
    {
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            Htmp(i, j) = h_mat.p[i * GlobalV::NLOCAL + j];
            //  Htmp(i,j) = (this->LM->Hloc2[i*GlobalV::NLOCAL+j] +this->LM->Hloc2_laststep[i*GlobalV::NLOCAL+j])/2.0;
            Stmp(i, j) = s_mat.p[i * GlobalV::NLOCAL + j];
        }
    }

    ModuleBase::ComplexMatrix wfc_tmp(GlobalV::NBANDS, GlobalV::NLOCAL, true);
    ModuleBase::ComplexMatrix wfc_laststep_tmp(GlobalV::NBANDS, GlobalV::NLOCAL, true);
    // wfc_laststep_tmp.c = wfc_k_laststep;

    for (int i = 0; i < GlobalV::NBANDS; i++)
    {
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            wfc_laststep_tmp.c[i * GlobalV::NLOCAL + j] = wfc_k_laststep[i * GlobalV::NLOCAL + j];
        }
    }

    /*
        GlobalV::ofs_running << " Htmp: " <<std::endl;
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            GlobalV::ofs_running << Htmp(i,j) <<"\t";
                    }
                    GlobalV::ofs_running <<std::endl;
            }
            GlobalV::ofs_running <<std::endl;

        GlobalV::ofs_running << " Stmp: " <<std::endl;
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            GlobalV::ofs_running << Stmp(i,j) <<"\t";
                    }
                    GlobalV::ofs_running <<std::endl;
            }
            GlobalV::ofs_running <<std::endl;
    */
    /*
        int INFO;

            int LWORK=3*GlobalV::NLOCAL-1; //tmp
            std::complex<double> * WORK = new std::complex<double>[LWORK];
            ModuleBase::GlobalFunc::ZEROS(WORK, LWORK);
            int IPIV[GlobalV::NLOCAL];

            LapackConnector::zgetrf( GlobalV::NLOCAL, GlobalV::NLOCAL, Stmp, GlobalV::NLOCAL, IPIV, &INFO);
            LapackConnector::zgetri( GlobalV::NLOCAL, Stmp, GlobalV::NLOCAL, IPIV, WORK, LWORK, &INFO);
    */
    /*
            std::cout << " S^-1: " <<std::endl;
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            std::cout << Stmp(i,j) <<"\t";
                    }
                    std::cout <<std::endl;
            }
            std::cout <<std::endl;
    */
    /*
        ModuleBase::ComplexMatrix S_plus_H(GlobalV::NLOCAL,GlobalV::NLOCAL);
        S_plus_H = Stmp*Htmp;
    */
    /*
            std::cout << " S^-1*H: " <<std::endl;
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            std::cout << S_plus_H(i,j) <<"\t";
                    }
                    std::cout <<std::endl;
            }
            std::cout <<std::endl;
    */

    ModuleBase::ComplexMatrix Denominator(GlobalV::NLOCAL, GlobalV::NLOCAL);
    /*
        for (int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for (int j=0; j<GlobalV::NLOCAL; j++)
                    {
                              Denominator(i,j) = std::complex<double>( -S_plus_H(i,j).imag(), S_plus_H(i,j).real() );
                    }
            }

            ModuleBase::ComplexMatrix Idmat(GlobalV::NLOCAL,GlobalV::NLOCAL);
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            if(i==j) Idmat(i,j) = std::complex<double>(1.0, 0.0);
                            else Idmat(i,j) = std::complex<double>(0.0, 0.0);
                    }
            }
    */
    ModuleBase::ComplexMatrix Numerator(GlobalV::NLOCAL, GlobalV::NLOCAL);
    //        Numerator = Idmat - 0.5*INPUT.mdp.md_dt*41.34*Denominator;
    //        Denominator = Idmat + 0.5*INPUT.mdp.md_dt*41.34*Denominator;

    complex<double> para = {0.0, 1.0};
    para = para * 0.25 * INPUT.mdp.md_dt;
    Numerator = Stmp - para * Htmp;
    Denominator = Stmp + para * Htmp;

    int info;
    int lwork = 3 * GlobalV::NLOCAL - 1; // tmp
    std::complex<double>* work = new std::complex<double>[lwork];
    ModuleBase::GlobalFunc::ZEROS(work, lwork);
    int ipiv[GlobalV::NLOCAL];

    LapackConnector::zgetrf(GlobalV::NLOCAL, GlobalV::NLOCAL, Denominator, GlobalV::NLOCAL, ipiv, &info);
    LapackConnector::zgetri(GlobalV::NLOCAL, Denominator, GlobalV::NLOCAL, ipiv, work, lwork, &info);

    ModuleBase::ComplexMatrix U_operator(GlobalV::NLOCAL, GlobalV::NLOCAL);
    /*
            GlobalV::ofs_running << " Numerator: " <<std::endl;
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            GlobalV::ofs_running << Numerator(i,j) <<"\t";
                    }
                    GlobalV::ofs_running <<std::endl;
            }
            GlobalV::ofs_running <<std::endl;

            GlobalV::ofs_running << " Denominator: " <<std::endl;
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            GlobalV::ofs_running << Denominator(i,j) <<"\t";
                    }
                    GlobalV::ofs_running <<std::endl;
            }
            GlobalV::ofs_running <<std::endl;
    */

    U_operator = Denominator * Numerator;
    /*
        GlobalV::ofs_running << "U_operator Success!!!" <<std::endl;
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            GlobalV::ofs_running<< U_operator(i,j).real()<<"+"<<U_operator(i,j).imag()<<"i ";
                    }
                    GlobalV::ofs_running <<std::endl;
            }
        GlobalV::ofs_running <<std::endl;
    */

    // Calculate wave function at t+delta t

    //	std::cout << "wave function coe at t+delta t !" << std::endl;

    /*	for(int i=0; i<GlobalV::NBANDS; i++)
        {
            for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                std::cout << wfc_k_grid[ik][i][j] << "\t";
            }
            std::cout <<std::endl;
        }
        std::cout << std::endl;
    */

    /*
            GlobalV::ofs_running<<"wfc_k_laststep "<<endl;
            for(int i=0; i<GlobalV::NBANDS; i++)
        {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                GlobalV::ofs_running << wfc_k_laststep.c[i*GlobalV::NLOCAL+j].real() <<
       "+"<<wfc_k_laststep.c[i*GlobalV::NLOCAL+j].imag()<<"i ";
            }
            GlobalV::ofs_running <<std::endl;
        }
        GlobalV::ofs_running <<std::endl;
    */

    const bool conjugate = false;
    wfc_tmp = wfc_laststep_tmp * transpose(U_operator, conjugate);

    ModuleBase::ComplexMatrix cmatrix(GlobalV::NBANDS, GlobalV::NBANDS);
    cmatrix = conj(wfc_tmp) * Stmp * transpose(wfc_tmp, conjugate);

    /*
           GlobalV::ofs_running<<"wfc_k before renomalization "<<endl;
           for(int i=0; i<GlobalV::NBANDS; i++)
       {
                   for(int j=0; j<GlobalV::NLOCAL; j++)
           {
               GlobalV::ofs_running << wfc_k.c[i*GlobalV::NLOCAL+j].real() <<
       "+"<<wfc_k.c[i*GlobalV::NLOCAL+j].imag()<<"i ";
           }
           GlobalV::ofs_running <<std::endl;
       }
       GlobalV::ofs_running <<std::endl;
   */

    for (int i = 0; i < GlobalV::NBANDS; i++)
    {
        double factor;
        factor = 1.0 / sqrt(cmatrix.c[i * GlobalV::NBANDS + i].real());
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            wfc_tmp.c[i * GlobalV::NLOCAL + j] *= factor;
        }
    }

    for (int i = 0; i < GlobalV::NBANDS; i++)
    {
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            wfc_k[i * GlobalV::NLOCAL + j] = wfc_tmp.c[i * GlobalV::NLOCAL + j];
            wfc_k_grid[i][j] = wfc_tmp.c[i * GlobalV::NLOCAL + j];
        }
    }
    /*
            GlobalV::ofs_running<<"wfc_k "<<endl;
            for(int i=0; i<GlobalV::NBANDS; i++)
        {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                GlobalV::ofs_running << wfc_k.c[i*GlobalV::NLOCAL+j].real() <<
       "+"<<wfc_k.c[i*GlobalV::NLOCAL+j].imag()<<"i ";
            }
            GlobalV::ofs_running <<std::endl;
        }
        GlobalV::ofs_running <<std::endl;
    */

    /*
            GlobalV::ofs_running<<"wfc_k_grid "<<endl;
            for(int i=0; i<GlobalV::NBANDS; i++)
        {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                GlobalV::ofs_running << wfc_k_grid[i][j].real() << "+"<<wfc_k_grid[i][j].imag()<<"i ";
            }
            GlobalV::ofs_running <<std::endl;
        }
        GlobalV::ofs_running <<std::endl;
    */

    ///*
    // calculate energy level
    ModuleBase::ComplexMatrix Ematrix(GlobalV::NLOCAL, GlobalV::NLOCAL);
    Ematrix = conj(wfc_tmp) * Htmp * transpose(wfc_tmp, conjugate);
    for (int i = 0; i < GlobalV::NBANDS; i++)
    {
        ekb[i] = Ematrix.c[i * GlobalV::NBANDS + i].real();
    }
    //*/

    /*
            GlobalV::ofs_running<<endl;
            GlobalV::ofs_running<<"print ekb : "<<endl;
            for(int ib=0; ib<GlobalV::NBANDS; ++ib)
            {
                    //GlobalV::ofs_running<<"ekb[" << ib+1 << "]  " << ekb[ib] << std::endl;
                    GlobalV::ofs_running<<"ekb[" << ib+1 << "]  " << ekb[ib]*13.605693 << std::endl;
            }
            GlobalV::ofs_running<<endl;
    */
    /*
    cout<<"E matrix"<<endl;
    for(int i=0; i<GlobalV::NBANDS; i++)
            {
                    for(int j=0; j<GlobalV::NBANDS; j++)
                    {
                            std::cout <<
    Ematrix.c[i*GlobalV::NBANDS+j].real()<<"+"<<Ematrix.c[i*GlobalV::NLOCAL+j].imag()<<"i  ";
                    }
                    std::cout <<std::endl;
            }
            std::cout <<std::endl;
    */

    /*	for(int i=0; i<GlobalV::NBANDS; i++)
        {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                std::cout << wfc_k_grid[ik][i][j] << "\t";
            }
            std::cout <<std::endl;
        }
    */

    /*      std::cout << " c: " <<std::endl;
            for(int i=0; i<GlobalV::NLOCAL; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            std::cout << wfc_k_grid[ik][i][j] <<"\t";
                    }
                    std::cout <<std::endl;
            }
            std::cout <<std::endl;
    */
    /*
        for(int i=0; i<GlobalV::NBANDS; i++)
            {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                            std::cout << wfc_k_grid[ik][i][j] << "\t";
                    }
                    std::cout <<std::endl;
            }
            std::cout << std::endl;
    */

    delete[] work;
    //	delete[] ipiv;

    return;
}

#ifdef __MPI
void Evolve_LCAO_Matrix::using_ScaLAPACK_complex(const int& ik,
                                                 hamilt::Hamilt* phami,
                                                 std::complex<double>** wfc_k_grid,
                                                 std::complex<double>* wfc_k,
                                                 std::complex<double>* wfc_k_laststep,
                                                 Local_Orbital_wfc& lowf,
                                                 double* ekb) const
{
    ModuleBase::TITLE("Evolve_LCAO_Matrix", "using_ScaLAPACK_complex");

    // inverse of matrix
    // pzgetrf (int *m, int *n, Complex16 *a, int ia, int ja, int *desca, int *ipiv, int info);
    // pzgetri (int *n, Complex16 *a, int *ia, int ja, int *desca, int *ipiv, Complex16 *Work, int *lwork, int *iwork,
    // int *liwork, int *info);

    // product of vector and matrix
    // pzgemv(const char *trans, const int *m, const int *n, const Complex16 *alpha, const Complex16 *a,
    //	const int *ia, const int *ja, const int *desca, const Complex16*x, const int *ix, const int *jx,
    //	const int *descx, const int *incx, const Complex16 *beta, Complex16 *y, const int *iy,
    //	const int *jy, const int *descy, const int *incy);

    // matrix-matrix sum
    // pzgeadd (const char *trans, const int *m, const int *n, const Complex16 *alpha, const Complex16 *a, const int
    // *ia, 	const int *ja, const int *desca, const Complex16 *beta, Complex16 *c, const int *ic, const int *jc,
    // const
    // int *descc);

    // matrix-matrix product
    // pzgemm

    // cout << "begin1: " <<endl;

    char uplo = 'U';
    const int inc = 1;

    const int one_int = 1;
    int print_matrix = 0;
    int nrow = this->ParaV->nrow;
    int ncol = this->ParaV->ncol;
    int ncol_bands = this->ParaV->ncol_bands;
    // cout<<"ncol_bands="<<ncol_bands<<" ncol="<<ncol<<" nrow="<<nrow<<endl;
    // int nprocs, myid;
    // MPI_status status;
    // MPI_Comm_size(comm_2D, &nprocs);
    // MPI_Comm_rank(comm_2D, &myid);

    int loc_pos;
    // complex<double>* Stmp = this->LM->Sdiag2;
    // complex<double>* Htmp1 = this->LM->Hdiag2;
    // complex<double>* Htmp2 = this->LM->Sdiag2;
    // complex<double>* Htmp3 = this->LM->Sdiag2;

    complex<double>* Stmp = new complex<double>[this->ParaV->nloc];
    complex<double>* Htmp1 = new complex<double>[this->ParaV->nloc];
    complex<double>* Htmp2 = new complex<double>[this->ParaV->nloc];
    complex<double>* Htmp3 = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Stmp, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(Htmp1, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(Htmp2, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(Htmp3, this->ParaV->nloc);

    // cout << "this->LM->Hloc2" << *this->LM->Hloc2 << endl;
    // cout << "*Htmp2: " << *Htmp2 << endl;

    double* eigen = new double[GlobalV::NLOCAL];
    ModuleBase::GlobalFunc::ZEROS(eigen, GlobalV::NLOCAL);

    // cout << "GlobalV::NLOCAL : " << GlobalV::NLOCAL << endl;
    // cout << "nloc : " << this->ParaV->nloc << endl;

    // cout << "begin02:" <<endl;

    // assert(loc_size > 0);
    // complex<double>* Z = new complex<double>[this->loc_size * NLOCAL];
    // ModuleBase::GlobalFunc::ZEROS(Z, this->loc_size * NLOCAL);

    hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
    phami->matrix(h_mat, s_mat);
    // cout<<"h_mat : "<<h_mat.p[0]<<" "<<h_mat.p[1]<<endl;

    zcopy_(&this->ParaV->nloc, s_mat.p, &inc, Stmp, &inc);
    zcopy_(&this->ParaV->nloc, h_mat.p, &inc, Htmp1, &inc);
    zcopy_(&this->ParaV->nloc, h_mat.p, &inc, Htmp2, &inc);
    zcopy_(&this->ParaV->nloc, h_mat.p, &inc, Htmp3, &inc);

    // cout << "*Htmp2: " << *Htmp2 << endl;
    complex<double> imag = {0.0, 1.0};

    if (print_matrix)
    {
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " S matrix :" << endl;
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                GlobalV::ofs_running << Stmp[i * ncol + j].real() << "+" << Stmp[i * ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " H matrix :" << endl;
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                GlobalV::ofs_running << Htmp1[i * ncol + j].real() << "+" << Htmp1[i * ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
    }

    complex<double> alpha = {1.0, 0.0};
    char transa = 'N';
    int desca = 0;
    complex<double> beta = {0.0, -0.25 * INPUT.mdp.md_dt}; // this need modify
    int descc = 0;

    pzgeadd_(&transa,
             &GlobalV::NLOCAL,
             &GlobalV::NLOCAL,
             &alpha,
             Stmp,
             &one_int,
             &one_int,
             this->ParaV->desc,
             &beta,
             Htmp1,
             &one_int,
             &one_int,
             this->ParaV->desc);

    // beta = (0.0, 0.5)*INPUT.md_dt;
    beta = {0.0, 0.25 * INPUT.mdp.md_dt}; // this need modify
    // cout<<"dt="<<INPUT.mdp.md_dt<<endl;

    // cout << "*Htmp1: " << *Htmp1 << endl;
    // cout << "Htmp2: " << *Htmp2 << endl;

    pzgeadd_(&transa,
             &GlobalV::NLOCAL,
             &GlobalV::NLOCAL,
             &alpha,
             Stmp,
             &one_int,
             &one_int,
             this->ParaV->desc,
             &beta,
             Htmp2,
             &one_int,
             &one_int,
             this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " fenmu:" << endl;
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                GlobalV::ofs_running << Htmp2[i * ncol + j].real() << "+" << Htmp2[i * ncol + j].imag() << "i ";
                // GlobalV::ofs_running<<i<<" "<<j<<" "<<Htmp3[i*GlobalV::NLOCAL+j]<<endl;
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
    }

    // Next, invert the denominator
    int* ipiv = new int[this->ParaV->nloc];
    int info;

    // cout << "*Htmp2: " << *Htmp2 << endl;
    // cout << "begin04:" << endl;

    pzgetrf_(&GlobalV::NLOCAL, &GlobalV::NLOCAL, Htmp2, &one_int, &one_int, this->ParaV->desc, ipiv, &info);

    int LWORK = -1, liWORK = -1;
    std::vector<std::complex<double>> WORK(1, 0);
    std::vector<int> iWORK(1, 0);

    // cout << "begin05:" << endl;

    pzgetri_(&GlobalV::NLOCAL,
             Htmp2,
             &one_int,
             &one_int,
             this->ParaV->desc,
             ipiv,
             WORK.data(),
             &LWORK,
             iWORK.data(),
             &liWORK,
             &info);

    LWORK = WORK[0].real();
    WORK.resize(LWORK, 0);
    liWORK = iWORK[0];
    iWORK.resize(liWORK, 0);

    pzgetri_(&GlobalV::NLOCAL,
             Htmp2,
             &one_int,
             &one_int,
             this->ParaV->desc,
             ipiv,
             WORK.data(),
             &LWORK,
             iWORK.data(),
             &liWORK,
             &info);

    // alpha = (1.0, 0.0);
    // beta = (0.0, 0.0);
    char transb = 'T'; // This place requires subsequent testing of different transb.
    int descb = 0;

    double alpha_1[2] = {1.0, 0.0};
    double beta_1[2] = {0.0, 0.0};

    // cout << "*Htmp2: " << *Htmp2 << endl;
    // cout << "begin06:" << endl;

    pzgemm_(&transa,
            &transb,
            &GlobalV::NLOCAL,
            &GlobalV::NLOCAL,
            &GlobalV::NLOCAL,
            &alpha_1[0],
            Htmp2,
            &one_int,
            &one_int,
            this->ParaV->desc,
            Htmp1,
            &one_int,
            &one_int,
            this->ParaV->desc,
            &beta_1[0],
            Htmp3,
            &one_int,
            &one_int,
            this->ParaV->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << " fenmu^-1:" << endl;
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                GlobalV::ofs_running << Htmp2[i * ncol + j].real() << "+" << Htmp2[i * ncol + j].imag() << "i ";
                // GlobalV::ofs_running<<i<<" "<<j<<" "<<Htmp3[i*GlobalV::NLOCAL+j]<<endl;
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " fenzi:" << endl;
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                GlobalV::ofs_running << Htmp1[i * ncol + j].real() << "+" << Htmp1[i * ncol + j].imag() << "i ";
                // GlobalV::ofs_running<<i<<" "<<j<<" "<<Htmp3[i*GlobalV::NLOCAL+j]<<endl;
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " U operator:" << endl;
        for (int i = 0; i < nrow; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                double aa, bb;
                aa = Htmp3[i * ncol + j].real();
                bb = Htmp3[i * ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
                // Htmp3[i*ncol+j]={1,0};
                //GlobalV::ofs_running << Htmp3[i * ncol + j].real() << "+" << Htmp3[i * ncol + j].imag() << "i ";
                // GlobalV::ofs_running<<i<<" "<<j<<" "<<Htmp3[i*GlobalV::NLOCAL+j]<<endl;
            }
            GlobalV::ofs_running << endl;
        }
    }

    // cout << "U_operator Success!!!" <<endl;

    transa = 'T';
    transb = 'T';
    pzgemm_(&transa,
            &transb,
            &GlobalV::NBANDS,
            &GlobalV::NLOCAL,
            &GlobalV::NLOCAL,
            &alpha_1[0],
            wfc_k_laststep,
            &one_int,
            &one_int,
            this->ParaV->desc_wfc,
            Htmp3,
            &one_int,
            &one_int,
            this->ParaV->desc,
            &beta_1[0],
            // wfc_k.c, &one_int, &one_int, this->ParaV->desc_wfc);
            Htmp2,
            &one_int,
            &one_int,
            this->ParaV->desc);

    pztranu_(&GlobalV::NLOCAL,
             &GlobalV::NLOCAL,
             &alpha_1[0],
             Htmp2,
             &one_int,
             &one_int,
             this->ParaV->desc,
             &beta_1[0],
             Htmp3,
             &one_int,
             &one_int,
             this->ParaV->desc);
    zcopy_(&this->ParaV->nloc_wfc, Htmp3, &inc, wfc_k, &inc);

    // renormalize wfc_k
    complex<double>* tmp1 = new complex<double>[this->ParaV->nloc];
    complex<double>* tmp2 = new complex<double>[this->ParaV->nloc];
    complex<double>* tmp3 = new complex<double>[this->ParaV->nloc_wfc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(tmp2, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(tmp3, this->ParaV->nloc_wfc);
    const char N_char = 'N', T_char = 'T', C_char = 'C';
    const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};
    pzgemm_(&T_char,
            &N_char,
            &GlobalV::NBANDS,
            &GlobalV::NLOCAL,
            &GlobalV::NLOCAL,
            &one_float[0],
            wfc_k,
            &one_int,
            &one_int,
            this->ParaV->desc_wfc,
            Stmp,
            &one_int,
            &one_int,
            this->ParaV->desc,
            &zero_float[0],
            tmp1,
            &one_int,
            &one_int,
            this->ParaV->desc);
    pztranu_(&GlobalV::NLOCAL,
             &GlobalV::NLOCAL,
             &one_float[0],
             tmp1,
             &one_int,
             &one_int,
             this->ParaV->desc,
             &zero_float[0],
             tmp2,
             &one_int,
             &one_int,
             this->ParaV->desc);
    zcopy_(&this->ParaV->nloc_wfc, tmp2, &inc, tmp3, &inc);
    // ModuleBase::ComplexMatrix tmp4 = conj(wfc_k);
    complex<double>* Cij = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Cij, this->ParaV->nloc);
    pzgemm_(&C_char,
            //&T_char,
            &N_char,
            &GlobalV::NBANDS,
            &GlobalV::NBANDS,
            &GlobalV::NLOCAL,
            &one_float[0],
            // tmp4.c,
            wfc_k,
            &one_int,
            &one_int,
            this->ParaV->desc_wfc,
            tmp3,
            &one_int,
            &one_int,
            this->ParaV->desc_wfc,
            &zero_float[0],
            Cij,
            &one_int,
            &one_int,
            this->ParaV->desc);

    int myid;
    MPI_Comm_rank(this->ParaV->comm_2D, &myid);
    int naroc[2]; // maximum number of row or column
    for (int iprow = 0; iprow < this->ParaV->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < this->ParaV->dim1; ++ipcol)
        {
            const int coord[2] = {iprow, ipcol};
            int src_rank;
            info = MPI_Cart_rank(this->ParaV->comm_2D, coord, &src_rank);
            if (myid == src_rank)
            {
                naroc[0] = this->ParaV->nrow;
                naroc[1] = this->ParaV->ncol;
                //}
                // info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, this->ParaV->comm_2D);
                // info = MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, uhm.LM->ParaV->comm_2D);
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, this->ParaV->nb, this->ParaV->dim1, ipcol);
                    if (igcol >= GlobalV::NBANDS)
                        continue;
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, this->ParaV->nb, this->ParaV->dim0, iprow);
                        if (igrow >= GlobalV::NBANDS)
                            continue;
                        if (igcol == igrow)
                        {
                            Cij[j * naroc[0] + i] = {1.0 / sqrt(Cij[j * naroc[0] + i].real()), 0.0};
                        }
                        else
                        {
                            Cij[j * naroc[0] + i] = {0.0, 0.0};
                        }
                        // info = MPI_Bcast(&GlobalC::wf.ekb[ik][igcol], 1, MPI_DOUBLE, src_rank,
                        // uhm.LM->ParaV->comm_2D);
                    }
                }
            }
        } // loop ipcol
    } // loop iprow
    ///*
    pzgemm_(&T_char,
            &N_char,
            &GlobalV::NLOCAL,
            &GlobalV::NBANDS,
            &GlobalV::NBANDS,
            &one_float[0],
            // wfc_k.c,&one_int,&one_int,this->ParaV->desc_wfc,
            Htmp2,
            &one_int,
            &one_int,
            this->ParaV->desc,
            Cij,
            &one_int,
            &one_int,
            this->ParaV->desc,
            &zero_float[0],
            wfc_k,
            &one_int,
            &one_int,
            this->ParaV->desc_wfc);
    // Htmp3,&one_int,&one_int,this->ParaV->desc);
    // zcopy_(&this->ParaV->nloc_wfc, Htmp3, &inc, wfc_k.c, &inc);
    //*/
    /*
    pzgemv_(
        &transa,
        &GlobalV::NLOCAL, &GlobalV::NLOCAL,
        &alpha_1[0],
        Htmp3, &one_int, &one_int, this->ParaV->desc,
        wfc_k_laststep.c, &one_int, &one_int, this->ParaV->desc_wfc, &one_int,
        &beta_1[0],
        wfc_k.c, &one_int, &one_int, this->ParaV->desc_wfc, &one_int
        );*/

    if (print_matrix)
    {
        GlobalV::ofs_running << " Cij:" << endl;
        for (int i = 0; i < ncol; i++)
        {
            for (int j = 0; j < nrow; j++)
            {
                GlobalV::ofs_running << Cij[i * ncol + j].real() << "+" << Cij[i * ncol + j].imag() << "i ";
                // GlobalV::ofs_running<<i<<" "<<j<<" "<<Htmp3[i*GlobalV::NLOCAL+j]<<endl;
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " wfc_k_laststep:" << endl;
        for (int i = 0; i < ncol_bands; i++)
        {
            for (int j = 0; j < nrow; j++)
            {
                double aa, bb;
                aa = wfc_k_laststep[i * ncol + j].real();
                bb = wfc_k_laststep[i * ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
                //GlobalV::ofs_running << wfc_k_laststep[i * ncol + j].real() << "+"
                //                     << wfc_k_laststep[i * ncol + j].imag() << "i ";
                // GlobalV::ofs_running<<i<<" "<<j<<" "<<Htmp3[i*GlobalV::NLOCAL+j]<<endl;
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " wfc_k:" << endl;
        for (int i = 0; i < ncol_bands; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                double aa, bb;
                aa = wfc_k[i * ncol + j].real();
                bb = wfc_k[i * ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
                //GlobalV::ofs_running << wfc_k[i * ncol + j].real() << "+" << wfc_k[i * ncol + j].imag() << "i ";
                // GlobalV::ofs_running<<i<<" "<<j<<" "<<Htmp3[i*GlobalV::NLOCAL+j]<<endl;
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " wfc_k nlocal*nlocal:" << endl;
        for (int i = 0; i < ncol; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                double aa, bb;
                aa = Htmp2[i * ncol + j].real();
                bb = Htmp2[i * ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
                //GlobalV::ofs_running << Htmp2[i * ncol + j].real() << "+" << Htmp2[i * ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " wfc_k nlocal*nlocal transpose:" << endl;
        for (int i = 0; i < ncol; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                double aa, bb;
                aa = Htmp3[i * ncol + j].real();
                bb = Htmp3[i * ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
                // GlobalV::ofs_running << Htmp3[i * ncol + j].real() << "+" << Htmp3[i * ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
    }

    // set out_wfc_lcao=0 temporarily
    int zero = 0;
    // lowf.wfc_2d_to_grid(zero, wfc_k.c, wfc_k_grid, zero);

    if (print_matrix)
    {
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " wfc_k_grid" << endl;
        for (int i = 0; i < ncol_bands; i++)
        {
            for (int j = 0; j < ncol; j++)
            {
                GlobalV::ofs_running << wfc_k_grid[i][j].real() << "+" << wfc_k_grid[i][j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
    }

    // calculate ekb
    complex<double>* Htmp = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(tmp2, this->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(tmp3, this->ParaV->nloc_wfc);
    ModuleBase::GlobalFunc::ZEROS(Htmp, this->ParaV->nloc);
    zcopy_(&this->ParaV->nloc, h_mat.p, &one_int, Htmp, &one_int);
    pzgemm_(&T_char,
            &N_char,
            &GlobalV::NBANDS,
            &GlobalV::NLOCAL,
            &GlobalV::NLOCAL,
            &one_float[0],
            wfc_k,
            &one_int,
            &one_int,
            this->ParaV->desc_wfc,
            Htmp,
            &one_int,
            &one_int,
            this->ParaV->desc,
            &zero_float[0],
            tmp1,
            &one_int,
            &one_int,
            this->ParaV->desc);
    pztranu_(&GlobalV::NLOCAL,
             &GlobalV::NLOCAL,
             &one_float[0],
             tmp1,
             &one_int,
             &one_int,
             this->ParaV->desc,
             &zero_float[0],
             tmp2,
             &one_int,
             &one_int,
             this->ParaV->desc);
    zcopy_(&this->ParaV->nloc_wfc, tmp2, &inc, tmp3, &inc);
    // ModuleBase::ComplexMatrix tmp5 = conj(wfc_k);
    complex<double>* Eij = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Eij, this->ParaV->nloc);
    pzgemm_(&C_char,
            &N_char,
            &GlobalV::NBANDS,
            &GlobalV::NBANDS,
            &GlobalV::NLOCAL,
            &one_float[0],
            // tmp5.c,
            wfc_k,
            &one_int,
            &one_int,
            this->ParaV->desc_wfc,
            tmp3,
            &one_int,
            &one_int,
            this->ParaV->desc_wfc,
            &zero_float[0],
            Eij,
            &one_int,
            &one_int,
            this->ParaV->desc);
    double* Eii = new double[GlobalV::NBANDS];
    for (int i = 0; i < GlobalV::NBANDS; i++)
        Eii[i] = 0.0;
    for (int iprow = 0; iprow < this->ParaV->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < this->ParaV->dim1; ++ipcol)
        {
            const int coord[2] = {iprow, ipcol};
            int src_rank;
            info = MPI_Cart_rank(this->ParaV->comm_2D, coord, &src_rank);
            if (myid == src_rank)
            {
                naroc[0] = this->ParaV->nrow;
                naroc[1] = this->ParaV->ncol;
                //}
                // info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, this->ParaV->comm_2D);
                // info = MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, uhm.LM->ParaV->comm_2D);
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, this->ParaV->nb, this->ParaV->dim1, ipcol);
                    if (igcol >= GlobalV::NBANDS)
                        continue;
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, this->ParaV->nb, this->ParaV->dim0, iprow);
                        if (igrow >= GlobalV::NBANDS)
                            continue;
                        if (igcol == igrow)
                        {
                            Eii[igcol] = Eij[j * naroc[0] + i].real();
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow
    info = MPI_Allreduce(Eii, GlobalC::wf.ekb[ik], GlobalV::NBANDS, MPI_DOUBLE, MPI_SUM, this->ParaV->comm_2D);

    // the eigenvalues.
    // dcopy_(&NBANDS, eigen, &inc, ekb, &inc);
    delete[] Stmp;
    delete[] Htmp1;
    delete[] Htmp2;
    delete[] Htmp3;
    delete[] tmp1;
    delete[] tmp2;
    delete[] tmp3;
    delete[] Htmp;
    delete[] Eij;
    delete[] Eii;
    delete[] Cij;
    delete[] eigen;
    delete[] ipiv;

    // Z is delete in gath_eig
    // ModuleBase::timer::tick("Evolve_LCAO_Matrix","gath_eig_complex",'G');
    return;
}
#endif
