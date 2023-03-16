#include "LCAO_evolve.h"

#include "module_io/input.h"
#include "module_base/scalapack_connector.h"
#include "module_base/lapack_connector.h"
#include "module_io/write_HS.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "ELEC_evolve.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"

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
                                               hamilt::Hamilt<double>* p_hamilt,
                                               psi::Psi<std::complex<double>>* psi_k,
                                               psi::Psi<std::complex<double>>* psi_k_laststep,
                                               double* ekb) const
{
    ModuleBase::TITLE("Evolve_LCAO_Matrix", "evolve_complex_matrix");
    time_t time_start = time(NULL);
    GlobalV::ofs_running << " Start Time : " << ctime(&time_start);

    if (GlobalV::ESOLVER_TYPE == "tddft")
    {
#ifdef __MPI
        this->using_ScaLAPACK_complex(GlobalV::NBANDS, GlobalV::NLOCAL, psi_k_laststep[0].get_pointer(), p_hamilt, psi_k[0].get_pointer(), ekb);
#else
        this->using_LAPACK_complex(ik, p_hamilt, psi_k[0].get_pointer(), psi_k_laststep[0].get_pointer(), ekb);
#endif
    }
    else
    {
        ModuleBase::WARNING_QUIT("Evolve_LCAO_Matrix::evolve_complex_matrix",
                                 "only esolver_type == tddft cando evolve");
    }

    time_t time_end = time(NULL);
    ModuleBase::GlobalFunc::OUT_TIME("evolve(std::complex)", time_start, time_end);

    return;
}

void Evolve_LCAO_Matrix::using_LAPACK_complex(const int& ik,
                                              hamilt::Hamilt<double>* p_hamilt,
                                              std::complex<double>* psi_k,
                                              std::complex<double>* psi_k_laststep,
                                              double* ekb) const
{
    ModuleBase::TITLE("Evolve_LCAO_Matrix", "using_LAPACK_complex");

    //	Calculate the U operator

    ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    ModuleBase::ComplexMatrix Stmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
    p_hamilt->matrix(h_mat, s_mat);
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
    // wfc_laststep_tmp.c = psi_k_laststep;

    for (int i = 0; i < GlobalV::NBANDS; i++)
    {
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            wfc_laststep_tmp.c[i * GlobalV::NLOCAL + j] = psi_k_laststep[i * GlobalV::NLOCAL + j];
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

            int lwork=3*GlobalV::NLOCAL-1; //tmp
            std::complex<double> * work = new std::complex<double>[lwork];
            ModuleBase::GlobalFunc::ZEROS(work, lwork);
            int IPIV[GlobalV::NLOCAL];

            LapackConnector::zgetrf( GlobalV::NLOCAL, GlobalV::NLOCAL, Stmp, GlobalV::NLOCAL, IPIV, &INFO);
            LapackConnector::zgetri( GlobalV::NLOCAL, Stmp, GlobalV::NLOCAL, IPIV, work, lwork, &INFO);
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
    int* ipiv = new int[GlobalV::NLOCAL];

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

    /*
            GlobalV::ofs_running<<"psi_k_laststep "<<endl;
            for(int i=0; i<GlobalV::NBANDS; i++)
        {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                GlobalV::ofs_running << psi_k_laststep.c[i*GlobalV::NLOCAL+j].real() <<
       "+"<<psi_k_laststep.c[i*GlobalV::NLOCAL+j].imag()<<"i ";
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
           GlobalV::ofs_running<<"psi_k before renomalization "<<endl;
           for(int i=0; i<GlobalV::NBANDS; i++)
       {
                   for(int j=0; j<GlobalV::NLOCAL; j++)
           {
               GlobalV::ofs_running << psi_k.c[i*GlobalV::NLOCAL+j].real() <<
       "+"<<psi_k.c[i*GlobalV::NLOCAL+j].imag()<<"i ";
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
            psi_k[i * GlobalV::NLOCAL + j] = wfc_tmp.c[i * GlobalV::NLOCAL + j];
        }
    }
    /*
            GlobalV::ofs_running<<"psi_k "<<endl;
            for(int i=0; i<GlobalV::NBANDS; i++)
        {
                    for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                GlobalV::ofs_running << psi_k.c[i*GlobalV::NLOCAL+j].real() <<
       "+"<<psi_k.c[i*GlobalV::NLOCAL+j].imag()<<"i ";
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

    delete[] work;
    delete[] ipiv;

    return;
}

#ifdef __MPI

void Evolve_LCAO_Matrix::using_ScaLAPACK_complex(
            const int nband,
            const int nlocal,         
            const std::complex<double>* psi_k_laststep,
            hamilt::Hamilt<double>* p_hamilt,
            std::complex<double>* psi_k,
            double* ekb) const
{
    ModuleBase::TITLE("Evolve_LCAO_Matrix", "using_ScaLAPACK_complex");

    int print_matrix = 0;
    hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
    p_hamilt->matrix(h_mat, s_mat);

    complex<double>* Stmp = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Stmp, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, s_mat.p, 1, Stmp, 1);

    complex<double>* Htmp = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Htmp, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, h_mat.p, 1, Htmp, 1);
    
    complex<double>* U_operator = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(U_operator, this->ParaV->nloc);

// (1)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief compute U_operator
    /// @input Stmp, Htmp, print_matrix
    /// @output U_operator
    compute_U_operator(nband, nlocal, Stmp, Htmp, U_operator, print_matrix);

// (2)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief apply U_operator to the wave function of the previous step for new wave function
    /// @input U_operator, psi_k_laststep
    /// @output psi_k
    U_to_wfc(nband, nlocal, U_operator, psi_k_laststep, psi_k);

    

// (3)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /// @brief normalize psi_k
    /// @input Stmp, psi_not_norm, psi_k, print_matrix
    /// @output psi_k
    norm_wfc(nband, nlocal, Stmp, psi_k, print_matrix);


// (4)->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    /// @brief compute ekb
    /// @input Htmp, psi_k
    /// @output ekb
    compute_ekb(nband, nlocal, Htmp, psi_k, ekb);

    delete[] Stmp;
    delete[] Htmp;
    delete[] U_operator;
    return;
}


void Evolve_LCAO_Matrix::compute_U_operator(
                const int nband,
                const int nlocal,   
                const std::complex<double>* Stmp,
                const std::complex<double>* Htmp,
                std::complex<double>* U_operator,
                const int print_matrix) const
{
    // (1) copy Htmp to Numerator & Denominator
    complex<double>* Numerator = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Numerator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Numerator, 1);

    complex<double>* Denominator = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Denominator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Denominator, 1);

    if (print_matrix)
    {
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " S matrix :" << endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Stmp[i * this->ParaV->ncol + j].real() << "+" << Stmp[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " H matrix :" << endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Numerator[i * this->ParaV->ncol + j].real() << "+" << Numerator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
    }

// ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute Numerator & Denominator by GEADD
    // Numerator = Stmp - i*para * Htmp;     beta1 = - para = -0.25 * INPUT.mdp.md_dt
    // Denominator = Stmp + i*para * Htmp;   beta2 = para = 0.25 * INPUT.mdp.md_dt
    complex<double> alpha = {1.0, 0.0};
    complex<double> beta1 = {0.0, -0.25 * INPUT.mdp.md_dt};
    complex<double> beta2 = {0.0, 0.25 * INPUT.mdp.md_dt};

    ScalapackConnector::geadd(
            'N',
            nlocal,
            nlocal,
            alpha,
            Stmp,
            1,
            1,
            this->ParaV->desc,
            beta1,
            Numerator,
            1,
            1,
            this->ParaV->desc
    );
    ScalapackConnector::geadd(
            'N',
            nlocal,
            nlocal,
            alpha,
            Stmp,
            1,
            1,
            this->ParaV->desc,
            beta2,
            Denominator,
            1,
            1,
            this->ParaV->desc
    );

    if (print_matrix)
    {
        GlobalV::ofs_running << " fenmu:" << endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Denominator[i * this->ParaV->ncol + j].real() << "+" << Denominator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
    }

//->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) Next, invert Denominator
    int* ipiv = new int[this->ParaV->nloc];
    int info = 0;
    // (3.1) compute ipiv
    ScalapackConnector::getrf(
        nlocal,
        nlocal,
        Denominator,
        1,
        1, 
        this->ParaV->desc,
        ipiv, 
        &info
    );
    int lwork = -1;
    int liwotk = -1;
    std::vector<std::complex<double>> work(1, 0);
    std::vector<int> iwork(1, 0);
    // (3.2) compute work
    ScalapackConnector::getri(
        nlocal,
        Denominator,
        1,
        1,
        this->ParaV->desc,
        ipiv,
        work.data(),
        &lwork,
        iwork.data(),
        &liwotk,
        &info
    );
    lwork = work[0].real();
    work.resize(lwork, 0);
    liwotk = iwork[0];
    iwork.resize(liwotk, 0);
    // (3.3) compute inverse matrix of matrix_A
    ScalapackConnector::getri(
        nlocal,
        Denominator,
        1,
        1,
        this->ParaV->desc,
        ipiv,
        work.data(),
        &lwork,
        iwork.data(),
        &liwotk,
        &info
    );
    assert(0 == info);

//->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // (4) U_operator = Denominator * Numerator;
    ScalapackConnector::gemm(
        'N',
        'N',
        nlocal,
        nlocal,
        nlocal,
        1.0,
        Denominator,
        1,
        1,
        this->ParaV->desc,
        Numerator,
        1,
        1,
        this->ParaV->desc,
        0.0,
        U_operator,
        1,
        1,
        this->ParaV->desc
    );

    if (print_matrix)
    {
        GlobalV::ofs_running << " fenmu^-1:" << endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Denominator[i * this->ParaV->ncol + j].real() << "+" << Denominator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " fenzi:" << endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                GlobalV::ofs_running << Numerator[i * this->ParaV->ncol + j].real() << "+" << Numerator[i * this->ParaV->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " U operator:" << endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa, bb;
                aa = U_operator[i * this->ParaV->ncol + j].real();
                bb = U_operator[i * this->ParaV->ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << endl;
        }
    }

// cout << "U_operator Success!!!" <<endl;

    delete[] Numerator;
    delete[] Denominator;
    delete[] ipiv;
}


void Evolve_LCAO_Matrix::U_to_wfc(
                const int nband,
                const int nlocal,   
                const std::complex<double>* U_operator,
                const std::complex<double>* psi_k_laststep,
                std::complex<double>* psi_k) const
{

    ScalapackConnector::gemm(
        'N',
        'N',
        nlocal,
        nband,
        nlocal,
        1.0,
        U_operator,
        1,
        1,
        this->ParaV->desc,
        psi_k_laststep,
        1,
        1,
        this->ParaV->desc_wfc,
        0.0,
        psi_k,
        1,
        1,
        this->ParaV->desc_wfc
    );
}


void Evolve_LCAO_Matrix::norm_wfc(
                const int nband,
                const int nlocal,   
                const std::complex<double>* Stmp,
                std::complex<double>* psi_k,
                const int print_matrix) const
{
    complex<double>* tmp1 = new complex<double>[this->ParaV->nloc_wfc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, this->ParaV->nloc_wfc);
   
    complex<double>* Cij = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Cij, this->ParaV->nloc);
    
    ScalapackConnector::gemm(
            'N',
            'N',
            nlocal,
            nband,
            nlocal,
            1.0,
            Stmp,
            1,
            1,
            this->ParaV->desc,
            psi_k,
            1,
            1,
            this->ParaV->desc_wfc,
            0.0,
            tmp1,
            1,
            1,
            this->ParaV->desc_wfc
    );

    
    ScalapackConnector::gemm(
        'C',
        'N',
        nband,
        nband,
        nlocal,
        1.0,
        psi_k,
        1,
        1,
        this->ParaV->desc_wfc,
        tmp1,
        1,
        1,
        this->ParaV->desc_wfc,
        0.0,
        Cij,
        1,
        1,
        this->ParaV->desc
    );

    int info;
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
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, this->ParaV->nb, this->ParaV->dim1, ipcol);
                    if (igcol >= nband)
                        continue;
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, this->ParaV->nb, this->ParaV->dim0, iprow);
                        if (igrow >= nband)
                            continue;
                        if (igcol == igrow)
                        {
                            Cij[j * naroc[0] + i] = {1.0 / sqrt(Cij[j * naroc[0] + i].real()), 0.0};
                        }
                        else
                        {
                            Cij[j * naroc[0] + i] = {0.0, 0.0};
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow

    // std::cout << "nlocal" << nlocal << std::endl;
    // std::cout << "GlobalV::NLOCAL" << GlobalV::NLOCAL << std::endl;
    BlasConnector::copy(this->ParaV->nloc_wfc, psi_k, 1, tmp1, 1);

    ScalapackConnector::gemm(
        'N',
        'N',
        nlocal,
        nband,
        nband,
        1.0,
        tmp1,
        1,
        1,
        this->ParaV->desc_wfc,
        Cij,
        1,
        1,
        this->ParaV->desc,
        0.0,
        psi_k,
        1,
        1,
        this->ParaV->desc_wfc
    );


    if (print_matrix)
    {
        GlobalV::ofs_running << " Cij:" << endl;
        for (int i = 0; i < this->ParaV->ncol; i++)
        {
            for (int j = 0; j < this->ParaV->nrow; j++)
            {
                GlobalV::ofs_running << Cij[i * this->ParaV->ncol + j].real() << "+" << Cij[i * this->ParaV->ncol_bands + j].imag() << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " psi_k:" << endl;
        for (int i = 0; i < this->ParaV->ncol_bands; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa, bb;
                aa = psi_k[i * this->ParaV->ncol + j].real();
                bb = psi_k[i * this->ParaV->ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << " psi_k nlocal*nlocal:" << endl;
        for (int i = 0; i < this->ParaV->ncol; i++)
        {
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa, bb;
                aa = tmp1[i * this->ParaV->ncol + j].real();
                bb = tmp1[i * this->ParaV->ncol + j].imag();
                if (abs(aa) < 1e-8)
                    aa = 0.0;
                if (abs(bb) < 1e-8)
                    bb = 0.0;
                GlobalV::ofs_running << aa << "+" << bb << "i ";
            }
            GlobalV::ofs_running << endl;
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running << endl;
    }


    delete[] tmp1;
    delete[] Cij;

}


void Evolve_LCAO_Matrix::compute_ekb(
                const int nband,
                const int nlocal,   
                const std::complex<double>* Htmp,
                const std::complex<double>* psi_k,
                double* ekb) const
{

    complex<double>* tmp1 = new complex<double>[this->ParaV->nloc_wfc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, this->ParaV->nloc_wfc);

    complex<double>* Eij = new complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Eij, this->ParaV->nloc);

    ScalapackConnector::gemm(
        'N',
        'N',
        nlocal,
        nband,
        nlocal,
        1.0,
        Htmp,
        1,
        1,
        this->ParaV->desc,
        psi_k,
        1,
        1,
        this->ParaV->desc_wfc,
        0.0,
        tmp1,
        1,
        1,
        this->ParaV->desc_wfc
    );

    ScalapackConnector::gemm(
            'C',
            'N',
            nband,
            nband,
            nlocal,
            1.0,
            psi_k,
            1,
            1,
            this->ParaV->desc_wfc,
            tmp1,
            1,
            1,
            this->ParaV->desc_wfc,
            0.0,
            Eij,
            1,
            1,
            this->ParaV->desc
    );


    if (ELEC_evolve::td_print_eij > 0.0)
    {
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << endl;
        GlobalV::ofs_running << " Eij:" << endl;
        for (int i = 0; i < this->ParaV->ncol; i++)
        {
            for (int j = 0; j < this->ParaV->nrow; j++)
            {
                double aa, bb;
                aa = Eij[i * this->ParaV->ncol + j].real();
                bb = Eij[i * this->ParaV->ncol + j].imag();
                if (abs(aa) < ELEC_evolve::td_print_eij)
                    aa = 0.0;
                if (abs(bb) < ELEC_evolve::td_print_eij)
                    bb = 0.0;
                if (aa > 0.0 || bb > 0.0)
                {
                    GlobalV::ofs_running << i << " " << j << " " << aa << "+" << bb << "i " << endl;
                }
            }
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << endl;
    }


    int info;
    int myid;
    int naroc[2];
    MPI_Comm_rank(this->ParaV->comm_2D, &myid);
    
    double* Eii = new double[nband];
    ModuleBase::GlobalFunc::ZEROS(Eii, nband);
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
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, this->ParaV->nb, this->ParaV->dim1, ipcol);
                    if (igcol >= nband)
                        continue;
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, this->ParaV->nb, this->ParaV->dim0, iprow);
                        if (igrow >= nband)
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
    info = MPI_Allreduce(Eii, ekb, nband, MPI_DOUBLE, MPI_SUM, this->ParaV->comm_2D);
    
    delete[] tmp1;
    delete[] Eij;
    delete[] Eii;
}


#endif
