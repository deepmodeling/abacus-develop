#include "LCAO_diago.h"

#include "../module_base/timer.h"
#include "../src_io/wf_local.h"
#include "../src_io/write_HS.h"
#include "../src_pdiag/pdiag_double.h"
#include "../src_pw/global.h"

Diago_LCAO_Matrix::Diago_LCAO_Matrix(LCAO_Matrix *lm) : LM(lm)
{
}
Diago_LCAO_Matrix::~Diago_LCAO_Matrix()
{
}

//wfc_k has been replaced by psi, this part needs rewriting
/*
void Diago_LCAO_Matrix::solve_complex_matrix(const int &ik, Local_Orbital_wfc &lowf)
{
    ModuleBase::TITLE("Diago_LCAO_Matrix", "solve_complex_matrix");
    time_t time_start = time(NULL);

    if (GlobalV::KS_SOLVER == "lapack")
    {
        this->using_LAPACK_complex(ik, lowf.wfc_k_grid[ik], lowf.wfc_k[ik]);
    }
    else
    {
#ifdef __MPI
        this->using_HPSEPS_complex(ik, lowf);
#else
        ModuleBase::WARNING_QUIT("Diago_LCAO_Matrix::solve_complex_matrix", "only lapack is available!");
#endif
    }

    time_t time_end = time(NULL);

    ModuleBase::GlobalFunc::OUT_TIME("diago(std::complex)", time_start, time_end);

    return;
}
*/


//wfc_k has been replaced by psi, this part needs rewriting
/*
void Diago_LCAO_Matrix::using_LAPACK_complex(const int &ik,
                                             std::complex<double> **wfc_k_grid,
                                             ModuleBase::ComplexMatrix &wfc_k) const
{
    ModuleBase::TITLE("Diago_LCAO_Matrix", "using_LAPACK_complex");

    assert(GlobalV::NPROC == 1);

    ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    ModuleBase::ComplexMatrix Stmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    for (int i = 0; i < GlobalV::NLOCAL; i++)
    {
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            Htmp(i, j) = this->LM->Hloc2[i * GlobalV::NLOCAL + j];
            Stmp(i, j) = this->LM->Sloc2[i * GlobalV::NLOCAL + j];
        }
    }

    //----------------------------
    // keep this for tests
    //    out.printcm_norm("Lapack_H", Htmp, 1.0e-5);
    //    out.printcm_norm("Lapack_S", Stmp, 1.0e-5);
    //----------------------------

    double *en = new double[GlobalV::NLOCAL];
    ModuleBase::GlobalFunc::ZEROS(en, GlobalV::NLOCAL);

    ModuleBase::ComplexMatrix hvec(GlobalV::NLOCAL, GlobalV::NBANDS);
    GlobalC::hm.diagH_LAPACK(GlobalV::NLOCAL, GlobalV::NBANDS, Htmp, Stmp, GlobalV::NLOCAL, en, hvec);

    wfc_k.create(GlobalV::NBANDS, GlobalV::NLOCAL);
    if (GlobalV::NSPIN != 4)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
            {
                wfc_k_grid[ib][iw] = hvec(iw, ib);
                wfc_k.c[ib * GlobalV::NLOCAL + iw] = wfc_k_grid[ib][iw];
            }
        }
    }
    else
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            for (int iw = 0; iw < GlobalV::NLOCAL / GlobalV::NPOL; iw++)
            {
                wfc_k_grid[ib][iw] = hvec(iw * GlobalV::NPOL, ib);
                wfc_k_grid[ib][iw + GlobalV::NLOCAL / GlobalV::NPOL] = hvec(iw * GlobalV::NPOL + 1, ib);
            }
        }
    }

    // energy for k-point ik
    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        GlobalC::wf.ekb[ik][ib] = en[ib];
    }

    if (Pdiag_Double::out_wfc_lcao)
    {
        std::stringstream ss;
        ss << GlobalV::global_out_dir << "LOWF_K_" << ik + 1 << ".dat";
        //WF_Local::write_lowf_complex(ss.str(), wfc_k_grid, ik, GlobalC::wf.ekb, GlobalC::wf.wg);
    }

    delete[] en;
    return;
}
*/

void Diago_LCAO_Matrix::using_LAPACK(const int &ik, Local_Orbital_wfc &lowf, double* ekb_ik) const
{
    ModuleBase::TITLE("Diago_LCAO_Matrix", "using_LAPACK");
    assert(GlobalV::NLOCAL > 0);

    // save H and S matrix to disk.
    //	bool bit = false;
    bool bit = true; // zhengdy-soc
    HS_Matrix::saving_HS(this->LM->Hloc.data(),
                         this->LM->Sloc.data(),
                         bit,
                         this->out_mat_hs,
                         "data-" + std::to_string(ik),
                         *lowf.ParaV);

    ModuleBase::matrix Htmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    ModuleBase::matrix Stmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    for (int i = 0; i < GlobalV::NLOCAL; i++)
    {
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            Htmp(i, j) = this->LM->Hloc[i * GlobalV::NLOCAL + j];
            Stmp(i, j) = this->LM->Sloc[i * GlobalV::NLOCAL + j];
        }
    }

    // @@@@@@@
    // test
    // @@@@@@@
    //    out.printrm("Lapack_H", Htmp);
    //    out.printrm("Lapack_S", Stmp);

    int itype = 1;
    int lwork = 3 * GlobalV::NLOCAL - 1; // tmp
    double *w = new double[GlobalV::NLOCAL];
    double *work = new double[lwork];
    ModuleBase::GlobalFunc::ZEROS(w, GlobalV::NLOCAL);
    ModuleBase::GlobalFunc::ZEROS(work, lwork);
    int info;

    clock_t clock_start, clock_end;
    clock_start = std::clock();
    LapackConnector::dsygv(itype,
                           'V',
                           'U',
                           GlobalV::NLOCAL,
                           Htmp,
                           GlobalV::NLOCAL,
                           Stmp,
                           GlobalV::NLOCAL,
                           w,
                           work,
                           lwork,
                           &info);
    clock_end = std::clock();
    double duration = (double)(clock_end - clock_start) / CLOCKS_PER_SEC;

    GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(20);
    //	GlobalV::ofs_running << " clock_start = " << clock_start << std::endl;
    //	GlobalV::ofs_running << " clock_end = " << clock_end << std::endl;
    GlobalV::ofs_running << " Time using dsygv in LAPACK (seconds) is " << duration << std::endl;
    GlobalV::ofs_running << std::resetiosflags(ios::fixed) << std::setprecision(10);

    for (int i = 0; i < GlobalV::NBANDS; i++)
    {
        // eigenvalues
        ekb_ik[i] = w[i];
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            lowf.wfc_k_grid[ik][i][j] = Htmp(j, i);
        }
    }

    // @@@@@@@
    // test
    // @@@@@@@
    /*
    std::cout << "\n Lapack, wfc after diago:" << std::endl;
    for(int i=0; i<GlobalV::NBANDS; i++)
    {
        std::cout << " Eigenvalue from LAPACK : " << std::setw(5) << std::setw(12) << GlobalC::wf.ekb[ik][i] <<
    std::endl; std::cout << " Eigenfunctions" << std::endl; for(int j=0; j<GlobalV::NLOCAL; j++)
        {
            std::cout << std::setw(12) << wfc[i][j];
        }
        std::cout << std::endl;
    }
    //exit(0);
    */

    delete[] w;
    delete[] work;
    return;
}
