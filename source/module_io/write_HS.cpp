#include "write_HS.h"

#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void ModuleIO::saving_HS(const int istep,
                         const double* Hloc,
                         const double* Sloc,
                         const bool bit,
                         const int& out_mat_hs,
                         const std::string& file_name,
                         const Parallel_Orbitals& pv,
                         bool tri)
{

    if (out_mat_hs == 1)
    {
        if (tri)
        {
            save_HS_triangle(istep, Hloc, Sloc, bit, file_name, pv);
        }
        else
        {
            save_HS_complete(istep, Hloc, Sloc, bit, file_name, pv);
        }
    }
    else if (out_mat_hs == 2)
    {
        if (tri)
        {
            save_HS_triangle(istep, Hloc, Sloc, bit, file_name, pv);
        }
        else
        {
            save_HS_complete(istep, Hloc, Sloc, bit, file_name, pv);
        }
    }
    else if (out_mat_hs == 3)
    {
        // please call individually
    }
    else if (out_mat_hs == 0)
    {
        // do nothing.
    }
    else
    {
        ModuleBase::WARNING("Diago_LCAO_Matrix", "unrecorganized out_mat_hs value.");
    }
    return;
}

/*
void ModuleIO::save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H,
        const double *nzval_H, const double *nzval_S, bool bit)
{
    ModuleBase::TITLE("ModuleIO","save_HS_ccf");

    if(GlobalV::DRANK!=0)return;

    std::stringstream ssh;
    std::stringstream sss;

    if(bit)
    {
        ssh << GlobalV::global_out_dir << "H_bit.ccf";
        sss << GlobalV::global_out_dir << "S_bit.ccf";
    }
    else
    {
        // mohan update 2021-02-10
        ssh << GlobalV::global_out_dir << "H" << ELEC_scf::iter << "_" << iter+1 << ".ccf";
        sss << GlobalV::global_out_dir << "S" << ELEC_scf::iter << "_" << iter+1 << ".ccf";
    }

    if(bit)
    {
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&GlobalV::NLOCAL,sizeof(int),1,g1);
        fwrite(&Hnnz,sizeof(int),1,g1);
        fwrite(&GlobalV::NLOCAL,sizeof(int),1,g2);
        fwrite(&Hnnz,sizeof(int),1,g2);

        fclose(g1);
        fclose(g2);
    }


    if(!bit)
    {
        std::ofstream g1(ssh.str().c_str());
        std::ofstream g2(sss.str().c_str());

        g1 << GlobalV::NLOCAL << " " << Hnnz << std::endl;
        g2 << GlobalV::NLOCAL << " " << Hnnz << std::endl;

        for(int i=0; i<GlobalV::NLOCAL+1; ++i)
        {
            g1 << colptr_H[i] << " ";
            g2 << colptr_H[i] << " ";
        }
        g1 << std::endl;
        g2 << std::endl;

        for(int i=0; i<Hnnz; ++i)
        {
            g1 << rowind_H[i] << " ";
            g2 << rowind_H[i] << " ";
        }
        g1 << std::endl;
        g2 << std::endl;

        for(int i=0; i<Hnnz; ++i)
        {
            g1 << nzval_H[i] << " ";
            g2 << nzval_S[i] << " ";
        }
        g1 << std::endl;
        g2 << std::endl;

        g1.close();
        g2.close();
    }
    return;
}
*/

// mohan add 2010/3/20, output H and S matrix, convinence for diagonalization
// test or save the middle information for next start.
void ModuleIO::save_HS_triangle(const int istep,
                                const double* H,
                                const double* S,
                                const bool bit,
                                const std::string& file_name,
                                const Parallel_Orbitals& pv)
{
    ModuleBase::TITLE("ModuleIO", "save_HS_bit");
    ModuleBase::timer::tick("ModuleIO", "save_HS_bit");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Dimension of H and S", GlobalV::NLOCAL);

    std::stringstream ssh;
    std::stringstream sss;

    if (bit)
    {
        ssh << GlobalV::global_out_dir << file_name + "-H-bit";
        sss << GlobalV::global_out_dir << file_name + "-S-bit";
    }
    else
    {
        if (GlobalV::out_app_flag)
        {
            ssh << GlobalV::global_out_dir << file_name + "-H";
            sss << GlobalV::global_out_dir << file_name + "-S";
        }
        else
        {
            ssh << GlobalV::global_out_dir << istep << "_" << file_name + "-H";
            sss << GlobalV::global_out_dir << istep << "_" << file_name + "-S";
        }
    }
    if (bit)
    {
#ifdef __MPI
        FILE *g1 = nullptr;
        FILE *g2 = nullptr;

        if (GlobalV::DRANK == 0)
        {
            g1 = fopen(ssh.str().c_str(), "wb");
            g2 = fopen(sss.str().c_str(), "wb");
            fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g1);
            fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g2);
        }

        int ir, ic;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            double *lineH = new double[GlobalV::NLOCAL - i];
            double *lineS = new double[GlobalV::NLOCAL - i];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL - i);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL - i);

            ir = pv.global2local_row[i];
            if (ir >= 0)
            {
                // data collection
                for (int j = i; j < GlobalV::NLOCAL; j++)
                {
                    ic = pv.global2local_col[j];
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        // lineH[j-i] = H[ir*pv.ncol+ic];
                        // lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j - i] = H[iic];
                        lineS[j - i] = S[iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH, GlobalV::NLOCAL - i);
            Parallel_Reduce::reduce_double_all(lineS, GlobalV::NLOCAL - i);

            if (GlobalV::DRANK == 0)
            {
                for (int j = i; j < GlobalV::NLOCAL; j++)
                {
                    fwrite(&lineH[j - i], sizeof(double), 1, g1);
                    fwrite(&lineS[j - i], sizeof(double), 1, g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (GlobalV::DRANK == 0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(), "wb");
        FILE *g2 = fopen(sss.str().c_str(), "wb");

        fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g1);
        fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g2);

        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = i; j < GlobalV::NLOCAL; j++)
            {
                fwrite(&H[i * GlobalV::NLOCAL + j], sizeof(double), 1, g1);
                fwrite(&S[i * GlobalV::NLOCAL + j], sizeof(double), 1, g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } // end bit
    else
    {
#ifdef __MPI
        std::ofstream g1;
        std::ofstream g2;

        if (GlobalV::DRANK == 0)
        {
            if (GlobalV::out_app_flag)
            {
                g1.open(ssh.str().c_str(), std::ofstream::app);
                g2.open(sss.str().c_str(), std::ofstream::app);
            }
            else
            {
                g1.open(ssh.str().c_str());
                g2.open(sss.str().c_str());
            }
            g1 << GlobalV::NLOCAL;
            g2 << GlobalV::NLOCAL;
        }

        int ir, ic;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            double *lineH = new double[GlobalV::NLOCAL - i];
            double *lineS = new double[GlobalV::NLOCAL - i];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL - i);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL - i);

            ir = pv.global2local_row[i];
            if (ir >= 0)
            {
                // data collection
                for (int j = i; j < GlobalV::NLOCAL; j++)
                {
                    ic = pv.global2local_col[j];
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        // lineH[j-i] = H[ir*pv.ncol+ic];
                        // lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j - i] = H[iic];
                        lineS[j - i] = S[iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH, GlobalV::NLOCAL - i);
            Parallel_Reduce::reduce_double_all(lineS, GlobalV::NLOCAL - i);

            if (GlobalV::DRANK == 0)
            {
                for (int j = i; j < GlobalV::NLOCAL; j++)
                {
                    g1 << " " << lineH[j - i];
                    g2 << " " << lineS[j - i];
                }
                g1 << std::endl;
                g2 << std::endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        // if (GlobalV::DRANK==0);
        if (GlobalV::DRANK == 0) // Peize Lin delete ; at 2020.01.31
        {
            g1.close();
            g2.close();
        }
#else
        if (GlobalV::out_app_flag)
        {
            std::ofstream g1(ssh.str().c_str(), ofstream::app);
            std::ofstream g2(sss.str().c_str(), ofstream::app);
        }
        else
        {
            std::ofstream g1(ssh.str().c_str());
            std::ofstream g2(sss.str().c_str());
        }

        g1 << GlobalV::NLOCAL;
        g2 << GlobalV::NLOCAL;

        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = i; j < GlobalV::NLOCAL; j++)
            {
                g1 << " " << H[i * GlobalV::NLOCAL + j];
                g2 << " " << S[i * GlobalV::NLOCAL + j];
            }
            g1 << std::endl;
            g2 << std::endl;
        }
        g1.close();
        g2.close();
#endif
    }

    ModuleBase::timer::tick("ModuleIO", "save_HS_bit");
    return;
}

// mohan add 2010/3/20, output H and S matrix, convinence for diagonalization
// test or save the middle information for next start.
void ModuleIO::save_HS_complete(const int istep,
                                const double* H,
                                const double* S,
                                const bool bit,
                                const std::string& file_name,
                                const Parallel_Orbitals& pv)
{
    ModuleBase::TITLE("ModuleIO", "save_HS_bit");
    ModuleBase::timer::tick("ModuleIO", "save_HS_bit");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Dimension of H and S", GlobalV::NLOCAL);

    std::stringstream ssh;
    std::stringstream sss;

    if (bit)
    {
        ssh << GlobalV::global_out_dir << file_name + "-H-bit";
        sss << GlobalV::global_out_dir << file_name + "-S-bit";
    }
    else
    {
        if (GlobalV::out_app_flag)
        {
            ssh << GlobalV::global_out_dir << file_name + "-H";
            sss << GlobalV::global_out_dir << file_name + "-S";
        }
        else
        {
            ssh << GlobalV::global_out_dir << istep << "_" << file_name + "-H";
            sss << GlobalV::global_out_dir << istep << "_" << file_name + "-S";
        }
    }
    if (bit)
    {
#ifdef __MPI
        FILE *g1 = nullptr;
        FILE *g2 = nullptr;

        if (GlobalV::DRANK == 0)
        {
            g1 = fopen(ssh.str().c_str(), "wb");
            g2 = fopen(sss.str().c_str(), "wb");
            fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g1);
            fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g2);
        }

        int ir, ic;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            double *lineH = new double[GlobalV::NLOCAL];
            double *lineS = new double[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL);

            ir = pv.global2local_row[i];
            if (ir >= 0)
            {
                // data collection
                for (int j = 0; j < GlobalV::NLOCAL; j++)
                {
                    ic = pv.global2local_col[j];
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        // lineH[j-i] = H[ir*pv.ncol+ic];
                        // lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j] = H[iic];
                        lineS[j] = S[iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH, GlobalV::NLOCAL);
            Parallel_Reduce::reduce_double_all(lineS, GlobalV::NLOCAL);

            if (GlobalV::DRANK == 0)
            {
                for (int j = 0; j < GlobalV::NLOCAL; j++)
                {
                    fwrite(&lineH[j], sizeof(double), 1, g1);
                    fwrite(&lineS[j], sizeof(double), 1, g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (GlobalV::DRANK == 0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(), "wb");
        FILE *g2 = fopen(sss.str().c_str(), "wb");

        fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g1);
        fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g2);

        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = 0; j < GlobalV::NLOCAL; j++)
            {
                fwrite(&H[i * GlobalV::NLOCAL + j], sizeof(double), 1, g1);
                fwrite(&S[i * GlobalV::NLOCAL + j], sizeof(double), 1, g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } // end bit
    else
    {
#ifdef __MPI
        std::ofstream g1;
        std::ofstream g2;

        if (GlobalV::DRANK == 0)
        {
            if (GlobalV::out_app_flag)
            {
                g1.open(ssh.str().c_str(), std::ofstream::app);
                g2.open(sss.str().c_str(), std::ofstream::app);
            }
            else
            {
                g1.open(ssh.str().c_str());
                g2.open(sss.str().c_str());
            }
            g1 << GlobalV::NLOCAL;
            g2 << GlobalV::NLOCAL;
        }

        int ir, ic;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            double *lineH = new double[GlobalV::NLOCAL];
            double *lineS = new double[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL);

            ir = pv.global2local_row[i];
            if (ir >= 0)
            {
                // data collection
                for (int j = 0; j < GlobalV::NLOCAL; j++)
                {
                    ic = pv.global2local_col[j];
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        // lineH[j-i] = H[ir*pv.ncol+ic];
                        // lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j] = H[iic];
                        lineS[j] = S[iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH, GlobalV::NLOCAL);
            Parallel_Reduce::reduce_double_all(lineS, GlobalV::NLOCAL);

            if (GlobalV::DRANK == 0)
            {
                for (int j = 0; j < GlobalV::NLOCAL; j++)
                {
                    g1 << " " << lineH[j];
                    g2 << " " << lineS[j];
                }
                g1 << std::endl;
                g2 << std::endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        // if (GlobalV::DRANK==0);
        if (GlobalV::DRANK == 0) // Peize Lin delete ; at 2020.01.31
        {
            g1.close();
            g2.close();
        }
#else
        if (GlobalV::out_app_flag)
        {
            std::ofstream g1(ssh.str().c_str(), ofstream::app);
            std::ofstream g2(sss.str().c_str(), ofstream::app);
        }
        else
        {
            std::ofstream g1(ssh.str().c_str());
            std::ofstream g2(sss.str().c_str());
        }

        g1 << GlobalV::NLOCAL;
        g2 << GlobalV::NLOCAL;

        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = 0; j < GlobalV::NLOCAL; j++)
            {
                g1 << " " << H[i * GlobalV::NLOCAL + j];
                g2 << " " << S[i * GlobalV::NLOCAL + j];
            }
            g1 << std::endl;
            g2 << std::endl;
        }
        g1.close();
        g2.close();
#endif
    }

    ModuleBase::timer::tick("ModuleIO", "save_HS_bit");
    return;
}

// LiuXh, 2017-03-21
void ModuleIO::saving_HS(const int istep,
                         std::complex<double>* Hloc,
                         std::complex<double>* Sloc,
                         const bool bit,
                         const int& out_mat_hs,
                         const std::string& file_name,
                         const Parallel_Orbitals& pv,
                         bool tri)
{
    if (out_mat_hs == 1)
    {
        if (tri)
        {
            save_HS_complex_triangle(istep, Hloc, Sloc, bit, file_name, pv);
        }
        else
        {
            save_HS_complex_complete(istep, Hloc, Sloc, bit, file_name, pv);
        }
    }
    else if (out_mat_hs == 0)
    {
        // do nothing.
    }
    else
    {
        ModuleBase::WARNING("Diago_LCAO_Matrix", "unrecorganized out_mat_hs value.");
    }
    return;
}

// LiuXh, 2017-03-21
void ModuleIO::save_HS_complex_triangle(const int istep,
                                        std::complex<double>* H,
                                        std::complex<double>* S,
                                        const bool bit,
                                        const std::string& file_name,
                                        const Parallel_Orbitals& pv)
{
    ModuleBase::TITLE("ModuleIO", "save_HS_bit");
    ModuleBase::timer::tick("ModuleIO", "save_HS_bit");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Dimension of H and S", GlobalV::NLOCAL);

    std::stringstream ssh;
    std::stringstream sss;

    if (bit)
    {
        ssh << GlobalV::global_out_dir << file_name + "-H-bit";
        sss << GlobalV::global_out_dir << file_name + "-S-bit";
    }
    else
    {
        if (GlobalV::out_app_flag)
        {
            ssh << GlobalV::global_out_dir << file_name + "-H";
            sss << GlobalV::global_out_dir << file_name + "-S";
        }
        else
        {
            ssh << GlobalV::global_out_dir << istep << "_" << file_name + "-H";
            sss << GlobalV::global_out_dir << istep << "_" << file_name + "-S";
        }
    }

    if (bit)
    {
#ifdef __MPI
        FILE *g1 = nullptr;
        FILE *g2 = nullptr;

        if (GlobalV::DRANK == 0)
        {
            g1 = fopen(ssh.str().c_str(), "wb");
            g2 = fopen(sss.str().c_str(), "wb");
            fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g1);
            fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g2);
        }

        int ir, ic;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            std::complex<double> *lineH = new std::complex<double>[GlobalV::NLOCAL - i];
            std::complex<double> *lineS = new std::complex<double>[GlobalV::NLOCAL - i];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL - i);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL - i);

            ir = pv.global2local_row[i];
            if (ir >= 0)
            {
                // data collection
                for (int j = i; j < GlobalV::NLOCAL; j++)
                {
                    ic = pv.global2local_col[j];
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        // lineH[j-i] = H[ir*pv.ncol+ic];
                        // lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j - i] = H[iic];
                        lineS[j - i] = S[iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_complex_double_pool(lineH, GlobalV::NLOCAL - i);
            Parallel_Reduce::reduce_complex_double_pool(lineS, GlobalV::NLOCAL - i);

            if (GlobalV::DRANK == 0)
            {
                for (int j = i; j < GlobalV::NLOCAL; j++)
                {
                    fwrite(&lineH[j - i], sizeof(std::complex<double>), 1, g1);
                    fwrite(&lineS[j - i], sizeof(std::complex<double>), 1, g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (GlobalV::DRANK == 0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(), "wb");
        FILE *g2 = fopen(sss.str().c_str(), "wb");

        fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g1);
        fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g2);

        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = i; j < GlobalV::NLOCAL; j++)
            {
                fwrite(&H[i * GlobalV::NLOCAL + j], sizeof(std::complex<double>), 1, g1);
                fwrite(&S[i * GlobalV::NLOCAL + j], sizeof(std::complex<double>), 1, g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } // end bit
    else
    {
#ifdef __MPI
        std::ofstream g1;
        std::ofstream g2;

        if (GlobalV::DRANK == 0)
        {
            if (GlobalV::out_app_flag)
            {
                g1.open(ssh.str().c_str(), std::ofstream::app);
                g2.open(sss.str().c_str(), std::ofstream::app);
            }
            else
            {
                g1.open(ssh.str().c_str());
                g2.open(sss.str().c_str());
            }
            g1 << GlobalV::NLOCAL;
            g2 << GlobalV::NLOCAL;
        }

        int ir, ic;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            std::complex<double> *lineH = new std::complex<double>[GlobalV::NLOCAL - i];
            std::complex<double> *lineS = new std::complex<double>[GlobalV::NLOCAL - i];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL - i);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL - i);

            ir = pv.global2local_row[i];
            if (ir >= 0)
            {
                // data collection
                for (int j = i; j < GlobalV::NLOCAL; j++)
                {
                    ic = pv.global2local_col[j];
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        // lineH[j-i] = H[ir*pv.ncol+ic];
                        // lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j - i] = H[iic];
                        lineS[j - i] = S[iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_complex_double_pool(lineH, GlobalV::NLOCAL - i);
            Parallel_Reduce::reduce_complex_double_pool(lineS, GlobalV::NLOCAL - i);

            if (GlobalV::DRANK == 0)
            {
                for (int j = i; j < GlobalV::NLOCAL; j++)
                {
                    g1 << " " << lineH[j - i];
                    g2 << " " << lineS[j - i];
                }
                g1 << std::endl;
                g2 << std::endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        // if (GlobalV::DRANK==0);
        if (GlobalV::DRANK == 0) // Peize Lin delete ; at 2020.01.31
        {
            g1.close();
            g2.close();
        }
#else

        if (GlobalV::out_app_flag)
        {
            std::ofstream g1(ssh.str().c_str(), ofstream::app);
            std::ofstream g2(sss.str().c_str(), ofstream::app);
        }
        else
        {
            std::ofstream g1(ssh.str().c_str());
            std::ofstream g2(sss.str().c_str());
        }

        g1 << GlobalV::NLOCAL;
        g2 << GlobalV::NLOCAL;

        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = i; j < GlobalV::NLOCAL; j++)
            {
                g1 << " " << H[i * GlobalV::NLOCAL + j];
                g2 << " " << S[i * GlobalV::NLOCAL + j];
            }
            g1 << std::endl;
            g2 << std::endl;
        }
        g1.close();
        g2.close();
#endif
    }

    ModuleBase::timer::tick("ModuleIO", "save_HS_bit");
    return;
}

// LiuXh, 2017-03-21
void ModuleIO::save_HS_complex_complete(const int istep,
                                        std::complex<double>* H,
                                        std::complex<double>* S,
                                        const bool bit,
                                        const std::string& file_name,
                                        const Parallel_Orbitals& pv)
{
    ModuleBase::TITLE("ModuleIO", "save_HS_bit");
    ModuleBase::timer::tick("ModuleIO", "save_HS_bit");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Dimension of H and S", GlobalV::NLOCAL);

    std::stringstream ssh;
    std::stringstream sss;

    if (bit)
    {
        ssh << GlobalV::global_out_dir << file_name + "-H-bit";
        sss << GlobalV::global_out_dir << file_name + "-S-bit";
    }
    else
    {
        if (GlobalV::out_app_flag)
        {
            ssh << GlobalV::global_out_dir << file_name + "-H";
            sss << GlobalV::global_out_dir << file_name + "-S";
        }
        else
        {
            ssh << GlobalV::global_out_dir << istep << "_" << file_name + "-H";
            sss << GlobalV::global_out_dir << istep << "_" << file_name + "-S";
        }
    }

    if (bit)
    {
#ifdef __MPI
        FILE *g1 = nullptr;
        FILE *g2 = nullptr;

        if (GlobalV::DRANK == 0)
        {
            g1 = fopen(ssh.str().c_str(), "wb");
            g2 = fopen(sss.str().c_str(), "wb");
            fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g1);
            fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g2);
        }

        int ir, ic;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            std::complex<double> *lineH = new std::complex<double>[GlobalV::NLOCAL];
            std::complex<double> *lineS = new std::complex<double>[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL);

            ir = pv.global2local_row[i];
            if (ir >= 0)
            {
                // data collection
                for (int j = 0; j < GlobalV::NLOCAL; j++)
                {
                    ic = pv.global2local_col[j];
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        // lineH[j-i] = H[ir*pv.ncol+ic];
                        // lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j] = H[iic];
                        lineS[j] = S[iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_complex_double_pool(lineH, GlobalV::NLOCAL);
            Parallel_Reduce::reduce_complex_double_pool(lineS, GlobalV::NLOCAL);

            if (GlobalV::DRANK == 0)
            {
                for (int j = 0; j < GlobalV::NLOCAL; j++)
                {
                    fwrite(&lineH[j], sizeof(std::complex<double>), 1, g1);
                    fwrite(&lineS[j], sizeof(std::complex<double>), 1, g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (GlobalV::DRANK == 0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(), "wb");
        FILE *g2 = fopen(sss.str().c_str(), "wb");

        fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g1);
        fwrite(&GlobalV::NLOCAL, sizeof(int), 1, g2);

        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = 0; j < GlobalV::NLOCAL; j++)
            {
                fwrite(&H[i * GlobalV::NLOCAL + j], sizeof(std::complex<double>), 1, g1);
                fwrite(&S[i * GlobalV::NLOCAL + j], sizeof(std::complex<double>), 1, g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } // end bit
    else
    {
#ifdef __MPI
        std::ofstream g1;
        std::ofstream g2;

        if (GlobalV::DRANK == 0)
        {
            if (GlobalV::out_app_flag)
            {
                g1.open(ssh.str().c_str(), std::ofstream::app);
                g2.open(sss.str().c_str(), std::ofstream::app);
            }
            else
            {
                g1.open(ssh.str().c_str());
                g2.open(sss.str().c_str());
            }
            g1 << GlobalV::NLOCAL;
            g2 << GlobalV::NLOCAL;
        }

        int ir, ic;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            std::complex<double> *lineH = new std::complex<double>[GlobalV::NLOCAL];
            std::complex<double> *lineS = new std::complex<double>[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL);
            ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL);

            ir = pv.global2local_row[i];
            if (ir >= 0)
            {
                // data collection
                for (int j = 0; j < GlobalV::NLOCAL; j++)
                {
                    ic = pv.global2local_col[j];
                    if (ic >= 0)
                    {
                        int iic;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * pv.nrow;
                        }
                        else
                        {
                            iic = ir * pv.ncol + ic;
                        }
                        // lineH[j-i] = H[ir*pv.ncol+ic];
                        // lineS[j-i] = S[ir*pv.ncol+ic];
                        lineH[j] = H[iic];
                        lineS[j] = S[iic];
                    }
                }
            }
            else
            {
                // do nothing
            }

            Parallel_Reduce::reduce_complex_double_pool(lineH, GlobalV::NLOCAL);
            Parallel_Reduce::reduce_complex_double_pool(lineS, GlobalV::NLOCAL);

            if (GlobalV::DRANK == 0)
            {
                for (int j = 0; j < GlobalV::NLOCAL; j++)
                {
                    g1 << " " << lineH[j];
                    g2 << " " << lineS[j];
                }
                g1 << std::endl;
                g2 << std::endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        // if (GlobalV::DRANK==0);
        if (GlobalV::DRANK == 0) // Peize Lin delete ; at 2020.01.31
        {
            g1.close();
            g2.close();
        }
#else
        if (GlobalV::out_app_flag)
        {
            std::ofstream g1(ssh.str().c_str(), ofstream::app);
            std::ofstream g2(sss.str().c_str(), ofstream::app);
        }
        else
        {
            std::ofstream g1(ssh.str().c_str());
            std::ofstream g2(sss.str().c_str());
        }

        g1 << GlobalV::NLOCAL;
        g2 << GlobalV::NLOCAL;

        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = 0; j < GlobalV::NLOCAL; j++)
            {
                g1 << " " << H[i * GlobalV::NLOCAL + j];
                g2 << " " << S[i * GlobalV::NLOCAL + j];
            }
            g1 << std::endl;
            g2 << std::endl;
        }
        g1.close();
        g2.close();
#endif
    }

    ModuleBase::timer::tick("ModuleIO", "save_HS_bit");
    return;
}

// void ModuleIO::save_HSR_tr(const int Rx, const int Ry, const int Rz, const double *H, const double *S)
void ModuleIO::save_HSR_tr(const int current_spin, LCAO_Matrix &lm)
// void ModuleIO::save_HSR_tr(void)
{
    ModuleBase::TITLE("ModuleIO", "save_HSR_tr");
    ModuleBase::timer::tick("ModuleIO", "save_HSR_tr");

    std::stringstream ssh;
    std::stringstream sss;

    ssh << GlobalV::global_out_dir << "data-HR-tr_SPIN" << current_spin;
    sss << GlobalV::global_out_dir << "data-SR-tr_SPIN" << current_spin;
    // ssh << GlobalV::global_out_dir << "data-HR-tr_SPIN";
    // sss << GlobalV::global_out_dir << "data-SR-tr_SPIN";

#ifdef __MPI
    std::ofstream g1;
    std::ofstream g2;

    if (GlobalV::DRANK == 0)
    {
        g1.open(ssh.str().c_str());
        g2.open(sss.str().c_str());
        g1 << "Matrix Dimension of H(R): " << GlobalV::NLOCAL << std::endl;
        g2 << "Matrix Dimension of S(R): " << GlobalV::NLOCAL << std::endl;
    }

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    // std::cout<<"R_x: "<<R_x<<std::endl;
    // std::cout<<"R_y: "<<R_y<<std::endl;
    // std::cout<<"R_z: "<<R_z<<std::endl;

    double R_minX = GlobalC::GridD.getD_minX();
    double R_minY = GlobalC::GridD.getD_minY();
    double R_minZ = GlobalC::GridD.getD_minZ();

    // int dRx, dRy, dRz;

    for (int ix = 0; ix < R_x; ix++)
    {
        int dRx = ix + R_minX;
        for (int iy = 0; iy < R_y; iy++)
        {
            int dRy = iy + R_minY;
            for (int iz = 0; iz < R_z; iz++)
            {
                int dRz = iz + R_minZ;
                // std::cout<<"dRx: "<<dRx<<std::endl;
                // std::cout<<"dRy: "<<dRy<<std::endl;
                // std::cout<<"dRz: "<<dRz<<std::endl;
                int ir, ic;
                for (int i = 0; i < GlobalV::NLOCAL; i++)
                {
                    // double* lineH = new double[GlobalV::NLOCAL-i];
                    // double* lineS = new double[GlobalV::NLOCAL-i];
                    double *lineH = nullptr;
                    double *lineS = nullptr;
                    std::complex<double> *lineH_soc = nullptr;
                    std::complex<double> *lineS_soc = nullptr;
                    if (GlobalV::NSPIN != 4)
                    {
                        lineH = new double[GlobalV::NLOCAL];
                        lineS = new double[GlobalV::NLOCAL];
                        ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL);
                        ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL);
                    }
                    else
                    {
                        lineH_soc = new std::complex<double>[GlobalV::NLOCAL];
                        lineS_soc = new std::complex<double>[GlobalV::NLOCAL];
                        ModuleBase::GlobalFunc::ZEROS(lineH_soc, GlobalV::NLOCAL);
                        ModuleBase::GlobalFunc::ZEROS(lineS_soc, GlobalV::NLOCAL);
                    }
                    // ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL-i);
                    // ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL-i);
                    // ModuleBase::GlobalFunc::ZEROS(lineH, GlobalV::NLOCAL);
                    // ModuleBase::GlobalFunc::ZEROS(lineS, GlobalV::NLOCAL);

                    ir = lm.ParaV->global2local_row[i];
                    if (ir >= 0)
                    {
                        // for(int j=i; j<GlobalV::NLOCAL; j++)
                        for (int j = 0; j < GlobalV::NLOCAL; j++)
                        {
                            ic = lm.ParaV->global2local_col[j];
                            if (ic >= 0)
                            {
                                // lineH[j-i] = H[ir*lm.ParaV->ncol+ic];
                                // lineS[j-i] = S[ir*lm.ParaV->ncol+ic];
                                int iic;
                                if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                                {
                                    iic = ir + ic * lm.ParaV->nrow;
                                }
                                else
                                {
                                    iic = ir * lm.ParaV->ncol + ic;
                                }
                                if (GlobalV::NSPIN != 4)
                                {
                                    lineH[j] = lm.HR_tr[ix][iy][iz][iic];
                                    lineS[j] = lm.SlocR_tr[ix][iy][iz][iic];
                                }
                                else
                                {
                                    lineH_soc[j] = lm.HR_tr_soc[ix][iy][iz][iic];
                                    lineS_soc[j] = lm.SlocR_tr_soc[ix][iy][iz][iic];
                                }
                            }
                        }
                    }
                    else
                    {
                        // do nothing
                    }

                    // Parallel_Reduce::reduce_double_all(lineH,GlobalV::NLOCAL-i);
                    // Parallel_Reduce::reduce_double_all(lineS,GlobalV::NLOCAL-i);
                    if (GlobalV::NSPIN != 4)
                    {
                        Parallel_Reduce::reduce_double_all(lineH, GlobalV::NLOCAL);
                        Parallel_Reduce::reduce_double_all(lineS, GlobalV::NLOCAL);
                    }
                    else
                    {
                        Parallel_Reduce::reduce_complex_double_all(lineH_soc, GlobalV::NLOCAL);
                        Parallel_Reduce::reduce_complex_double_all(lineS_soc, GlobalV::NLOCAL);
                    }

                    if (GlobalV::DRANK == 0)
                    {
                        // for(int j=i; j<GlobalV::NLOCAL; j++)
                        for (int j = 0; j < GlobalV::NLOCAL; j++)
                        {
                            if (i == 0 && j == 0)
                            {
                                g1 << dRx << " " << dRy << " " << dRz
                                   << "    //R std::vector(R2 - R1,unit: lattice vector)" << std::endl;
                                g2 << dRx << " " << dRy << " " << dRz
                                   << "    //R std::vector(R2 - R1,unit: lattice vector)" << std::endl;
                            }
                            // g1 << " " << lineH[j-i];
                            // g2 << " " << lineS[j-i];
                            if (GlobalV::NSPIN != 4)
                            {
                                if (std::abs(lineH[j]) < 1.0e-12)
                                    lineH[j] = 0.0;
                                if (std::abs(lineS[j]) < 1.0e-12)
                                    lineS[j] = 0.0;
                                g1 << " " << lineH[j];
                                g2 << " " << lineS[j];
                            }
                            else
                            {
                                if (std::abs(lineH_soc[j].real()) < 1.0e-12)
                                    lineH_soc[j] = std::complex<double>(0.0, lineH_soc[j].imag());
                                if (std::abs(lineH_soc[j].imag()) < 1.0e-12)
                                    lineH_soc[j] = std::complex<double>(lineH_soc[j].real(), 0.0);
                                if (std::abs(lineS_soc[j].real()) < 1.0e-12)
                                    lineS_soc[j] = std::complex<double>(0.0, lineS_soc[j].imag());
                                if (std::abs(lineS_soc[j].imag()) < 1.0e-12)
                                    lineS_soc[j] = std::complex<double>(lineS_soc[j].real(), 0.0);
                                g1 << " " << lineH_soc[j];
                                g2 << " " << lineS_soc[j];
                            }
                        }
                        g1 << std::endl;
                        g2 << std::endl;
                    }
                    if (GlobalV::NSPIN != 4)
                    {
                        delete[] lineH;
                        delete[] lineS;
                    }
                    else
                    {
                        delete[] lineH_soc;
                        delete[] lineS_soc;
                    }
                }
                /*
                                if(GlobalV::DRANK==0);
                                {
                                    g1.close();
                                    g2.close();
                                }
                */
            }
        }
    }
    // if(GlobalV::DRANK==0);
    if (GlobalV::DRANK == 0) // Peize Lin delete ; at 2020.01.31
    {
        g1.close();
        g2.close();
    }

#else
    std::ofstream g1(ssh.str().c_str());
    std::ofstream g2(sss.str().c_str());

    g1 << GlobalV::NLOCAL;
    g2 << GlobalV::NLOCAL;

    for (int i = 0; i < GlobalV::NLOCAL; i++)
    {
        for (int j = i; j < GlobalV::NLOCAL; j++)
        {
            // not correct for sequential version; need change
            // g1 << " " << H[i*GlobalV::NLOCAL+j];
            // g2 << " " << S[i*GlobalV::NLOCAL+j];
        }
    }
    g1.close();
    g2.close();
#endif

    ModuleBase::timer::tick("ModuleIO", "save_HSR_tr");
    return;
}
