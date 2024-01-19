#include "write_HS.h"

#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/formatter_physfmt.h"

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

// output a square matrix
template <typename T>
void ModuleIO::save_mat(const int istep,
    const T* mat,
    const int dim,
    const bool bit,
    const int precision,
    const bool tri,
    const bool app,
    const std::string label,
    const std::string& file_name,
    const Parallel_2D& pv,
    const int drank)
{
    ModuleBase::TITLE("ModuleIO", "save_mat");
    ModuleBase::timer::tick("ModuleIO", "save_mat");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Dimension of " + label + " : ", dim);

    std::stringstream ss;

    formatter::PhysicalFmt physfmt;
    physfmt.adjust_formatter_flexible(precision, -1, true);
    if (bit)ss << GlobalV::global_out_dir << file_name + "-" + label + "-bit";
    else
    {
        if (app || istep < 0)
            ss << GlobalV::global_out_dir << file_name + "-" + label;
        else
            ss << GlobalV::global_out_dir << istep << "_" << file_name + "-" + label;
    }
    if (bit)
    {
#ifdef __MPI
        FILE* g = nullptr;

        if (drank == 0)
        {
            g = fopen(ss.str().c_str(), "wb");
            fwrite(&dim, sizeof(int), 1, g);
        }

        int ir, ic;
        for (int i = 0; i < dim; ++i)
        {
            T* line = new T[tri ? dim - i : dim];
            ModuleBase::GlobalFunc::ZEROS(line, tri ? dim - i : dim);

            ir = pv.global2local_row(i);
            if (ir >= 0)
            {
                // data collection
                for (int j = (tri ? i : 0); j < dim; ++j)
                {
                    ic = pv.global2local_col(j);
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
                        line[tri ? j - i : j] = mat[iic];
                    }
                }
            }

            Parallel_Reduce::reduce_all(line, tri ? dim - i : dim);

            if (drank == 0)
            {
                for (int j = (tri ? i : 0); j < dim; ++j)
                {
                    fwrite(&line[tri ? j - i : j], sizeof(T), 1, g);
                }
            }
            delete[] line;

            MPI_Barrier(DIAG_WORLD);
        }

        if (drank == 0)
            fclose(g);
#else
        FILE* g = fopen(ss.str().c_str(), "wb");

        fwrite(&dim, sizeof(int), 1, g);

        for (int i = 0; i < dim; i++)
        {
            for (int j = (tri ? i : 0); j < dim; j++)
            {
                fwrite(&mat[i * dim + j], sizeof(T), 1, g);
            }
        }
        fclose(g);
#endif
    } // end bit
    else
    {
        std::ofstream g;
#ifdef __MPI
        if (drank == 0)
        {
            if (app)
                g.open(ss.str().c_str(), std::ofstream::app);
            else
                g.open(ss.str().c_str());
            g << dim;
        }

        int ir, ic;
        for (int i = 0; i < dim; i++)
        {
            T* line = new T[tri ? dim - i : dim];
            ModuleBase::GlobalFunc::ZEROS(line, tri ? dim - i : dim);

            ir = pv.global2local_row(i);
            if (ir >= 0)
            {
                // data collection
                for (int j = (tri ? i : 0); j < dim; ++j)
                {
                    ic = pv.global2local_col(j);
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
                        line[tri ? j - i : j] = mat[iic];
                    }
                }
            }

            Parallel_Reduce::reduce_all(line, tri ? dim - i : dim);

            if (drank == 0)
            {
                for (int j = (tri ? i : 0); j < dim; j++) g << " " << physfmt.get_p_formatter()->format(line[tri ? j - i : j]);
                g << std::endl;
            }
            delete[] line;
        }

        if (drank == 0) // Peize Lin delete ; at 2020.01.31
            g.close();
#else
        if (app)
            std::ofstream g(ss.str().c_str(), std::ofstream::app);
        else
            std::ofstream g(ss.str().c_str());

        g << dim;

        for (int i = 0; i < dim; i++)
        {
            for (int j = (tri ? i : 0); j < dim; j++)
            {
                g << " " << physfmt.get_p_formatter()->format(mat[i * dim + j]);
            }
            g << std::endl;
        }
        g.close();
#endif
    }
    ModuleBase::timer::tick("ModuleIO", "save_mat");
    return;
}


/* comment out this function for not used
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

                    ir = lm.ParaV->global2local_row(i);
                    if (ir >= 0)
                    {
                        // for(int j=i; j<GlobalV::NLOCAL; j++)
                        for (int j = 0; j < GlobalV::NLOCAL; j++)
                        {
                            ic = lm.ParaV->global2local_col(j);
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

                    // Parallel_Reduce::reduce_all(lineH,GlobalV::NLOCAL-i);
                    // Parallel_Reduce::reduce_all(lineS,GlobalV::NLOCAL-i);
                    if (GlobalV::NSPIN != 4)
                    {
                        Parallel_Reduce::reduce_all(lineH, GlobalV::NLOCAL);
                        Parallel_Reduce::reduce_all(lineS, GlobalV::NLOCAL);
                    }
                    else
                    {
                        Parallel_Reduce::reduce_all(lineH_soc, GlobalV::NLOCAL);
                        Parallel_Reduce::reduce_all(lineS_soc, GlobalV::NLOCAL);
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

                                // if(GlobalV::DRANK==0);
                                // {
                                //     g1.close();
                                //     g2.close();
                                // }

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
*/
