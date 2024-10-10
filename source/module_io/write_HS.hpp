#include "write_HS.h"

#include "module_parameter/parameter.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
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
        ssh << PARAM.globalv.global_out_dir << "H_bit.ccf";
        sss << PARAM.globalv.global_out_dir << "S_bit.ccf";
    }
    else
    {
        // mohan update 2021-02-10
        ssh << PARAM.globalv.global_out_dir << "H" << ELEC_scf::iter << "_" << iter+1 << ".ccf";
        sss << PARAM.globalv.global_out_dir << "S" << ELEC_scf::iter << "_" << iter+1 << ".ccf";
    }

    if(bit)
    {
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&PARAM.globalv.nlocal,sizeof(int),1,g1);
        fwrite(&Hnnz,sizeof(int),1,g1);
        fwrite(&PARAM.globalv.nlocal,sizeof(int),1,g2);
        fwrite(&Hnnz,sizeof(int),1,g2);

        fclose(g1);
        fclose(g2);
    }


    if(!bit)
    {
        std::ofstream g1(ssh.str().c_str());
        std::ofstream g2(sss.str().c_str());

        g1 << PARAM.globalv.nlocal << " " << Hnnz << std::endl;
        g2 << PARAM.globalv.nlocal << " " << Hnnz << std::endl;

        for(int i=0; i<PARAM.globalv.nlocal+1; ++i)
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
    const int drank,
    const bool reduce)
{
    ModuleBase::TITLE("ModuleIO", "save_mat");
    ModuleBase::timer::tick("ModuleIO", "save_mat");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Dimension of " + label + " : ", dim);

    std::stringstream ss;
    if (bit) {ss << PARAM.globalv.global_out_dir << file_name + "-" + label + "-bit";
    } else
    {
        if (app || istep < 0) {
            ss << PARAM.globalv.global_out_dir << file_name + "-" + label;
        } else {
            ss << PARAM.globalv.global_out_dir << istep << "_" << file_name + "-" + label;
}
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
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
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

            if (reduce) Parallel_Reduce::reduce_all(line, tri ? dim - i : dim);

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

        if (drank == 0) {
            fclose(g);
}
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
        g << std::setprecision(precision);
#ifdef __MPI
        if (drank == 0)
        {
            if (app && istep > 0) {
                g.open(ss.str().c_str(), std::ofstream::app);
            } else {
                g.open(ss.str().c_str());
}
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
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
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

            if (reduce) Parallel_Reduce::reduce_all(line, tri ? dim - i : dim);

            if (drank == 0)
            {
                for (int j = (tri ? i : 0); j < dim; j++) { g << " " << line[tri ? j - i : j];
}
                g << std::endl;
            }
            delete[] line;
        }

        if (drank == 0) { // Peize Lin delete ; at 2020.01.31
            g.close();
}
#else
        if (app)
            std::ofstream g(ss.str().c_str(), std::ofstream::app);
        else
            std::ofstream g(ss.str().c_str());

        g << dim;
        g << std::setprecision(precision);
        for (int i = 0; i < dim; i++)
        {
            for (int j = (tri ? i : 0); j < dim; j++)
            {
                g << " " << mat[i * dim + j];
            }
            g << std::endl;
        }
        g.close();
#endif
    }
    ModuleBase::timer::tick("ModuleIO", "save_mat");
    return;
}