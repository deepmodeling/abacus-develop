#include "write_HS.h"

#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_io/module_output/filename.h" // use filename_output function


template <typename T>
void ModuleIO::write_hsk(
		const int nspin,
        const int nks, 
        const int nkstot,
        const std::vector<int> &ik2iktot,
        const std::vector<int> &isk,
		hamilt::Hamilt<T>* p_hamilt,
		const Parallel_Orbitals &pv,
		const int istep,
		const ModuleContext::FileSystemLayout& files,
		const ModuleContext::ParallelTopology& parallel,
		const ModuleContext::LogStreams& logs,
		const ModuleContext::BasisInfo& basis,
		const ModuleContext::SolverConfig& solver,
		const ModuleContext::MatrixOutputConfig& output)
{

    *logs.running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        ">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    *logs.running << " |                                            "
        "                        |" << std::endl;
    *logs.running << " | Write Hamiltonian matrix H(k) or overlap matrix S(k) in numerical  |" << std::endl;
    *logs.running << " | atomic orbitals at each k-point.                                   |" << std::endl;
    *logs.running << " |                                            "
        "                        |" << std::endl;
    *logs.running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        ">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    *logs.running << "\n WRITE H(k) OR S(k)" << std::endl;

    for (int ik = 0; ik < nks; ++ik)
    {
        p_hamilt->updateHk(ik);
        bool bit = false; // LiuXh, 2017-03-21
        // if set bit = true, there would be error in soc-multi-core
        // calculation, noted by zhengdy-soc

        hamilt::MatrixBlock<T> h_mat;
        hamilt::MatrixBlock<T> s_mat;

        p_hamilt->matrix(h_mat, s_mat);

        const int out_label=1; // 1: .txt, 2: .dat

        std::string h_fn = ModuleIO::filename_output(files.output_directory,
                "hk","nao",ik,ik2iktot,nspin,nkstot,
                out_label,output.append,basis.gamma_only_local,istep);

        ModuleIO::save_mat(istep,
                h_mat.p,
                basis.nlocal,
                bit,
                output.hs_k.precision,
                1,
                output.append,
                h_fn,
                pv,
                parallel.diagonalization_rank,
                solver);

        // mohan note 2025-06-02
        // for overlap matrix, the two spin channels yield the same matrix
        // so we only need to print matrix from one spin channel.
        const int current_spin = isk[ik];
        if(current_spin == 1)
        {
            continue;
        }

        std::string s_fn = ModuleIO::filename_output(files.output_directory,
                "sk","nao",ik,ik2iktot,nspin,nkstot,
                out_label,output.append,basis.gamma_only_local,istep);

        *logs.running << " The output filename is " << s_fn << std::endl;

        ModuleIO::save_mat(istep,
                s_mat.p,
                basis.nlocal,
                bit,
                output.hs_k.precision,
                1,
                output.append,
                s_fn,    
                pv,
                parallel.diagonalization_rank,
                solver);
    } // end ik
}


// output a square matrix
template <typename T>
void ModuleIO::save_mat(const int istep,
    const T* mat,
    const int dim,
    const bool bit,
    const int precision,
    const bool tri,
    const bool app,
    const std::string& filename,
    const Parallel_2D& pv,
    const int drank,
    const ModuleContext::SolverConfig& solver,
    const bool reduce)
{
    ModuleBase::TITLE("ModuleIO", "save_mat");
    ModuleBase::timer::start("ModuleIO", "save_mat");

    const bool gamma_only = std::is_same<T, double>::value;

    // write .dat file
    if (bit)
    {
// write .dat file with MPI
#ifdef __MPI
        FILE* out_matrix = nullptr;

        if (drank == 0)
        {
            out_matrix = fopen(filename.c_str(), "wb");
            if (out_matrix == nullptr)
            {
                ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
            }
            fwrite(&dim, sizeof(int), 1, out_matrix);
        }

        int ir=0;
        int ic=0;
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
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(solver.ks_solver))
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

            if (reduce) 
            {
                Parallel_Reduce::reduce_all(line, tri ? dim - i : dim);
            }

            if (drank == 0)
            {
                for (int j = (tri ? i : 0); j < dim; ++j)
                {
                    fwrite(&line[tri ? j - i : j], sizeof(T), 1, out_matrix);
                }
            }
            delete[] line;

            MPI_Barrier(DIAG_WORLD);
        }

        if (drank == 0) 
        {
            fclose(out_matrix);
        }
// write .dat file without MPI
#else
        FILE* out_matrix = fopen(filename.c_str(), "wb");
        if (out_matrix == nullptr)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
        }

        fwrite(&dim, sizeof(int), 1, out_matrix);

        for (int i = 0; i < dim; i++)
        {
            for (int j = (tri ? i : 0); j < dim; j++)
            {
                fwrite(&mat[i * dim + j], sizeof(T), 1, out_matrix);
            }
        }
        fclose(out_matrix);
#endif
    } // end writing .dat file
    else // write .txt file
    {
        std::ofstream out_matrix;
        out_matrix << std::scientific << std::setprecision(precision);
#ifdef __MPI
        if (drank == 0)
        {
            if (app && istep > 0) 
            {
                out_matrix.open(filename.c_str(), std::ofstream::app);
            } 
            else 
            {
                out_matrix.open(filename.c_str());
            }
            if (!out_matrix.is_open())
            {
                ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
            }
            out_matrix << "#------------------------------------------------------------------------" << std::endl;
            out_matrix << "# ionic step " << istep+1 << std::endl; // istep starts from 0 
            out_matrix << "# filename " << filename << std::endl;
            out_matrix << "# gamma only " << gamma_only << std::endl;
            out_matrix << "# rows " << dim << std::endl;
            out_matrix << "# columns " << dim << std::endl;
            out_matrix << "#------------------------------------------------------------------------" << std::endl;

        }

        int ir=0;
        int ic=0;
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
                        int iic=0;
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(solver.ks_solver))
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

            if (reduce) 
            {
                Parallel_Reduce::reduce_all(line, tri ? dim - i : dim);
            }

            if (drank == 0)
            {
                out_matrix << "Row " << i+1 << std::endl;
                size_t count = 0;
                for (int j = (tri ? i : 0); j < dim; j++) 
                {
                    out_matrix << " " << line[tri ? j - i : j];
                    ++count;
                    if(count%8==0)
                    {
                        if(j!=dim-1)
                        {
                            out_matrix << std::endl;
                        }
                    }
                }
                out_matrix << std::endl;
            }
            delete[] line;
        }

        if (drank == 0) 
        {
            out_matrix.close();
        }
#else
        if (app)
        {
            out_matrix.open(filename.c_str(), std::ofstream::app);
        }
        else
        {
            out_matrix.open(filename.c_str());
        }
        if (!out_matrix.is_open())
        {
            ModuleBase::WARNING_QUIT("ModuleIO::save_mat", "Cannot open matrix file: " + filename);
        }

        out_matrix << dim;
        out_matrix << std::setprecision(precision);
        for (int i = 0; i < dim; i++)
        {
            for (int j = (tri ? i : 0); j < dim; j++)
            {
                out_matrix << " " << mat[i * dim + j];
            }
            out_matrix << std::endl;
        }
        out_matrix.close();
#endif
    }
    ModuleBase::timer::end("ModuleIO", "save_mat");
    return;
}
