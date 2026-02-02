#include "write_HS.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_io/module_output/filename.h" // use filename_output function


template <typename T>
void ModuleIO::write_hsk(
        const std::string &global_out_dir,
        const int nspin,
        const int nks, 
        const int nkstot,
		const std::vector<int> &ik2iktot,
        const std::vector<int> &isk,
		hamilt::Hamilt<T>* p_hamilt,
	    const Parallel_Orbitals &pv,
        const bool gamma_only,
        const bool out_app_flag,
        const int istep,
        std::ofstream &ofs_running)	
{

	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
		">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	ofs_running << " |                                            "
		"                        |" << std::endl;
	ofs_running << " | Write Hamiltonian matrix H(k) or overlap matrix S(k) in numerical  |" << std::endl; 
	ofs_running << " | atomic orbitals at each k-point.                                   |" << std::endl; 
	ofs_running << " |                                            "
		"                        |" << std::endl;
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
		">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	ofs_running << "\n WRITE H(k) OR S(k)" << std::endl;

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

		std::string h_fn = ModuleIO::filename_output(global_out_dir,
				"hk","nao",ik,ik2iktot,nspin,nkstot,
				out_label,out_app_flag,gamma_only,istep);

		ModuleIO::save_mat(istep,
				h_mat.p,
				PARAM.globalv.nlocal,
				bit,
				PARAM.inp.out_mat_hs[1],
				1,
				out_app_flag,
				h_fn,
				pv,
				GlobalV::DRANK);

        // mohan note 2025-06-02
        // for overlap matrix, the two spin channels yield the same matrix
        // so we only need to print matrix from one spin channel.
		const int current_spin = isk[ik];
		if(current_spin == 1)
		{
			continue;
		}

		std::string s_fn = ModuleIO::filename_output(global_out_dir,
				"sk","nao",ik,ik2iktot,nspin,nkstot,
				out_label,out_app_flag,gamma_only,istep);

		ofs_running << " The output filename is " << s_fn << std::endl;

		ModuleIO::save_mat(istep,
				s_mat.p,
				PARAM.globalv.nlocal,
				bit,
				PARAM.inp.out_mat_hs[1],
				1,
				out_app_flag,
			    s_fn,	
				pv,
				GlobalV::DRANK);
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
    const bool reduce)
{
    ModuleBase::TITLE("ModuleIO", "save_mat");
    ModuleBase::timer::tick("ModuleIO", "save_mat");

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
			std::ofstream out_matrix(filename.c_str(), std::ofstream::app);
		}
		else
		{
			std::ofstream out_matrix(filename.c_str());
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
    ModuleBase::timer::tick("ModuleIO", "save_mat");
    return;
}
