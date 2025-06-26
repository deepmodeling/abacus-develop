#include "write_HS.h"

#include "module_parameter/parameter.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_pw/hamilt_pwdft/global.h"
#include "module_io/filename.h" // use filename_output function


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

    // print out .dat file
	if (bit)
	{
#ifdef __MPI
        FILE* g = nullptr;

        if (drank == 0)
        {
            g = fopen(filename.c_str(), "wb");
            fwrite(&dim, sizeof(int), 1, g);
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
                    fwrite(&line[tri ? j - i : j], sizeof(T), 1, g);
                }
            }
            delete[] line;

            MPI_Barrier(DIAG_WORLD);
        }

		if (drank == 0) 
		{
			fclose(g);
		}
#else
        FILE* g = fopen(filename.c_str(), "wb");

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
    } // end .dat file
    else // .txt file
    {
        std::ofstream g;
        g << std::setprecision(precision);
#ifdef __MPI
        if (drank == 0)
        {
			if (app && istep > 0) 
			{
				g.open(filename.c_str(), std::ofstream::app);
			} 
			else 
			{
				g.open(filename.c_str());
			}
			g << dim;
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
				for (int j = (tri ? i : 0); j < dim; j++) 
				{
					g << " " << line[tri ? j - i : j];
				}
				g << std::endl;
            }
            delete[] line;
        }

		if (drank == 0) 
		{ // Peize Lin delete ; at 2020.01.31
			g.close();
		}
#else
		if (app)
		{
			std::ofstream g(filename.c_str(), std::ofstream::app);
		}
		else
		{
			std::ofstream g(filename.c_str());
		}

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
