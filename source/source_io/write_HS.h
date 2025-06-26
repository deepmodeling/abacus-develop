#ifndef WRITE_HS_H
#define WRITE_HS_H

#include <string>
#include <vector>

//#include "source_base/global_function.h"
//#include "source_base/global_variable.h"
#include "source_basis/module_ao/parallel_orbitals.h" // use Parallel_Orbitals


// mohan add this file 2010-09-10
namespace ModuleIO
{
	template<typename T>
		void write_hsk(
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
				std::ofstream &ofs_running);	

    /// @brief save a square matrix, such as H(k) and S(k)
    /// @param[in] istep : the step of the calculation
    /// @param[in] mat : the local matrix
    /// @param[in] bit : true for binary, false for decimal
    /// @param[in] tri : true for upper triangle, false for full matrix
    /// @param[in] app : true for append, false for overwrite
    /// @param[in] file_name : the name of the output file
    /// @param[in] pv : the 2d-block parallelization information
    /// @param[in] drank : the rank of the current process
    template<typename T>
    void save_mat(const int istep,
        const T* mat,
        const int dim,
        const bool bit,
        const int precision,
        const bool tri,
        const bool app,
        const std::string& file_name,
        const Parallel_2D& pv,
        const int drank,
        const bool reduce = true);

}
#include "write_HS.hpp"
#endif
