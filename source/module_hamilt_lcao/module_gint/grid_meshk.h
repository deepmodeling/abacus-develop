#ifndef GRID_MESHK_H
#define GRID_MESHK_H
#include "module_base/global_function.h"
#include "module_base/global_variable.h"

class Grid_MeshK
{
	public:
		Grid_MeshK();
		~Grid_MeshK();
		// from 1D index to unitcell.
		std::vector<int> ucell_index2x;
		std::vector<int> ucell_index2y;
		std::vector<int> ucell_index2z;

		// the unitcell parameters.
		std::vector<int> max_ucell_para;
		std::vector<int> min_ucell_para;
		std::vector<int> num_ucell_para;

		// calculate the index of unitcell.
        int cal_Rindex(const int& u1, const int& u2, const int& u3)const;

        /// move operator for the next ESolver to directly use its infomation
        Grid_MeshK& operator=(Grid_MeshK&& rhs) = default;

      protected:
		// the max and the min unitcell.
		int maxu1;
		int maxu2;
		int maxu3;

		int minu1;
		int minu2;
		int minu3;

		// the number of unitcells.
		int nu1;
		int nu2;
		int nu3;
		int nutot;

		// calculate the extended unitcell.
		void cal_extended_cell(const int &dxe, const int &dye, const int &dze,
								const int& nbx, const int& nby, const int& nbz);
		// initialize the unitcell parameters.
		void init_ucell_para(void);
};

#endif
