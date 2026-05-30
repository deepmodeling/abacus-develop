#ifndef ABFS_CONSTRUCT_PCA_H
#define ABFS_CONSTRUCT_PCA_H

#include "../../source_basis/module_ao/ORB_read.h"
#include "source_cell/unitcell.h"
#include <vector>
#include <RI/global/Tensor.h>

//	training data: lcaos[i] * lcaos[j]
//	old basis:     abfs
//	new basis:     to be constructed
//	( all lcaos and abfs on same atom )

namespace ABFs_Construct
{
namespace PCA
{
extern std::vector<std::vector<std::pair<std::vector<double>, RI::Tensor<double>>>> cal_PCA(
		const UnitCell& ucell,
    const LCAO_Orbitals& orb,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& lcaos,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs, // abfs must be orthonormal
    const double kmesh_times);

template <typename T>
void tensor_syev(char jobz, char uplo, RI::Tensor<T>& a, RI::Global_Func::To_Real_t<T>* w, int& info);
} // namespace PCA
} // namespace ABFs_Construct

#endif // ABFS_CONSTRUCT_PCA_H
