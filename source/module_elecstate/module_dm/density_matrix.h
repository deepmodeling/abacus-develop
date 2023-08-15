#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_cell/klist.h"
#include <string>
//#include "module_io/csr_reader.h"

namespace elecstate
{
/**
 * class: DensityMatrix
 */
template <typename T>
class DensityMatrix
{
  public:
  // Destructor of class DensityMatrix
  ~DensityMatrix();

  // Constructor of class DensityMatrix
  DensityMatrix(const int _nlocal, const K_Vectors* _kv, const Parallel_Orbitals* _paraV = nullptr);

  // write density matrix dm(k) into *.dmk
	void write_dmk(const std::string directory,const int ik);

	// read *.dmk into density matrix dm(k)
	void read_dmk(const std::string directory,const int ik);

  // write density matrix dm(k) into *.dmk with all k-points
  void write_all_dmk(const std::string directory);

  // read *.dmk into density matrix dm(k) with all k-points
  void read_all_dmk(const std::string directory);

  // get a matrix element of density matrix dm(k)
  T get_dmK(const int ik, const int i, const int j) const;

  private:
    // HContainer for density matrix in real space
    hamilt::HContainer<T>* _dmR = nullptr;

    // density matrix in k space
    std::vector<ModuleBase::ComplexMatrix> _dmK;

    // kvector item
    const K_Vectors* _kv;

    // number of total orbitals
    int _nlocal = 0;

    //pointer of Parallel_Orbitals, which is used to get atom-pair information
    const Parallel_Orbitals* _paraV = nullptr;
};

} // namespace elecstate

#endif
