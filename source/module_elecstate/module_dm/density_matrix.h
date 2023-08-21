#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

#include <string>

#include "module_cell/klist.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
// #include "module_io/csr_reader.h"

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
    DensityMatrix(const K_Vectors* _kv, const Parallel_Orbitals* _paraV);

    // initialize density matrix DMR from UnitCell
    void init_DMR(Grid_Driver* GridD_in, const UnitCell* ucell);

    // initialize density matrix DMR from another HContainer
    void init_DMR(const hamilt::HContainer<T>& _DMR_in);

    // set _DMK element directly
    void set_DMK(const int ik, const int i, const int j, const T value);

    // read *.dmk into density matrix dm(k) with all k-points
    void set_DMK_files(const std::string directory);

    // write density matrix dm(k) into *.dmk
    void write_DMK(const std::string directory, const int ik);

    // read *.dmk into density matrix dm(k)
    void read_DMK(const std::string directory, const int ik);

    // write density matrix dm(k) into *.dmk with all k-points
    void output_DMK(const std::string directory);

    // get a matrix element of density matrix dm(k)
    T get_DMK(const int ik, const int i, const int j) const;

    // get size information of DMK
    int get_DMK_nks() const;
    int get_DMK_nrow() const;
    int get_DMK_ncol() const;

    // get DMR pointer
    hamilt::HContainer<T>* get_DMR_pointer();

    // calculate DMR from DMK
    void cal_DMR();

  private:
    // HContainer for density matrix in real space
    hamilt::HContainer<T>* _DMR = nullptr;

    // density matrix in k space
    std::vector<ModuleBase::ComplexMatrix> _DMK;

    // kvector item
    const K_Vectors* _kv;

    // pointer of Parallel_Orbitals, which is used to get atom-pair information
    const Parallel_Orbitals* _paraV = nullptr;
};

} // namespace elecstate

#endif
