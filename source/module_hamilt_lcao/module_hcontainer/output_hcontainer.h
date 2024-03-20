#ifndef OUTPUT_HCONTAINER_H
#define OUTPUT_HCONTAINER_H

#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

namespace hamilt
{

/**
 * @brief A class to output the HContainer
 */
template <typename T>
class Output_HContainer
{
  public:
    Output_HContainer(hamilt::HContainer<T>* hcontainer,
                      const Parallel_Orbitals* ParaV,
                      const UnitCell& ucell,
                      std::ostream& ofs,
                      double sparse_threshold,
                      int precision);
    // write the matrices of all R vectors to the output stream
    void write();

    /**
     * write the matrix of a single R vector to the output stream
     * rx_in, ry_in, rz_in: the R vector from the input
     */
    void write(int rx_in, int ry_in, int rz_in);

    /**
     * write the matrix of a single R vector to the output stream
     * rx, ry, rz: the R vector from the HContainer
     */
    void write_single_R(int rx, int ry, int rz);

  private:
    hamilt::HContainer<T>* _hcontainer;
    const UnitCell& _ucell;
    const Parallel_Orbitals* _ParaV;
    std::ostream& _ofs;
    double _sparse_threshold;
    int _precision;
};

} // namespace ModuleIO

#endif // OUTPUT_HCONTAINER_H