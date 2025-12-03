#ifndef READ_HCONTAINER_H
#define READ_HCONTAINER_H

#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_cell/unitcell.h"

namespace hamilt
{

/**
 * @brief A class to read the HContainer
 */
template <typename T>
class Read_HContainer
{
  public:
    Read_HContainer(
        hamilt::HContainer<T>* hcontainer, 
        const std::string& filename,
        const int nlocal,
        const UnitCell* ucell
    );
    // read the matrices of all R vectors to the read stream
    void read();

    /**
     * read the matrix of a single R vector to the output stream
     * rx_in, ry_in, rz_in: the R vector from the input
     */
    void read(int rx_in, int ry_in, int rz_in);

    /**
     * read the matrix of a single R vector to the output stream
     * rx, ry, rz: the R vector from the HContainer
     */
    void read_single_R(int rx, int ry, int rz);

  private:
    hamilt::HContainer<T>* _hcontainer;
    std::string _filename;
    int _nlocal;
    const UnitCell* _ucell;
};

} // namespace hamilt

#endif // OUTPUT_HCONTAINER_H