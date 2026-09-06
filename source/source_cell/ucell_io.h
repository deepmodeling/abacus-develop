#ifndef UCELL_IO_H
#define UCELL_IO_H

#include "source_cell/unitcell.h"

#include <fstream>

namespace ModuleIO {

/**
 * @brief A class for unit cell I/O operations
 * 
 * This class provides methods to write and read unit cell information
 * to/from files, particularly for DMK files.
 */
class UcellIO {
public:
    /**
     * @brief Writes the unit cell information to a file.
     *
     * @param ofs The output file stream.
     * @param ucell A pointer to the UnitCell object.
     */
    static void write_ucell(std::ofstream& ofs, const UnitCell* ucell);

    /**
     * @brief Reads the unit cell information lines in a file.
     *
     * @param ifs The input file stream.
     */
    static void read_ucell(std::ifstream& ifs);
};

} // namespace ModuleIO

#endif // UCELL_IO_H
