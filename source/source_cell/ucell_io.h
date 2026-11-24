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
     * @brief Skips the UnitCell block that follows the DMK header.
     *
     * The function consumes the latName, lat0, three lattice vectors, atom
     * labels, per-type atom counts, the "Direct" marker, and every coordinate
     * line in order. On success the stream is positioned at the first line of
     * the trailing matrix body.
     *
     * @param ifs The input file stream positioned at the start of the block.
     * @return true when the full block is consumed without parse errors;
     *         false when any mandatory line cannot be read, or when the atom
     *         count line yields no positive integer (leaving the stream at
     *         the coordinates, which would otherwise be misinterpreted as
     *         the matrix body by the caller).
     */
    static bool read_ucell(std::ifstream& ifs);
};

} // namespace ModuleIO

#endif // UCELL_IO_H
