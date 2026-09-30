#ifndef EXX_SOURCE_IO_H
#define EXX_SOURCE_IO_H

#include <iosfwd>
#include <string>

class K_Vectors;
namespace ModuleBase
{
class matrix;
}

namespace ModuleIO
{
// A versioned companion to WAVEFUNC*.dat, containing the frozen SCF ensemble.
// Streams keep serialization independent of process ownership and file paths.
void write_exx_source(std::ostream& out,
                      const K_Vectors& points,
                      const ModuleBase::matrix& weights,
                      int nspin,
                      int nx,
                      int ny,
                      int nz,
                      const std::string& configuration);

void read_exx_source(std::istream& in,
                     K_Vectors& points,
                     ModuleBase::matrix& weights,
                     int nspin,
                     int nx,
                     int ny,
                     int nz,
                     const std::string& configuration);
}
#endif
