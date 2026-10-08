#ifndef EXX_SOURCE_IO_H
#define EXX_SOURCE_IO_H

#include <iosfwd>
#include <string>
#include <vector>

class K_Vectors;
namespace ModuleBase
{
class matrix;
}

namespace ModuleIO
{
// Metadata already present in the existing binary PW wavefunction format.
struct ExxSourceHeader
{
    int ik;
    int nks;
    int npw;
    int nbands;
    double kvec_c[3];
    double weight;
    double ecutwfc;
    double lat0;
    double tpiba;
};

ExxSourceHeader read_exx_source_header(std::istream& in);

// Read the existing, k-weighted occupation column without changing its format.
void read_exx_source_occupations(std::istream& in,
                                 const std::vector<ExxSourceHeader>& headers,
                                 int nspin,
                                 ModuleBase::matrix& weights);

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
