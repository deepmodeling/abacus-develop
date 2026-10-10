#ifndef PW_HYBRID_SOURCE_H
#define PW_HYBRID_SOURCE_H

#include <iosfwd>
#include <vector>

namespace ModuleBase
{
class matrix;
}

namespace ModuleESolver
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

}
#endif
