#ifndef PW_HYBRID_SOURCE_H
#define PW_HYBRID_SOURCE_H

#include <iosfwd>
#include <string>
#include <vector>

namespace ModuleBase
{
class matrix;
}

namespace pw
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

/**
 * @brief Locate the EXX source file written by the latest output step.
 * @details The fixed file name is preferred (modes 0 and 1). When it does not
 *          exist, the highest g{ionic step} variant is selected, and within
 *          that step the highest e{electronic iteration} variant. This is the
 *          rule used for mode-2 output directories, where every ionic step
 *          gets its own file.
 * @param directories Candidate directories, tried in order.
 * @param fixed_name File name without the directory, e.g. "wfk1s1_pw.dat".
 * @param tail Fixed suffix that follows the optional step tokens, e.g.
 *        "_pw.dat" or ".txt".
 * @param allow_electronic_suffix Whether g{step}e{iter} names are accepted.
 * @param ionic_step When non-null, receives the parsed one-based ionic step,
 *        or -1 for a fixed (step-less) file name.
 * @param electronic_step When non-null, receives the parsed one-based
 *        electronic iteration, or -1.
 * @return Full path of the selected file.
 * @throws std::runtime_error When no candidate exists.
 */
std::string resolve_source_file(const std::vector<std::string>& directories,
                                const std::string& fixed_name,
                                const std::string& tail,
                                const bool allow_electronic_suffix,
                                int* ionic_step = nullptr,
                                int* electronic_step = nullptr);

}
#endif
