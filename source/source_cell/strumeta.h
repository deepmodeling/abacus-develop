#ifndef STRUMETA_H
#define STRUMETA_H

#include <string>
#include <vector>

/**
 * @brief STRU fields retained only to reproduce an MD restart STRU file.
 *
 * Atom labels, masses, and atom counts belong to MDCell because they are
 * physical topology data. This type deliberately contains only input/output
 * information that is not used by MD integration or force evaluation.
 */
struct StruSpecies
{
    std::string pseudo_file;
    std::string pseudo_type;
    std::string orbital_file;
    double start_mag = 0.0;
};

struct StruMeta
{
    std::vector<StruSpecies> species;
    std::string descriptor_file;
};

#endif
