#ifndef ATOMIC_RADIALS_H_
#define ATOMIC_RADIALS_H_

#include "module_basis/module_nao/radial_set.h"

//! The radial part of numerical atomic orbitals of a single element
/*!
 *  This class represents the radial part of all numerical atomic orbitals
 *  of a single element as read from an orbital file.
 *
 *  @see RadialSet
 *
 *  Usage:
 *
 *      int element_index = 1;
 *      std::ofstream ofs_log("/path/to/log/file");
 *      std::string orb_file = "/path/to/orbital/file";
 *
 *      AtomicRadials O_radials;
 *      O_radials.build(orb_file, element_index, ofs_log, GlobalV::MY_RANK);
 *
 *                                                                          */
class AtomicRadials : public RadialSet
{
  public:
    AtomicRadials() {}
    AtomicRadials(const AtomicRadials& other) : RadialSet(other), orb_ecut_(other.orb_ecut_) {}

    AtomicRadials& operator=(const AtomicRadials& rhs);
    AtomicRadials* clone() const { return new AtomicRadials(*this); } // covariant return type

    ~AtomicRadials() {} // ~RadialSet() is called automatically

    //! Build the class from an orbital file
    void build(const std::string& file,          //!< orbital file name
               const int itype = 0,              //!< element index in calculation
               std::ofstream* ptr_log = nullptr, //!< output file stream for logging
               const int rank = 0                //!< MPI rank
    );

    void build(RadialSet* const other, const int itype, const double rcut);

    //! Get the energy cutoff as given by the orbital file
    double orb_ecut() const { return orb_ecut_; }

    /**
     * @brief static version of read_abacus_orb. A delete-new operation may cause the memory leak, 
     * it is better to use std::vector to replace the raw pointer.
     * 
     * @param ifs [in] ifstream from the orbital file, via `std::ifstream ifs(forb);`
     * @param elem [out] element symbol
     * @param ecut [out] planewave energy cutoff
     * @param nr [out] number of radial grid points
     * @param dr [out] radial grid spacing
     * @param nzeta [out] number of zeta functions for each angular momentum
     * @param radials [out] radial orbitals
     * @param rank [in] MPI rank
     */
    static void read_abacus_orb(std::ifstream& ifs,
                                std::string& elem,
                                double& ecut,
                                int& nr,
                                double& dr,
                                std::vector<int>& nzeta,
                                std::vector<std::vector<double>>& radials,
                                const int rank = 0);

  private:
    double orb_ecut_; //!< energy cutoff as given by the orbital file

    //! Read the orbital file in the ABACUS format
    void read_abacus_orb(std::ifstream& ifs,               //!< input file stream from orbital file
                         std::ofstream* ptr_log = nullptr, //!< output file stream for logging
                         const int rank = 0                //!< MPI rank
    );

};

#endif
