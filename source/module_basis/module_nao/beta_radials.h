#ifndef BETA_RADIALS_H_
#define BETA_RADIALS_H_

#include "module_basis/module_nao/radial_set.h"

//! The radial part of beta functions of a single element
/*!
 *  This class represents the radial part of all Kleinman-Bylander beta
 *  functions of a single element as read from a pseudopotential file.
 *
 *  @see RadialSet
 *
 *  Usage:
 *
 *      int element_index = 1;
 *      std::ofstream ofs_log("/path/to/log/file");
 *      std::string upf_file = "/path/to/pseudopotential/file";
 *
 *      BetaRadials O_beta;
 *      O_beta.build(orb_file, element_index, ofs_log, GlobalV::MY_RANK);
 *
 *                                                                          */
class BetaRadials : public RadialSet
{
  public:
    BetaRadials() {}
    BetaRadials(const BetaRadials& other) : RadialSet(other) {} //!< deep copy

    using RadialSet::operator=;
    BetaRadials* clone() const { return new BetaRadials(*this); } // covariant return type

    ~BetaRadials() {}

    //! Build the class from a pseudopotential file
    void build(const std::string& file,          //!< pseudopotential file name
               const int itype = 0,              //!< element index in calculation
               std::ofstream* ptr_log = nullptr, //!< output file stream for logging
               const int rank = 0                //!< MPI rank
    );

  private:
    //! Read beta projectors from a pseudopotential file of UPF 1.0.0 format
    void read_beta_upf100(std::ifstream& ifs,               //!< input file stream from orbital file
                          std::ofstream* ptr_log = nullptr, //!< output file stream for logging
                          const int rank = 0                //!< MPI rank
    );

    //! Read beta projectors from a pseudopotential file of UPF 2.0.1 format
    void read_beta_upf201(std::ifstream& ifs,               //!< input file stream from orbital file
                          std::ofstream* ptr_log = nullptr, //!< output file stream for logging
                          const int rank = 0                //!< MPI rank
    );

    //! extract the substring between a pair of quotation marks (for UPF v2.0.1)
    std::string trim201(std::string const& str);
};

#endif
