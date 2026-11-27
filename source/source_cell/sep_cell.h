/**
 * @brief The Sep_Cell class is container for Sep potential.
 */

#ifndef SEP_CELL
#define SEP_CELL

#include "source_cell/sep.h"

#include <fstream>
#include <string>
#include <vector>

class Sep_Cell
{
  public:
    Sep_Cell() noexcept;
    ~Sep_Cell() noexcept;

    /**
     * @brief Sets the number of atom types and initializes internal vectors.
     *
     * @param ntype_in number of atom types
     */
    void init(const int ntype_in);

    /**
     * @brief Sets omega and tpiba2.
     *
     * @param omega_in unit cell volume
     * @param tpiba2_in tpiba squared
     */
    void set_omega(const double omega_in, const double tpiba2_in);

    /**
     * @brief Reads self potentials from STRU file and xx.sep files.
     *
     * @param ifpos input file stream
     * @param pp_dir pseudopotential directory
     * @param ofs_running output file stream for running log
     * @param ucell_atom_label atom labels from unit cell
     * @return true if successful, false otherwise
     */
    int read_sep_potentials(std::ifstream& ifpos,
                            const std::string& pp_dir,
                            std::ofstream& ofs_running,
                            std::vector<std::string>& ucell_atom_label);

#ifdef __MPI
    /**
     * @brief Broadcasts the Sep_Cell object to all processes.
     */
    void bcast_sep_cell();
#endif // __MPI

    /// @brief Getter methods
    const std::vector<SepPot>& get_seps() const
    {
        return seps;
    }
    int get_ntype() const
    {
        return ntype;
    }
    const std::vector<bool>& get_sep_enable() const
    {
        return sep_enable;
    }

    double get_omega() const
    {
        return omega;
    }

    double get_tpiba2() const
    {
        return tpiba2;
    }

  private:
    std::vector<SepPot> seps;     ///< Self potentials for each atom type
    int ntype;                    ///< number of atom types
    std::vector<bool> sep_enable; ///< Whether self potential is enabled for each atom type

    /// @brief unit cell data for VSep
    double omega;  ///< unit cell Volume
    double tpiba2; ///< tpiba ^ 2
};

#endif // SEP_CEll
