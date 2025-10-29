// The Sep_Cell class is container for Sep potential.

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

    // Sets the number of atom types and initializes internal vectors
    void init(const int ntype_in);

    void set_omega(const double omega_in, const double tpiba2_in);

    // Reads self potentials from STRU file and xx.sep files
    // Returns true if successful, false otherwise
    int read_sep_potentials(std::ifstream& ifpos,
                            const std::string& pp_dir,
                            std::ofstream& ofs_running,
                            std::vector<std::string>& ucell_atom_label);

#ifdef __MPI
    // Broadcasts the Sep_Cell object to all processes
    void bcast_sep_cell();
#endif // __MPI

    // Getter methods
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
    std::vector<SepPot> seps;     // Self potentials for each atom type
    int ntype;                    // number of atom types
    std::vector<bool> sep_enable; // Whether self potential is enabled for each atom type

    // unit cell data for VSep
    double omega;  // unit cell Volume
    double tpiba2; // tpiba ^ 2
};

// namespace GlobalC
// {
// extern Sep_Cell sep_cell;
// }

#endif // SEP_CEll
