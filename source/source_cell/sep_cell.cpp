#include "sep_cell.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_common.h"
#include "source_base/tool_title.h"

#include <algorithm>
#include <string>
#include <vector>

// namespace GlobalC
// {
// Sep_Cell sep_cell;
// }

Sep_Cell::Sep_Cell() noexcept : ntype(0), omega(0.0), tpiba2(0.0)
{
}

Sep_Cell::~Sep_Cell() noexcept = default;

void Sep_Cell::init(const int ntype_in)
{
    this->ntype = ntype_in;
    this->seps.resize(ntype);
    this->sep_enable.resize(ntype);
    std::fill(this->sep_enable.begin(), this->sep_enable.end(), false);
}

void Sep_Cell::set_omega(const double omega_in, const double tpiba2_in)
{
    this->omega = omega_in;
    this->tpiba2 = tpiba2_in;
}

/**
 * read sep potential files
 *
 * need to add following lines in STRU file, and order of elements must match ATOMIC_SPECIES.
 * SEP_FILES
 * symbol is_enable r_in r_out r_power enhence_a
 *
 * example
 * Li 0
 * F  1 F_pbe_50.sep 0.0 2.0 20.0 1.0
 */
int Sep_Cell::read_sep_potentials(std::ifstream& ifpos,
                                  const std::string& pp_dir,
                                  std::ofstream& ofs_running,
                                  std::vector<std::string>& ucell_atom_label)
{
    ModuleBase::TITLE("Sep_Cell", "read_sep_potentials");

    if (!ModuleBase::GlobalFunc::SCAN_BEGIN(ifpos, "SEP_FILES"))
    {
        GlobalV::ofs_running << "Cannot find SEP_FILES section in STRU" << std::endl;
        return false;
    }

    ifpos.ignore(300, '\n');

    for (int i = 0; i < this->ntype; ++i)
    {
        std::string one_line, atom_label;
        std::getline(ifpos, one_line);
        std::stringstream ss(one_line);

        // read the label of the atom
        bool enable_tmp;
        ss >> atom_label >> enable_tmp;

        // Validate atom label
        if (atom_label != ucell_atom_label[i])
        {
            GlobalV::ofs_running << "Sep potential and atom order do not match. "
                                 << "Expected: " << ucell_atom_label[i] << ", Got: " << atom_label << std::endl;
            return false;
        }
        this->sep_enable[i] = enable_tmp;
        if (this->sep_enable[i])
        {
            this->seps[i].is_enable = this->sep_enable[i];
            std::string sep_filename;
            ss >> sep_filename;
            ss >> this->seps[i].r_in >> this->seps[i].r_out >> this->seps[i].r_power >> this->seps[i].enhence_a;
            std::string sep_addr = pp_dir + sep_filename;
            std::ifstream sep_ifs(sep_addr.c_str(), std::ios::in);
            if (!sep_ifs)
            {
                GlobalV::ofs_running << "Cannot find sep potential file: " << sep_addr << std::endl;
                return false;
            }
            this->seps[i].read_sep(sep_ifs);
        }
    }

    return true;
}

#ifdef __MPI
void Sep_Cell::bcast_sep_cell()
{
    ModuleBase::TITLE("Sep_Cell", "bcast_sep_cell");
    Parallel_Common::bcast_int(this->ntype);

    if (GlobalV::MY_RANK != 0)
    {
        this->seps.resize(this->ntype);
        this->sep_enable.resize(this->ntype);
    }
    for (int i = 0; i < this->ntype; ++i)
    {
        bool tmp = false;
        if (GlobalV::MY_RANK == 0)
        {
            tmp = this->sep_enable[i];
        }
        Parallel_Common::bcast_bool(tmp);
        if (GlobalV::MY_RANK != 0)
        {
            this->sep_enable[i] = tmp;
        }
        this->seps[i].bcast_sep();
    }
}
#endif // __MPI
