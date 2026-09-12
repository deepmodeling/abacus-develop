#pragma once

#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include <array>
#include <string>
#include <vector>

namespace LR_IO
{

class RI_kRlist
{
  public:
    K_Vectors* klist = nullptr;
    K_Vectors klist_coarse;
    std::array<int, 3> period;
    std::vector<std::array<int, 3>> Rlist;
    int nspin;
    RI_kRlist() = default;
    RI_kRlist(const UnitCell& ucell, K_Vectors* pkv, const int npspin_in,
              const std::string& in_dir, const std::string& out_dir, const int use_fine_kgrid);
    ~RI_kRlist() = default;
    void read_kpts_coarse(const std::string& file, const UnitCell& ucell,
                          K_Vectors* const klist, const std::string& out_dir);
    void read_kpts_fine(const std::string& file, const UnitCell& ucell,
                        K_Vectors* const klist, const bool is_weighted,
                        const std::string& out_dir);
};

} // namespace LR_IO
