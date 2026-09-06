/**
 * @file klist_io.cpp
 * @brief this-free helpers extracted from K_Vectors (IBZ table formatting and
 *        line-mode k-point interpolation). Kept separate from klist.cpp so the
 *        logic is testable in isolation; klist.cpp only keeps thin wrappers.
 */
#include "klist_io.h"

#include "source_base/formatter.h"
#include "source_base/global_function.h"
#include "source_base/parallel_common.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/reciprocal_grid.h"
#include "source_cell/unitcell.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace KListIO
{

std::string ibz_kpt_table(const int nkstot,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                          const std::vector<int>& ibz_index,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz)
{
    std::stringstream ss;
    ss << " " << std::setw(40) << "nkstot"
       << " = " << nkstot << std::setw(66) << "ibzkpt" << std::endl;
    std::string table;
    table += "K-POINTS REDUCTION ACCORDING TO SYMMETRY\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s%12s%12s%12s\n",
                             "KPT",
                             "DIRECT_X",
                             "DIRECT_Y",
                             "DIRECT_Z",
                             "IBZ",
                             "DIRECT_X",
                             "DIRECT_Y",
                             "DIRECT_Z");
    for (int i = 0; i < nkstot; ++i)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8d%12.8f%12.8f%12.8f\n",
                                 i + 1,
                                 kvec_d[i].x,
                                 kvec_d[i].y,
                                 kvec_d[i].z,
                                 ibz_index[i] + 1,
                                 kvec_d_ibz[ibz_index[i]].x,
                                 kvec_d_ibz[ibz_index[i]].y,
                                 kvec_d_ibz[ibz_index[i]].z);
    }
    ss << table << std::endl;
    return ss.str();
}

std::string ibz_wk_table(const int nkstot_ibz,
                         const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                         const std::vector<double>& wk_ibz,
                         const std::vector<int>& ibz2bz)
{
    std::string table;
    table += "\n K-POINTS REDUCTION ACCORDING TO SYMMETRY\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s%8s\n", "IBZ", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT", "ibz2bz");
    for (int ik = 0; ik < nkstot_ibz; ik++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f%8d\n",
                                 ik + 1,
                                 kvec_d_ibz[ik].x,
                                 kvec_d_ibz[ik].y,
                                 kvec_d_ibz[ik].z,
                                 wk_ibz[ik],
                                 ibz2bz[ik]);
    }
    return table;
}

bool find_kpoints_header(std::ifstream& ifk)
{
    std::string word;
    while (ifk.good())
    {
        ifk >> word;
        // LiuXh add 20180416, fix bug in k-point file when the first line with comments
        ifk.ignore(150, '\n');
        if (word == "K_POINTS" || word == "KPOINTS" || word == "K")
        {
            return true;
        }
    }
    return false;
}

void read_kpt_list(std::ifstream& ifk,
                   const int nkstot,
                   std::vector<ModuleBase::Vector3<double>>& kvec,
                   std::vector<double>& wk)
{
    for (int i = 0; i < nkstot; i++)
    {
        ifk >> kvec[i].x >> kvec[i].y >> kvec[i].z;
        ModuleBase::GlobalFunc::READ_VALUE(ifk, wk[i]);
    }
}

LineK interp_line(std::ifstream& ifk, const int nks_special)
{
    // number of points to the next k points
    std::vector<int> nkl(nks_special, 0);

    // coordinates of special points.
    std::vector<ModuleBase::Vector3<double>> ks(nks_special);

    LineK out;
    std::vector<int> kpt_segids;
    int kpt_segid = 0;
    for (int iks = 0; iks < nks_special; iks++)
    {
        ifk >> ks[iks].x;
        ifk >> ks[iks].y;
        ifk >> ks[iks].z;
        ModuleBase::GlobalFunc::READ_VALUE(ifk, nkl[iks]);

        if (nkl[iks] <= 0)
        {
            ModuleBase::WARNING_QUIT("KListIO::interp_line",
                                     "Line-mode interpolation counts must be positive.");
        }
        out.nks_total += nkl[iks];
        /* ISSUE#3482: to distinguish different kline segments */
        if ((nkl[iks] == 1) && (iks != (nks_special - 1))) {
            kpt_segid++;
        }
        kpt_segids.push_back(kpt_segid);
    }
    if (nkl[nks_special - 1] != 1)
    {
        ModuleBase::WARNING_QUIT("KListIO::interp_line",
                                 "The final line-mode k-point must have an interpolation count of 1.");
    }

    out.kpts.resize(out.nks_total);
    out.segids.reserve(out.nks_total);

    int count = 0;
    for (int iks = 1; iks < nks_special; iks++)
    {
        double dxs = (ks[iks].x - ks[iks - 1].x) / nkl[iks - 1];
        double dys = (ks[iks].y - ks[iks - 1].y) / nkl[iks - 1];
        double dzs = (ks[iks].z - ks[iks - 1].z) / nkl[iks - 1];
        for (int is = 0; is < nkl[iks - 1]; is++)
        {
            out.kpts[count].x = ks[iks - 1].x + is * dxs;
            out.kpts[count].y = ks[iks - 1].y + is * dys;
            out.kpts[count].z = ks[iks - 1].z + is * dzs;
            out.segids.push_back(kpt_segids[iks - 1]); /* ISSUE#3482 */
            ++count;
        }
    }

    // deal with the last special k point.
    out.kpts[count].x = ks[nks_special - 1].x;
    out.kpts[count].y = ks[nks_special - 1].y;
    out.kpts[count].z = ks[nks_special - 1].z;
    out.segids.push_back(kpt_segids[nks_special - 1]); /* ISSUE#3482 */
    ++count;

    assert(count == out.nks_total);
    assert(out.segids.size() == static_cast<size_t>(out.nks_total)); /* ISSUE#3482 */
    return out;
}

void build_kstars(const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                  const std::vector<ModuleBase::Matrix3>& kgmatrix,
                  const int nrotkm,
                  const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                  const double epsilon,
                  const std::function<bool(double, double)>& equal,
                  std::vector<std::map<int, ModuleBase::Vector3<double>>>& kstars)
{
    const int nkstot = static_cast<int>(kvec_d.size());
    const int nkstot_ibz = static_cast<int>(kvec_d_ibz.size());
    kstars.resize(nkstot_ibz);

    ModuleBase::Vector3<double> kvec_rot;
    for (int i = 0; i < nkstot; ++i)
    {
        int exist_number = -1;
        int isym = 0;
        for (int j = 0; j < nrotkm; ++j)
        {
            kvec_rot = kvec_d[i] * kgmatrix[j];
            ModuleCell::restrict_kpt(kvec_rot, epsilon);
            for (int k = 0; k < nkstot_ibz; ++k)
            {
                if (equal(kvec_rot.x, kvec_d_ibz[k].x) && equal(kvec_rot.y, kvec_d_ibz[k].y)
                    && equal(kvec_rot.z, kvec_d_ibz[k].z))
                {
                    isym = j;
                    exist_number = k;
                    break;
                }
            }
            if (exist_number != -1)
            {
                break;
            }
        }
        kstars[exist_number].insert(std::make_pair(isym, kvec_d[i]));
    }
}

int append_time_reversal_ops(const ModuleSymmetry::Symmetry& symm,
                             std::vector<ModuleBase::Matrix3>& kgmatrix,
                             const int nrotkm)
{
    const ModuleBase::Matrix3 inv{-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0};

    bool include_inv = false;
    for (int i = 0; i < nrotkm; ++i)
    {
        if (kgmatrix[i] == inv)
        {
            include_inv = true;
        }
    }

    if (symm.magnetic_nspin4)
    {
        // (nspin=4, magnetic) Time reversal Theta reverses the magnetization,
        // so Theta alone is NOT a symmetry; only the antiunitary Theta*g
        // elements with g in the moment-reversing coset belong to the
        // Shubnikov group. The same index convention j + nrotk is decoded
        // in restore_dm. (nspin=2 is unaffected: there the antiunitary
        // operation is plain conjugation K, which leaves D_s(-k)=D_s*(k).)
        for (int j = 0; j < symm.nrotk_anti; ++j)
        {
            kgmatrix[j + symm.nrotk] = inv * symm.kgmatrix_anti[j];
        }
        return symm.nrotk + symm.nrotk_anti;
    }
    if (!include_inv)
    {
        for (int i = 0; i < symm.nrotk; ++i)
        {
            kgmatrix[i + symm.nrotk] = inv * symm.kgmatrix[i];
        }
        return 2 * symm.nrotk;
    }
    return nrotkm;
}

void pack_kpts(const std::vector<int>& isk,
               const std::vector<double>& wk,
               const std::vector<ModuleBase::Vector3<double>>& kvec_c,
               const std::vector<ModuleBase::Vector3<double>>& kvec_d,
               const std::vector<ModuleBase::Vector3<double>>& kvec_c_full,
               const int nkstot,
               std::vector<int>& isk_aux,
               std::vector<double>& wk_aux,
               std::vector<double>& kvec_c_aux,
               std::vector<double>& kvec_d_aux,
               std::vector<double>& kvec_c_full_aux)
{
    for (int ik = 0; ik < nkstot; ik++)
    {
        isk_aux[ik] = isk[ik];
        wk_aux[ik] = wk[ik];
        kvec_c_aux[3 * ik] = kvec_c[ik].x;
        kvec_c_aux[3 * ik + 1] = kvec_c[ik].y;
        kvec_c_aux[3 * ik + 2] = kvec_c[ik].z;
        kvec_d_aux[3 * ik] = kvec_d[ik].x;
        kvec_d_aux[3 * ik + 1] = kvec_d[ik].y;
        kvec_d_aux[3 * ik + 2] = kvec_d[ik].z;
        kvec_c_full_aux[3 * ik] = kvec_c_full[ik].x;
        kvec_c_full_aux[3 * ik + 1] = kvec_c_full[ik].y;
        kvec_c_full_aux[3 * ik + 2] = kvec_c_full[ik].z;
    }
}

void unpack_kpts(const std::vector<int>& isk_aux,
                 const std::vector<double>& wk_aux,
                 const std::vector<double>& kvec_c_aux,
                 const std::vector<double>& kvec_d_aux,
                 const std::vector<double>& kvec_c_full_aux,
                 const int nks,
                 const int startk,
                 std::vector<int>& isk,
                 std::vector<double>& wk,
                 std::vector<ModuleBase::Vector3<double>>& kvec_c,
                 std::vector<ModuleBase::Vector3<double>>& kvec_d,
                 std::vector<ModuleBase::Vector3<double>>& kvec_c_full)
{
    for (int i = 0; i < nks; i++)
    {
        // 3 is because each k point has three value:kx, ky, kz
        const int k_index = i + startk;
        kvec_c[i].x = kvec_c_aux[k_index * 3];
        kvec_c[i].y = kvec_c_aux[k_index * 3 + 1];
        kvec_c[i].z = kvec_c_aux[k_index * 3 + 2];
        kvec_d[i].x = kvec_d_aux[k_index * 3];
        kvec_d[i].y = kvec_d_aux[k_index * 3 + 1];
        kvec_d[i].z = kvec_d_aux[k_index * 3 + 2];
        kvec_c_full[i].x = kvec_c_full_aux[k_index * 3];
        kvec_c_full[i].y = kvec_c_full_aux[k_index * 3 + 1];
        kvec_c_full[i].z = kvec_c_full_aux[k_index * 3 + 2];
        wk[i] = wk_aux[k_index];
        isk[i] = isk_aux[k_index];
    }
}

void bcast_kstars(std::vector<std::map<int, ModuleBase::Vector3<double>>>& kstars,
                  const int nkstot,
                  const int my_rank)
{
    kstars.resize(nkstot);
    for (int ikibz = 0; ikibz < nkstot; ++ikibz)
    {
        int starsize = kstars[ikibz].size();
        Parallel_Common::bcast_int(starsize);
        auto ks = kstars[ikibz].begin();
        for (int ik = 0; ik < starsize; ++ik)
        {
            int isym = 0;
            ModuleBase::Vector3<double> ks_vec(0, 0, 0);
            if (my_rank == 0)
            {
                isym = ks->first;
                ks_vec = ks->second;
                ++ks;
            }
            Parallel_Common::bcast_int(isym);
            Parallel_Common::bcast_double(ks_vec.x);
            Parallel_Common::bcast_double(ks_vec.y);
            Parallel_Common::bcast_double(ks_vec.z);
            if (my_rank != 0)
            {
                kstars[ikibz].insert(std::make_pair(isym, ks_vec));
            }
        }
    }
}

void fill_full_kvec(const bool kc_done,
                    const bool kd_done,
                    const int nkstot_nospin,
                    const ModuleBase::Matrix3& reciprocal_vec,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                    std::vector<ModuleBase::Vector3<double>>& kvec_c_full)
{
    if (!kc_done && kd_done)
    {
        for (int ik = 0; ik < nkstot_nospin; ++ik)
        {
            kvec_c_full[ik] = kvec_d[ik] * reciprocal_vec;
        }
    }
    else if (kc_done && !kd_done)
    {
        for (int ik = 0; ik < nkstot_nospin; ++ik)
        {
            kvec_c_full[ik] = kvec_c[ik];
        }
    }
}

void build_ik2iktot(const int my_pool,
                    const std::vector<int>& startk_pool,
                    const int spin_mult,
                    const int nks,
                    const int nkstot,
                    std::vector<int>& ik2iktot)
{
    ik2iktot.resize(nks);
#ifdef __MPI
    if (spin_mult == 2)
    {
        for (int ik = 0; ik < nks / 2; ++ik)
        {
            ik2iktot[ik] = startk_pool[my_pool] + ik;
            ik2iktot[ik + nks / 2] = nkstot / 2 + startk_pool[my_pool] + ik;
        }
    }
    else
    {
        for (int ik = 0; ik < nks; ++ik)
        {
            ik2iktot[ik] = startk_pool[my_pool] + ik;
        }
    }
#else
    for (int ik = 0; ik < nks; ++ik)
    {
        ik2iktot[ik] = ik;
    }
#endif
}

void expand_spin_kpoints(const int spin_mult,
                         std::vector<ModuleBase::Vector3<double>>& kvec_c,
                         std::vector<ModuleBase::Vector3<double>>& kvec_d,
                         std::vector<double>& wk,
                         std::vector<int>& isk,
                         int& nks,
                         int& nkstot)
{
    //=========================================================================
    // on output: the number of points is doubled and xk and wk in the
    // first (nks/2) positions correspond to up spin
    // those in the second (nks/2) ones correspond to down spin
    // spin_mult can only be 1 or 2 here: K_Vectors::set() maps nspin=4
    // (non-collinear) to 1 before the k-list is built.
    //=========================================================================
    switch (spin_mult)
    {
    case 1:
        for (int ik = 0; ik < nks; ik++)
        {
            isk[ik] = 0;
        }
        break;

    case 2:
        for (int ik = 0; ik < nks; ik++)
        {
            kvec_c[ik + nks] = kvec_c[ik];
            kvec_d[ik + nks] = kvec_d[ik];
            wk[ik + nks] = wk[ik];
            isk[ik] = 0;
            isk[ik + nks] = 1;
        }

        nks *= 2;
        nkstot *= 2;
        break;
    }

    return;
}

void write_auto_kfile(const UnitCell& ucell,
                      const std::string& fn,
                      const bool gamma_only_local,
                      const double kspacing[3],
                      const std::string& kmesh_type,
                      const double koffset[3],
                      std::ofstream& ofs_warning)
{
    if (gamma_only_local)
    {
        ofs_warning << " Auto generating k-points file: " << fn << std::endl;
        std::ofstream ofs(fn.c_str());
        ofs << "K_POINTS" << std::endl;
        ofs << "0" << std::endl;
        ofs << "Gamma" << std::endl;
        ofs << "1 1 1 0 0 0" << std::endl;
        ofs.close();
    }
    else if (kspacing[0] > 0.0)
    {
        if (kspacing[1] <= 0 || kspacing[2] <= 0)
        {
            ModuleBase::WARNING_QUIT("K_Vectors", "kspacing should > 0");
        };
        // number of K points = max(1,int(|bi|/KSPACING+1))
        ModuleBase::Matrix3 btmp = ucell.G;
        double b1 = sqrt(btmp.e11 * btmp.e11 + btmp.e12 * btmp.e12 + btmp.e13 * btmp.e13);
        double b2 = sqrt(btmp.e21 * btmp.e21 + btmp.e22 * btmp.e22 + btmp.e23 * btmp.e23);
        double b3 = sqrt(btmp.e31 * btmp.e31 + btmp.e32 * btmp.e32 + btmp.e33 * btmp.e33);
        int nk1 = std::max(1, static_cast<int>(b1 * ModuleBase::TWO_PI / kspacing[0] / ucell.lat0 + 1));
        int nk2 = std::max(1, static_cast<int>(b2 * ModuleBase::TWO_PI / kspacing[1] / ucell.lat0 + 1));
        int nk3 = std::max(1, static_cast<int>(b3 * ModuleBase::TWO_PI / kspacing[2] / ucell.lat0 + 1));

        ofs_warning << " Generate k-points file according to KSPACING: " << fn << std::endl;
        std::ofstream ofs(fn.c_str());
        ofs << "K_POINTS" << std::endl;
        ofs << "0" << std::endl;
        if (kmesh_type == "mp")
        {
            ofs << "Monkhorst-Pack" << std::endl;
        }
        else
        {
            ofs << "Gamma" << std::endl;
        }
        ofs << nk1 << " " << nk2 << " " << nk3 << " " << koffset[0] << " " << koffset[1] << " "
            << koffset[2] << std::endl;
        ofs.close();
    }
}

} // namespace KListIO
