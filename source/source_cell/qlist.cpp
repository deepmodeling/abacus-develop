// ============================================================
// QList: q-point mesh generation and star reduction.
// ============================================================

#include "qlist.h"

#include "module_symmetry/symmetry.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/formatter.h"
#include "source_base/tool_quit.h"
#include "unitcell.h"
#include <algorithm>
#include <cassert>
#include <iomanip>

namespace ModuleCell {

QList::QList() {}

QList::~QList() {}

void QList::generate_mesh(UnitCell& ucell, ModuleSymmetry::Symmetry& symm,
                          const std::vector<int>& mp_grid, bool use_irreps) {
    if (mp_grid.size() != 3)
    {
        ModuleBase::WARNING_QUIT("QList::generate_mesh", "mp_grid must have three components.");
    }

    this->is_mp = true;
    this->nmp[0] = mp_grid[0];
    this->nmp[1] = mp_grid[1];
    this->nmp[2] = mp_grid[2];

    // Gamma-centered Monkhorst-Pack q mesh (k_type = 0), zero offset.
    const double offset[3] = {0.0, 0.0, 0.0};
    this->Monkhorst_Pack(this->nmp, offset, 0);

    this->nkstot_nospin = this->nkstot;
    this->nks = this->nkstot;

    // Star reduction: always use symmetry, always include the -q partner.
    bool match = true;
    std::string skpt;
    this->reduce_by_symmetry(ucell, symm, true, skpt, match, GlobalV::MY_RANK, GlobalV::ofs_running);
    if (!match)
    {
        ModuleBase::WARNING("QList::generate_mesh",
                            "Reciprocal lattice is incompatible with the real-space lattice. "
                            "Falling back to the unreduced q-point mesh.");
        this->nkstot = this->nks = this->nkstot_nospin;
    }

    // weights sum to 1 (average over the full Brillouin zone)
    this->normalize_wk(1);

    // Cartesian coordinates of the reduced q-point list (from the direct ones).
    // The reciprocal lattice is stored in ucell.G (columns are the reciprocal
    // primitive vectors), so q_cart = q_direct * G.
    this->kvec_d2c(ucell.G);
    this->kc_done = true;

    // little-group irreducible-representation data (opt-in)
    if (use_irreps)
    {
        this->get_irreps(ucell, symm);
    }
    else
    {
        this->nirr_.clear();
        this->irrep_modes_.clear();
    }
}

void QList::read_from_file(const std::string& filename, UnitCell& ucell) {
    std::ifstream ifq(filename.c_str());
    if (!ifq)
    {
        ModuleBase::WARNING("QList::read_from_file", "Can not find the q-points file.");
        this->nkstot = this->nks = 0;
        return;
    }

    ifq >> std::setiosflags(std::ios::uppercase);
    ifq.clear();
    ifq.seekg(0);

    // find the "Q_POINTS" (or "QPOINTS" / "Q") header, mirroring read_kpoints
    std::string word;
    std::string qword;
    int ierr = 0;
    while (ifq.good())
    {
        ifq >> word;
        ifq.ignore(150, '\n');
        if (word == "Q_POINTS" || word == "QPOINTS" || word == "Q")
        {
            ierr = 1;
            break;
        }
        ifq.rdstate();
    }
    if (ierr == 0)
    {
        ModuleBase::WARNING("QList::read_from_file", "symbol Q_POINTS not found.");
        this->nkstot = this->nks = 0;
        return;
    }

    ModuleBase::GlobalFunc::READ_VALUE(ifq, this->nkstot);
    this->k_nkstot = this->nkstot;
    ModuleBase::GlobalFunc::READ_VALUE(ifq, qword);
    this->k_kword = qword;

    const int max_qpoints = 100000;
    if (this->nkstot < 0 || this->nkstot > max_qpoints)
    {
        ModuleBase::WARNING("QList::read_from_file", "nkstot is negative or greater than MAX_QPOINTS.");
        this->nkstot = this->nks = 0;
        return;
    }

    int q_type = 0;
    if (this->nkstot == 0) // Monkhorst-Pack mesh
    {
        this->is_mp = true;
        if (qword == "Gamma")
        {
            q_type = 0;
        }
        else if (qword == "Monkhorst-Pack" || qword == "MP" || qword == "mp")
        {
            q_type = 1;
        }
        else
        {
            ModuleBase::WARNING("QList::read_from_file", "neither Gamma nor Monkhorst-Pack.");
            this->nkstot = this->nks = 0;
            return;
        }

        ifq >> this->nmp[0] >> this->nmp[1] >> this->nmp[2];
        double koffset[3] = {0.0, 0.0, 0.0};
        if (!(ifq >> koffset[0] >> koffset[1] >> koffset[2]))
        {
            ModuleBase::WARNING("QList::read_from_file", "Missing q-point offsets in the q-points file.");
        }
        this->Monkhorst_Pack(this->nmp, koffset, q_type);
    }
    else // explicit list or line path
    {
        if (qword == "Cartesian" || qword == "C")
        {
            this->renew(this->nkstot);
            for (int i = 0; i < this->nkstot; ++i)
            {
                ifq >> kvec_c[i].x >> kvec_c[i].y >> kvec_c[i].z;
                ModuleBase::GlobalFunc::READ_VALUE(ifq, wk[i]);
            }
            this->kc_done = true;
        }
        else if (qword == "Direct" || qword == "D")
        {
            this->renew(this->nkstot);
            for (int i = 0; i < this->nkstot; ++i)
            {
                ifq >> kvec_d[i].x >> kvec_d[i].y >> kvec_d[i].z;
                ModuleBase::GlobalFunc::READ_VALUE(ifq, wk[i]);
            }
            this->kd_done = true;
        }
        else if (qword == "Line_Direct" || qword == "L" || qword == "Line")
        {
            interpolate_q_between(ifq, kvec_d);
            std::for_each(wk.begin(), wk.end(), [](double& d) { d = 1.0; });
            this->kd_done = true;
        }
        else if (qword == "Line_Cartesian")
        {
            interpolate_q_between(ifq, kvec_c);
            std::for_each(wk.begin(), wk.end(), [](double& d) { d = 1.0; });
            this->kc_done = true;
        }
        else
        {
            ModuleBase::WARNING("QList::read_from_file", "neither Cartesian nor Direct qpoint.");
            this->nkstot = this->nks = 0;
            return;
        }
    }

    this->nkstot_nospin = this->nks = this->nkstot;

    // complement the coordinates: fill the missing representation
    if (!this->kc_done && this->kd_done)
    {
        this->kvec_d2c(ucell.G);
        this->kc_done = true;
    }
    else if (this->kc_done && !this->kd_done)
    {
        this->kvec_c2d(ucell.latvec);
        this->kd_done = true;
    }

    // weights: a mesh or explicit list is normalized to sum 1; a line path
    // keeps its unnormalized weights (each point weight 1)
    if (this->k_kword != "Line_Direct" && this->k_kword != "L" && this->k_kword != "Line"
        && this->k_kword != "Line_Cartesian")
    {
        this->normalize_wk(1);
    }

    // no symmetry reduction (no symmetry object in this interface), so no
    // little-group irrep decomposition either; the DFPT driver requires at
    // least the fallback fully-symmetric placeholder (nirr = 1, empty mode
    // basis -> solve the full 3N displacement basis)
    this->nirr_.assign(this->nkstot, 1);
    this->irrep_modes_.assign(this->nkstot, std::vector<std::vector<int>>(1));
}

void QList::interpolate_q_between(std::ifstream& ifq, std::vector<ModuleBase::Vector3<double>>& qvec) {
    const int nqs_special = this->nkstot;
    std::vector<int> nql(nqs_special, 0);
    std::vector<ModuleBase::Vector3<double>> qs(nqs_special);

    // recalculate nkstot
    this->nkstot = 0;
    this->kl_segids.clear();
    this->kl_segids.shrink_to_fit();
    int qpt_segid = 0;
    for (int iqs = 0; iqs < nqs_special; ++iqs)
    {
        ifq >> qs[iqs].x;
        ifq >> qs[iqs].y;
        ifq >> qs[iqs].z;
        ModuleBase::GlobalFunc::READ_VALUE(ifq, nql[iqs]);
        if (nql[iqs] <= 0)
        {
            ModuleBase::WARNING_QUIT("QList::interpolate_q_between",
                                     "Line-mode interpolation counts must be positive.");
        }
        this->nkstot += nql[iqs];
        if ((nql[iqs] == 1) && (iqs != (nqs_special - 1)))
        {
            ++qpt_segid;
        }
        this->kl_segids.push_back(qpt_segid);
    }
    if (nql[nqs_special - 1] != 1)
    {
        ModuleBase::WARNING_QUIT("QList::interpolate_q_between",
                                 "The final line-mode q-point must have an interpolation count of 1.");
    }

    this->renew(this->nkstot);

    int count = 0;
    for (int iqs = 1; iqs < nqs_special; ++iqs)
    {
        double dxs = (qs[iqs].x - qs[iqs - 1].x) / nql[iqs - 1];
        double dys = (qs[iqs].y - qs[iqs - 1].y) / nql[iqs - 1];
        double dzs = (qs[iqs].z - qs[iqs - 1].z) / nql[iqs - 1];
        for (int is = 0; is < nql[iqs - 1]; ++is)
        {
            qvec[count].x = qs[iqs - 1].x + is * dxs;
            qvec[count].y = qs[iqs - 1].y + is * dys;
            qvec[count].z = qs[iqs - 1].z + is * dzs;
            ++count;
        }
    }
    qvec[count].x = qs[nqs_special - 1].x;
    qvec[count].y = qs[nqs_special - 1].y;
    qvec[count].z = qs[nqs_special - 1].z;
    ++count;
    assert(count == this->nkstot);
}

void QList::print_qlists(std::ofstream& ofs) const {
    ModuleBase::TITLE("QList", "print_qlists");
    const int nq = this->nks;
    if (this->nkstot < nq)
    {
        ModuleBase::WARNING_QUIT("QList::print_qlists", "nkstot < nks");
    }
    std::string table;
    table += " Q-POINTS CARTESIAN COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "QPOINTS", "CARTESIAN_X", "CARTESIAN_Y", "CARTESIAN_Z", "WEIGHT");
    for (int i = 0; i < nq; i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 this->kvec_c[i].x,
                                 this->kvec_c[i].y,
                                 this->kvec_c[i].z,
                                 this->wk[i]);
    }
    ofs << "\n" << table << std::endl;

    table.clear();
    table += " Q-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "QPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < nq; i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 this->kvec_d[i].x,
                                 this->kvec_d[i].y,
                                 this->kvec_d[i].z,
                                 this->wk[i]);
    }
    ofs << "\n" << table << std::endl;
    return;
}

std::vector<int> QList::get_irrep_modes(int q_idx, int irrep_idx) const {
    if (q_idx < 0 || q_idx >= static_cast<int>(this->nirr_.size()))
    {
        return std::vector<int>();
    }
    if (irrep_idx < 0 || irrep_idx >= this->nirr_[q_idx])
    {
        return std::vector<int>();
    }
    return this->irrep_modes_[q_idx][irrep_idx];
}

void QList::reduce_by_symmetry(const UnitCell& ucell,
                               const ModuleSymmetry::Symmetry& symm,
                               bool use_symm,
                               std::string& skpt,
                               bool& match,
                               const int my_rank,
                               std::ofstream& ofs_running) {
    (void)skpt;
    (void)my_rank;
    (void)ofs_running;
    // q-points are spin-free: build the point-group operations and always
    // double them by the time-reversal operation -q (no magnetic group).
    std::vector<ModuleBase::Matrix3> kgmatrix(48 * 2);
    ModuleBase::Matrix3 inv(-1, 0, 0, 0, -1, 0, 0, 0, -1);

    ModuleBase::Matrix3 q_vec; // k-lattice basis of the q mesh
    int nrotkm = 0;
    if (!this->build_star_ops(ucell, symm, use_symm, q_vec, kgmatrix, nrotkm))
    {
        match = false;
        return;
    }
    if (nrotkm == 0)
    {
        // no operations to apply: the mesh stays unreduced
        match = true;
        return;
    }

    bool include_inv = false;
    for (int i = 0; i < nrotkm; ++i)
    {
        if (kgmatrix[i] == inv)
        {
            include_inv = true;
            break;
        }
    }
    if (!include_inv)
    {
        for (int i = 0; i < nrotkm; ++i)
        {
            kgmatrix[i + nrotkm] = inv * kgmatrix[i];
        }
        nrotkm *= 2;
    }

    std::vector<ModuleBase::Matrix3> kkmatrix(nrotkm);
    symm.gmatrix_convert(kgmatrix.data(), kkmatrix.data(), nrotkm, ucell.G, q_vec);

    std::vector<ModuleBase::Vector3<double>> qvec_ibz;
    std::vector<double> wk_ibz;
    std::vector<int> ibz_index;
    std::vector<int> ibz2bz;
    this->reduce_ibz(kgmatrix.data(), nrotkm, ucell.G, q_vec, kkmatrix.data(), symm.epsilon, qvec_ibz, wk_ibz, ibz_index, ibz2bz);

    // update the reduced q-point list (no spin expansion)
    const int nq_ibz = qvec_ibz.size();
    this->nkstot = this->nks = nq_ibz;
    this->kvec_d.resize(this->nkstot);
    this->wk.resize(this->nkstot);
    for (int i = 0; i < this->nkstot; ++i)
    {
        this->kvec_d[i] = qvec_ibz[i];
        this->wk[i] = wk_ibz[i];
    }
    this->kd_done = true;
    this->kc_done = false;

    match = true;
    return;
}

void QList::get_irreps(const UnitCell& ucell, const ModuleSymmetry::Symmetry& symm) {
    (void)ucell;

    // Decompose each q-point via its little group (placeholder: one
    // fully-symmetric A1 irrep per q-point; the LittleGroup projection-operator
    // basis is filled in a later iteration).
    this->nirr_.assign(this->nkstot, 0);
    this->irrep_modes_.assign(this->nkstot, std::vector<std::vector<int>>());
    for (int iq = 0; iq < this->nkstot; ++iq)
    {
        this->little_group_.set_q(this->kvec_d[iq], symm);
        this->nirr_[iq] = this->little_group_.get_nirr();
        this->irrep_modes_[iq].resize(this->nirr_[iq]);
        for (int iirr = 0; iirr < this->nirr_[iq]; ++iirr)
        {
            this->irrep_modes_[iq][iirr] = this->little_group_.get_mode_basis(iirr);
        }
    }
}

} // namespace ModuleCell