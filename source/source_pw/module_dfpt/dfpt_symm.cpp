#include "dfpt_symm.h"

#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

std::vector<int> little_group_of_q(
    const ModuleBase::Matrix3 &latvec,
    const ModuleBase::Matrix3 *gmatrix,
    int nops,
    const ModuleBase::Vector3<double> &q_cart)
{
    ModuleBase::TITLE("DFPT_Symm", "little_group_of_q");
    ModuleBase::timer::start("DFPT_Symm", "little_group_of_q");
    std::vector<int> group;
    if (gmatrix == nullptr || nops <= 0)
    {
        ModuleBase::timer::end("DFPT_Symm", "little_group_of_q");
        return group;
    }
    // cartesian rotation images of the fractional operations, following
    // the symmetrize_mat3 construction (R = latvec^-1 G latvec)
    const ModuleBase::Matrix3 inv_lat = latvec.Inverse();
    // direct primitive vectors (rows of latvec) for the reciprocal test
    const ModuleBase::Vector3<double> a1(latvec.e11, latvec.e12, latvec.e13);
    const ModuleBase::Vector3<double> a2(latvec.e21, latvec.e22, latvec.e23);
    const ModuleBase::Vector3<double> a3(latvec.e31, latvec.e32, latvec.e33);
    const double tol = 1.0e-6; ///< empirical parameter: integrality roundoff of the double dot products
    for (int iop = 0; iop < nops; ++iop)
    {
        const ModuleBase::Matrix3 rot = inv_lat * gmatrix[iop] * latvec;
        const ModuleBase::Vector3<double> rq = rot * q_cart;
        for (int s = 1; s >= -1; s -= 2)
        {
            const ModuleBase::Vector3<double> diff = (s == 1) ? (rq - q_cart) : (rq + q_cart);
            const double d1 = diff * a1 / ModuleBase::TWO_PI;
            const double d2 = diff * a2 / ModuleBase::TWO_PI;
            const double d3 = diff * a3 / ModuleBase::TWO_PI;
            const bool recip = std::abs(d1 - std::round(d1)) < tol && std::abs(d2 - std::round(d2)) < tol &&
                               std::abs(d3 - std::round(d3)) < tol;
            if (recip)
            {
                group.push_back(iop);
                break; // either branch suffices; keep each operation once
            }
        }
    }
    ModuleBase::timer::end("DFPT_Symm", "little_group_of_q");
    return group;
}

} // namespace ModuleDFPT

namespace
{

/// wrap a fractional coordinate component to [0, 1) with an exact-1 snap
ModuleBase::Vector3<double> wrap01(const ModuleBase::Vector3<double> &v)
{
    const double c[3] = {v.x, v.y, v.z};
    double r[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
    {
        r[i] = c[i] - std::floor(c[i] + 1.0e-9);
        if (std::abs(r[i] - 1.0) < 1.0e-9)
        {
            r[i] = 0.0;
        }
    }
    return ModuleBase::Vector3<double>(r[0], r[1], r[2]);
}

/// lexicographic (x, y, z) ordering for representative selection
bool k_less(const ModuleBase::Vector3<double> &a, const ModuleBase::Vector3<double> &b)
{
    if (a.x != b.x)
    {
        return a.x < b.x;
    }
    if (a.y != b.y)
    {
        return a.y < b.y;
    }
    return a.z < b.z;
}

} // namespace

namespace ModuleDFPT
{

std::vector<KQPair> kq_pair_reps(
    const ModuleBase::Matrix3 &latvec,
    const ModuleBase::Matrix3 *gmatrix,
    int nops,
    const ModuleBase::Vector3<double> &q_frac,
    const std::vector<ModuleBase::Vector3<double>> &mesh_frac)
{
    ModuleBase::TITLE("DFPT_Symm", "kq_pair_reps");
    ModuleBase::timer::start("DFPT_Symm", "kq_pair_reps");
    std::vector<KQPair> pairs;
    const int nmesh = static_cast<int>(mesh_frac.size());
    if (gmatrix == nullptr || nops <= 0 || nmesh == 0)
    {
        ModuleBase::timer::end("DFPT_Symm", "kq_pair_reps");
        return pairs;
    }
    // reciprocal map k_cart = 2pi latvec^-1 k_frac (the columns of
    // 2pi latvec^-1 are the reciprocal primitive vectors b_i, satisfying
    // b_i . a_j = 2pi delta_ij with a_j the latvec rows); the transposed
    // variant 2pi (latvec^-1)^T differs for non-symmetric basis matrices
    // and was only accidentally right for the fcc tests
    const ModuleBase::Matrix3 inv_lat = latvec.Inverse();
    const ModuleBase::Matrix3 b_mat = inv_lat * ModuleBase::TWO_PI;
    const ModuleBase::Matrix3 lat = latvec;
    const ModuleBase::Vector3<double> q_cart = b_mat * q_frac;
    const ModuleBase::Vector3<double> a1(lat.e11, lat.e12, lat.e13);
    const ModuleBase::Vector3<double> a2(lat.e21, lat.e22, lat.e23);
    const ModuleBase::Vector3<double> a3(lat.e31, lat.e32, lat.e33);
    // unitary little co-group (s = +1 branch): R q = q + G
    std::vector<ModuleBase::Matrix3> unitary;
    const double tol = 1.0e-6;
    for (int iop = 0; iop < nops; ++iop)
    {
        const ModuleBase::Matrix3 rot = inv_lat * gmatrix[iop] * lat;
        const ModuleBase::Vector3<double> rq = rot * q_cart;
        const ModuleBase::Vector3<double> diff = rq - q_cart;
        const double d1 = diff * a1 / ModuleBase::TWO_PI;
        const double d2 = diff * a2 / ModuleBase::TWO_PI;
        const double d3 = diff * a3 / ModuleBase::TWO_PI;
        const bool recip = std::abs(d1 - std::round(d1)) < tol && std::abs(d2 - std::round(d2)) < tol &&
                           std::abs(d3 - std::round(d3)) < tol;
        if (recip)
        {
            unitary.push_back(rot);
        }
    }
    // cartesian images of the mesh for the orbit partition; membership is
    // compared as wrapped FRACTIONAL coordinates because a rotation image
    // may differ from a mesh point by a reciprocal-lattice vector (the
    // wrapping orbits): a raw cartesian comparison silently splits those
    // orbits, which the fcc 2x2x2 tests missed only because their star
    // images land exactly on mesh points
    std::vector<ModuleBase::Vector3<double>> mesh_cart(nmesh);
    std::vector<ModuleBase::Vector3<double>> mesh_wrap(nmesh);
    for (int ik = 0; ik < nmesh; ++ik)
    {
        mesh_cart[ik] = b_mat * mesh_frac[ik];
        mesh_wrap[ik] = wrap01(mesh_frac[ik]);
    }
    const ModuleBase::Matrix3 inv_b = b_mat.Inverse();
    std::vector<bool> seen(nmesh, false);
    const double match_tol = 1.0e-6;
    for (int ik = 0; ik < nmesh; ++ik)
    {
        if (seen[ik])
        {
            continue;
        }
        // orbit of ik; membership against the mesh modulo the reciprocal
        // lattice through the wrapped fractional images
        std::vector<int> members;
        for (int idx = 0; idx < nmesh; ++idx)
        {
            if (seen[idx] && idx != ik)
            {
                continue;
            }
            for (const ModuleBase::Matrix3 &rot : unitary)
            {
                const ModuleBase::Vector3<double> image = inv_b * (rot * mesh_cart[idx]);
                const ModuleBase::Vector3<double> image_wrap = wrap01(image);
                if (std::abs(image_wrap.x - mesh_wrap[ik].x) < match_tol &&
                    std::abs(image_wrap.y - mesh_wrap[ik].y) < match_tol &&
                    std::abs(image_wrap.z - mesh_wrap[ik].z) < match_tol)
                {
                    members.push_back(idx);
                    break;
                }
            }
        }
        // canonical representative: lexicographically smallest k_frac
        int rep = members.front();
        for (const int idx : members)
        {
            if (k_less(mesh_frac[idx], mesh_frac[rep]))
            {
                rep = idx;
            }
        }
        for (const int idx : members)
        {
            seen[idx] = true;
        }
        KQPair pair;
        pair.k_frac = wrap01(mesh_frac[rep]);
        pair.kq_frac = wrap01(mesh_frac[rep] + q_frac);
        pair.orbit_size = static_cast<int>(members.size());
        pairs.push_back(pair);
    }
    ModuleBase::timer::end("DFPT_Symm", "kq_pair_reps");
    return pairs;
}

} // namespace ModuleDFPT
