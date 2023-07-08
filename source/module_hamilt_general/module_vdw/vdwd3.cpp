//==========================================================
// AUTHOR : Yuyang Ji
// DATE : 2019-04-22
// UPDATE : 2021-4-19
//==========================================================

#include "vdwd3.h"

#include "module_base/constants.h"
#include "module_base/element_name.h"
#include "module_base/global_function.h"
#include "module_base/timer.h"

namespace vdw
{

void Vdwd3::init()
{
    lat_.resize(3);
    lat_[0] = ucell_.a1 * ucell_.lat0;
    lat_[1] = ucell_.a2 * ucell_.lat0;
    lat_[2] = ucell_.a3 * ucell_.lat0;

    std::vector<double> at_kind = atom_kind();
    iz_.reserve(ucell_.nat);
    xyz_.reserve(ucell_.nat);
    for (size_t it = 0; it != ucell_.ntype; it++)
        for (size_t ia = 0; ia != ucell_.atoms[it].na; ia++)
        {
            iz_.emplace_back(at_kind[it]);
            xyz_.emplace_back(ucell_.atoms[it].tau[ia] * ucell_.lat0);
        }

    std::vector<double> tau_max(3);
    if (para_.model() == "radius")
    {
        rep_vdw_.resize(3);
        set_criteria(para_.rthr2(), lat_, tau_max);
        for (size_t i = 0; i < 3; i++)
            rep_vdw_[i] = std::ceil(tau_max[i]);
    }
    else if (para_.model() == "period")
        rep_vdw_ = {para_.period().x, para_.period().y, para_.period().z};

    rep_cn_.resize(3);
    set_criteria(para_.cn_thr2(), lat_, tau_max);
    for (size_t i = 0; i < 3; i++)
        rep_cn_[i] = ceil(tau_max[i]);
}

void Vdwd3::set_criteria(double rthr, const std::vector<ModuleBase::Vector3<double>> &lat, std::vector<double> &tau_max)
{
    tau_max.resize(3);
    double r_cutoff = std::sqrt(rthr);
    ModuleBase::Vector3<double> norm1 = (lat_[1] ^ lat_[2]).normalize();
    ModuleBase::Vector3<double> norm2 = (lat_[2] ^ lat_[0]).normalize();
    ModuleBase::Vector3<double> norm3 = (lat_[0] ^ lat_[1]).normalize();
    double cos10 = norm1 * lat_[0];
    double cos21 = norm2 * lat_[1];
    double cos32 = norm3 * lat_[2];
    tau_max[0] = std::abs(r_cutoff / cos10);
    tau_max[1] = std::abs(r_cutoff / cos21);
    tau_max[2] = std::abs(r_cutoff / cos32);
}

std::vector<double> Vdwd3::atom_kind()
{
    std::vector<double> atom_kind(ucell_.ntype);
    for (size_t i = 0; i != ucell_.ntype; i++)
        for (int j = 0; j != ModuleBase::element_name.size(); j++)
            if (ucell_.atoms[i].ncpp.psd == ModuleBase::element_name[j])
            {
                atom_kind[i] = j;
                break;
            }
    return atom_kind;
}

void Vdwd3::cal_energy()
{
    ModuleBase::TITLE("Vdwd3", "cal_energy");
    ModuleBase::timer::tick("Vdwd3", "cal_energy");
    init();

    int ij;
    double c6 = 0.0, c8 = 0.0, r2 = 0.0, r6 = 0.0, r8 = 0.0, rr = 0.0, damp6 = 0.0, damp8 = 0.0;
    double e6 = 0.0, e8 = 0.0, eabc = 0.0;
    std::vector<double> cc6ab(ucell_.nat * ucell_.nat), cn(ucell_.nat);
    pbc_ncoord(cn);
    ModuleBase::Vector3<double> tau;
    if (para_.version() == "d3_0") // DFT-D3(zero-damping)
    {
        double tmp;
        for (int iat = 0; iat != ucell_.nat - 1; iat++)
            for (int jat = iat + 1; jat != ucell_.nat; jat++)
            {
                get_c6(iz_[iat], iz_[jat], cn[iat], cn[jat], c6);
                if (para_.abc())
                {
                    ij = lin(iat, jat);
                    cc6ab[ij] = std::sqrt(c6);
                }
                for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                    for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                        for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                        {
                            tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                                  + static_cast<double>(tauz) * lat_[2];
                            r2 = (xyz_[iat] - xyz_[jat] + tau).norm2();
                            if (r2 > para_.rthr2())
                                continue;
                            rr = para_.r0ab()[iz_[iat]][iz_[jat]] / std::sqrt(r2);
                            // zero-damping function
                            tmp = para_.rs6() * rr;
                            damp6 = 1.0 / (1.0 + 6.0 * std::pow(tmp, para_.alp6()));
                            tmp = para_.rs18() * rr;
                            damp8 = 1.0 / (1.0 + 6.0 * std::pow(tmp, para_.alp8()));

                            r6 = std::pow(r2, 3);
                            e6 += damp6 / r6 * c6;

                            c8 = 3.0 * para_.r2r4()[iz_[iat]] * para_.r2r4()[iz_[jat]] * c6;
                            r8 = r6 * r2;
                            e8 += c8 * damp8 / r8;
                        } // end tau
            } // end jat

        for (int iat = 0; iat != ucell_.nat; iat++)
        {
            int jat = iat;
            get_c6(iz_[iat], iz_[jat], cn[iat], cn[jat], c6);
            if (para_.abc())
            {
                ij = lin(iat, jat);
                cc6ab[ij] = std::sqrt(c6);
            }
            for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                    for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                    {
                        if (taux == 0 && tauy == 0 && tauz == 0)
                            continue;
                        tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                              + static_cast<double>(tauz) * lat_[2];
                        r2 = tau.norm2();
                        if (r2 > para_.rthr2())
                            continue;
                        rr = para_.r0ab()[iz_[iat]][iz_[jat]] / std::sqrt(r2);

                        // zero-damping function
                        tmp = para_.rs6() * rr;
                        damp6 = 1.0 / (1.0 + 6.0 * std::pow(tmp, para_.alp6()));
                        tmp = para_.rs18() * rr;
                        damp8 = 1.0 / (1.0 + 6.0 * std::pow(tmp, para_.alp8()));

                        r6 = std::pow(r2, 3);
                        e6 += damp6 / r6 * c6 * 0.5;

                        c8 = 3.0 * para_.r2r4()[iz_[iat]] * para_.r2r4()[iz_[jat]] * c6;
                        r8 = r6 * r2;
                        e8 += c8 * damp8 / r8 * 0.5;
                    } // end tau
        } // end iat
    } // end d3_0
    else if (para_.version() == "d3_bj") // DFT-D3(BJ-damping)
    {
        double r42;
        for (int iat = 0; iat != ucell_.nat; iat++)
        {
            for (int jat = iat + 1; jat != ucell_.nat; jat++)
            {
                get_c6(iz_[iat], iz_[jat], cn[iat], cn[jat], c6);
                if (para_.abc())
                {
                    ij = lin(iat, jat);
                    cc6ab[ij] = std::sqrt(c6);
                }
                for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                    for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                        for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                        {
                            tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                                  + static_cast<double>(tauz) * lat_[2];
                            r2 = (xyz_[iat] - xyz_[jat] + tau).norm2();
                            if (r2 > para_.rthr2())
                                continue;
                            rr = para_.r0ab()[iz_[iat]][iz_[jat]] / std::sqrt(r2);

                            // BJ-damping function
                            r42 = para_.r2r4()[iz_[iat]] * para_.r2r4()[iz_[jat]];
                            damp6 = std::pow((para_.rs6() * std::sqrt(3.0 * r42) + para_.rs18()), 6);
                            damp8 = std::pow((para_.rs6() * std::sqrt(3.0 * r42) + para_.rs18()), 8);

                            r6 = std::pow(r2, 3);
                            e6 += c6 / (r6 + damp6);

                            c8 = 3.0 * c6 * r42;
                            r8 = r6 * r2;
                            e8 += c8 / (r8 + damp8);
                        } // end tau
            } // end jat
            int jat = iat;
            get_c6(iz_[iat], iz_[jat], cn[iat], cn[jat], c6);
            r42 = para_.r2r4()[iz_[iat]] * para_.r2r4()[iz_[jat]];
            damp6 = std::pow((para_.rs6() * std::sqrt(3.0 * r42) + para_.rs18()), 6);
            damp8 = std::pow((para_.rs6() * std::sqrt(3.0 * r42) + para_.rs18()), 8);
            if (para_.abc())
            {
                ij = lin(iat, jat);
                cc6ab[ij] = std::sqrt(c6);
            }
            for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                    for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                    {
                        if (taux == 0 && tauy == 0 && tauz == 0)
                            continue;
                        tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                              + static_cast<double>(tauz) * lat_[2];
                        r2 = tau.norm2();
                        if (r2 > para_.rthr2())
                            continue;
                        rr = para_.r0ab()[iz_[iat]][iz_[jat]] / std::sqrt(r2);

                        r6 = std::pow(r2, 3);
                        e6 += c6 / (r6 + damp6) * 0.5;

                        c8 = 3.0 * c6 * r42;
                        r8 = r6 * r2;
                        e8 += c8 / (r8 + damp8) * 0.5;
                    } // end tau
        } // end iat
    } // end d3_bj

    if (para_.abc())
    {
        pbc_three_body(iz_, lat_, xyz_, rep_cn_, cc6ab, eabc);
    }
    energy_ = (-para_.s6() * e6 - para_.s18() * e8 - eabc) * 2;
    ModuleBase::timer::tick("Vdwd3", "cal_energy");
}

void Vdwd3::cal_force()
{
    ModuleBase::TITLE("Vdwd3", "cal_force");
    ModuleBase::timer::tick("Vdwd3", "cal_force");
    init();

    force_.clear();
    force_.resize(ucell_.nat);

    std::vector<ModuleBase::Vector3<double>> g;
    g.clear();
    g.resize(ucell_.nat);
    ModuleBase::matrix smearing_sigma(3, 3);

    pbc_gdisp(g, smearing_sigma);

    for (size_t iat = 0; iat != ucell_.nat; iat++)
        force_[iat] = -2.0 * g[iat];

    ModuleBase::timer::tick("Vdwd3", "cal_force");
}

void Vdwd3::cal_stress()
{
    ModuleBase::TITLE("Vdwd3", "cal_stress");
    ModuleBase::timer::tick("Vdwd3", "cal_stress");
    init();

    std::vector<ModuleBase::Vector3<double>> g;
    g.clear();
    g.resize(ucell_.nat);
    ModuleBase::matrix smearing_sigma(3, 3);

    pbc_gdisp(g, smearing_sigma);

    stress_ = ModuleBase::Matrix3(2.0 * smearing_sigma(0, 0), 2.0 * smearing_sigma(0, 1), 2.0 * smearing_sigma(0, 2),
                                  2.0 * smearing_sigma(1, 0), 2.0 * smearing_sigma(1, 1), 2.0 * smearing_sigma(1, 2),
                                  2.0 * smearing_sigma(2, 0), 2.0 * smearing_sigma(2, 1), 2.0 * smearing_sigma(2, 2))
              / ucell_.omega;
    ModuleBase::timer::tick("Vdwd3", "cal_stress");
}

void Vdwd3::get_c6(int iat, int jat, double nci, double ncj, double &c6)
{
    double c6mem = -1e99, rsum = 0.0, csum = 0.0, r_save = 1e99;
    double cn1, cn2, r;
    for (size_t i = 0; i != para_.mxc()[iat]; i++)
        for (size_t j = 0; j != para_.mxc()[jat]; j++)
        {
            c6 = para_.c6ab()[0][j][i][jat][iat];
            if (c6 > 0)
            {
                cn1 = para_.c6ab()[1][j][i][jat][iat];
                cn2 = para_.c6ab()[2][j][i][jat][iat];
                r = std::pow((cn1 - nci), 2) + std::pow((cn2 - ncj), 2);
                if (r < r_save)
                {
                    r_save = r;
                    c6mem = c6;
                }
                double tmp1 = exp(para_.k3() * r);
                rsum += tmp1;
                csum += tmp1 * c6;
            }
        }
    c6 = (rsum > 1e-99) ? csum / rsum : c6mem;
}

void Vdwd3::pbc_ncoord(std::vector<double> &cn)
{
    for (size_t i = 0; i != ucell_.nat; i++)
    {
        double xn = 0.0;
        ModuleBase::Vector3<double> tau;
        double r2, rr;
        for (size_t iat = 0; iat != ucell_.nat; iat++)
            for (int taux = -rep_cn_[0]; taux <= rep_cn_[0]; taux++)
                for (int tauy = -rep_cn_[1]; tauy <= rep_cn_[1]; tauy++)
                    for (int tauz = -rep_cn_[2]; tauz <= rep_cn_[2]; tauz++)
                    {
                        if (iat == i && taux == 0 && tauy == 0 && tauz == 0)
                            continue;
                        tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                              + static_cast<double>(tauz) * lat_[2];
                        r2 = (xyz_[iat] - xyz_[i] + tau).norm2();
                        if (r2 > para_.cn_thr2())
                            continue;
                        rr = (para_.rcov()[iz_[i]] + para_.rcov()[iz_[iat]]) / std::sqrt(r2);
                        xn += 1.0 / (1.0 + exp(-para_.k1() * (rr - 1.0)));
                    }
        cn[i] = xn;
    }
}

void Vdwd3::pbc_three_body(const std::vector<int> &iz,
                           const std::vector<ModuleBase::Vector3<double>> &lat,
                           const std::vector<ModuleBase::Vector3<double>> &xyz,
                           const std::vector<int> &rep_cn,
                           const std::vector<double> &cc6ab,
                           double &eabc)
{
    double sr9 = 0.75, alp9 = -16.0;
    int ij, ik, jk;
    double r0ij, r0ik, r0jk, c9, rij2, rik2, rjk2, rr0ij, rr0ik, rr0jk, geomean, fdamp, tmp1, tmp2, tmp3, tmp4, ang;
    ModuleBase::Vector3<double> ijvec, ikvec, jkvec, jtau, ktau;
    std::vector<double> repmin(3), repmax(3);
    for (int iat = 2; iat != ucell_.nat; iat++)
        for (int jat = 1; jat != iat; jat++)
        {
            ijvec = xyz_[jat] - xyz_[iat];
            ij = lin(iat, jat);
            r0ij = para_.r0ab()[iz_[jat]][iz_[iat]];
            for (int kat = 0; kat != jat; kat++)
            {
                ik = lin(iat, kat);
                jk = lin(jat, kat);
                ikvec = xyz_[kat] - xyz_[iat];
                jkvec = xyz_[kat] - xyz_[jat];
                c9 = -cc6ab[ij] * cc6ab[ik] * cc6ab[jk];

                r0ik = para_.r0ab()[iz_[kat]][iz_[iat]];
                r0jk = para_.r0ab()[iz_[kat]][iz_[jat]];

                for (int jtaux = -rep_cn_[0]; jtaux <= rep_cn_[0]; jtaux++)
                {
                    repmin[0] = std::max(-rep_cn_[0], jtaux - rep_cn_[0]);
                    repmax[0] = std::min(rep_cn_[0], jtaux + rep_cn_[0]);
                    for (int jtauy = -rep_cn_[1]; jtauy <= rep_cn_[1]; jtauy++)
                    {
                        repmin[1] = std::max(-rep_cn_[1], jtauy - rep_cn_[1]);
                        repmax[1] = std::min(rep_cn_[1], jtauy + rep_cn_[1]);
                        for (int jtauz = -rep_cn_[2]; jtauz <= rep_cn_[2]; jtauz++)
                        {
                            repmin[2] = std::max(-rep_cn_[2], jtauz - rep_cn_[2]);
                            repmax[2] = std::min(rep_cn_[2], jtauz + rep_cn_[2]);
                            jtau = static_cast<double>(jtaux) * lat_[0] + static_cast<double>(jtauy) * lat_[1]
                                   + static_cast<double>(jtauz) * lat_[2];
                            rij2 = (ijvec + jtau).norm2();
                            if (rij2 > para_.cn_thr2())
                                continue;
                            rr0ij = std::sqrt(rij2) / r0ij;

                            for (int ktaux = repmin[0]; ktaux <= repmax[0]; ktaux++)
                                for (int ktauy = repmin[1]; ktauy <= repmax[1]; ktauy++)
                                    for (int ktauz = repmin[2]; ktauz <= repmax[2]; ktauz++)
                                    {
                                        ktau = static_cast<double>(ktaux) * lat_[0]
                                               + static_cast<double>(ktauy) * lat_[1]
                                               + static_cast<double>(ktauz) * lat_[2];
                                        rik2 = (ikvec + ktau).norm2();
                                        if (rik2 > para_.cn_thr2())
                                            continue;
                                        rr0ik = std::sqrt(rik2) / r0ik;

                                        rjk2 = (jkvec + ktau - jtau).norm2();
                                        if (rjk2 > para_.cn_thr2())
                                            continue;
                                        rr0jk = std::sqrt(rjk2) / r0jk;

                                        geomean = std::pow(rr0ij * rr0ik * rr0jk, 1.0 / 3.0);
                                        fdamp = 1.0 / (1.0 + 6.0 * std::pow(sr9 * geomean, alp9));
                                        tmp1 = (rij2 + rjk2 - rik2);
                                        tmp2 = (rij2 + rik2 - rjk2);
                                        tmp3 = (rik2 + rjk2 - rij2);
                                        tmp4 = rij2 * rjk2 * rik2;

                                        ang = (0.375 * tmp1 * tmp2 * tmp3 / tmp4 + 1.0) / std::pow(tmp4, 1.5);

                                        eabc += ang * c9 * fdamp;
                                    } // end ktau
                        } // end jtauz
                    } // end jtauy
                } // end jtaux
            } // end kat
        } // end jat
    // end iat

    for (int iat = 1; iat != ucell_.nat; iat++)
    {
        int jat = iat;
        ij = lin(iat, jat);
        ijvec.set(0, 0, 0);
        r0ij = para_.r0ab()[iz_[jat]][iz_[iat]];
        for (int kat = 0; kat != iat; kat++)
        {
            jk = lin(jat, kat);
            ik = jk;
            ikvec = xyz_[kat] - xyz_[iat];
            jkvec = ikvec;
            c9 = -cc6ab[ij] * cc6ab[ik] * cc6ab[jk];

            r0ik = para_.r0ab()[iz_[kat]][iz_[iat]];
            r0jk = para_.r0ab()[iz_[kat]][iz_[jat]];
            for (int jtaux = -rep_cn_[0]; jtaux <= rep_cn_[0]; jtaux++)
            {
                repmin[0] = std::max(-rep_cn_[0], jtaux - rep_cn_[0]);
                repmax[0] = std::min(rep_cn_[0], jtaux + rep_cn_[0]);
                for (int jtauy = -rep_cn_[1]; jtauy <= rep_cn_[1]; jtauy++)
                {
                    repmin[1] = std::max(-rep_cn_[1], jtauy - rep_cn_[1]);
                    repmax[1] = std::min(rep_cn_[1], jtauy + rep_cn_[1]);
                    for (int jtauz = -rep_cn_[2]; jtauz <= rep_cn_[2]; jtauz++)
                    {
                        repmin[2] = std::max(-rep_cn_[2], jtauz - rep_cn_[2]);
                        repmax[2] = std::min(rep_cn_[2], jtauz + rep_cn_[2]);
                        if (jtaux == 0 && jtauy == 0 && jtauz == 0)
                            continue;
                        jtau = static_cast<double>(jtaux) * lat_[0] + static_cast<double>(jtauy) * lat_[1]
                               + static_cast<double>(jtauz) * lat_[2];
                        rij2 = (ijvec + jtau).norm2();
                        if (rij2 > para_.cn_thr2())
                            continue;
                        rr0ij = std::sqrt(rij2) / r0ij;

                        for (int ktaux = repmin[0]; ktaux <= repmax[0]; ktaux++)
                            for (int ktauy = repmin[1]; ktauy <= repmax[1]; ktauy++)
                                for (int ktauz = repmin[2]; ktauz <= repmax[2]; ktauz++)
                                {
                                    ktau = static_cast<double>(ktaux) * lat_[0] + static_cast<double>(ktauy) * lat_[1]
                                           + static_cast<double>(ktauz) * lat_[2];
                                    rik2 = (ikvec + ktau).norm2();
                                    if (rik2 > para_.cn_thr2())
                                        continue;
                                    rr0ik = std::sqrt(rik2) / r0ik;

                                    rjk2 = (jkvec + ktau - jtau).norm2();
                                    if (rjk2 > para_.cn_thr2())
                                        continue;
                                    rr0jk = std::sqrt(rjk2) / r0jk;

                                    geomean = std::pow(rr0ij * rr0ik * rr0jk, 1.0 / 3.0);
                                    fdamp = 1.0 / (1.0 + 6.0 * std::pow(sr9 * geomean, alp9));
                                    tmp1 = (rij2 + rjk2 - rik2);
                                    tmp2 = (rij2 + rik2 - rjk2);
                                    tmp3 = (rik2 + rjk2 - rij2);
                                    tmp4 = rij2 * rjk2 * rik2;

                                    ang = (0.375 * tmp1 * tmp2 * tmp3 / tmp4 + 1.0) / std::pow(tmp4, 1.5);

                                    eabc += ang * c9 * fdamp / 2.0;
                                } // end ktau
                    } // end jtauz
                } // end jtauy
            } // end jtaux
        } // end kat
    } // end iat

    for (int iat = 1; iat != ucell_.nat; iat++)
        for (int jat = 0; jat != iat; jat++)
        {
            int kat = jat;
            ij = lin(iat, jat);
            jk = lin(jat, kat);
            ik = ij;
            ikvec = xyz_[kat] - xyz_[iat];
            ijvec = ikvec;
            jkvec.set(0, 0, 0);
            c9 = -cc6ab[ij] * cc6ab[ik] * cc6ab[jk];

            r0ij = para_.r0ab()[iz_[jat]][iz_[iat]];
            r0ik = r0ij;
            r0jk = para_.r0ab()[iz_[kat]][iz_[jat]];

            for (int jtaux = -rep_cn_[0]; jtaux <= rep_cn_[0]; jtaux++)
            {
                repmin[0] = std::max(-rep_cn_[0], jtaux - rep_cn_[0]);
                repmax[0] = std::min(rep_cn_[0], jtaux + rep_cn_[0]);
                for (int jtauy = -rep_cn_[1]; jtauy <= rep_cn_[1]; jtauy++)
                {
                    repmin[1] = std::max(-rep_cn_[1], jtauy - rep_cn_[1]);
                    repmax[1] = std::min(rep_cn_[1], jtauy + rep_cn_[1]);
                    for (int jtauz = -rep_cn_[2]; jtauz <= rep_cn_[2]; jtauz++)
                    {
                        repmin[2] = std::max(-rep_cn_[2], jtauz - rep_cn_[2]);
                        repmax[2] = std::min(rep_cn_[2], jtauz + rep_cn_[2]);
                        jtau = static_cast<double>(jtaux) * lat_[0] + static_cast<double>(jtauy) * lat_[1]
                               + static_cast<double>(jtauz) * lat_[2];
                        rij2 = (ijvec + jtau).norm2();
                        if (rij2 > para_.cn_thr2())
                            continue;
                        rr0ij = std::sqrt(rij2) / r0ij;

                        for (int ktaux = repmin[0]; ktaux <= repmax[0]; ktaux++)
                            for (int ktauy = repmin[1]; ktauy <= repmax[1]; ktauy++)
                                for (int ktauz = repmin[2]; ktauz <= repmax[2]; ktauz++)
                                {
                                    if (jtaux == ktaux && jtauy == ktauy && jtauz == ktauz)
                                        continue;
                                    ktau = static_cast<double>(ktaux) * lat_[0] + static_cast<double>(ktauy) * lat_[1]
                                           + static_cast<double>(ktauz) * lat_[2];
                                    rik2 = (ikvec + ktau).norm2();
                                    if (rik2 > para_.cn_thr2())
                                        continue;
                                    rr0ik = std::sqrt(rik2) / r0ik;

                                    rjk2 = (jkvec + ktau - jtau).norm2();
                                    if (rjk2 > para_.cn_thr2())
                                        continue;
                                    rr0jk = std::sqrt(rjk2) / r0jk;

                                    geomean = std::pow(rr0ij * rr0ik * rr0jk, 1.0 / 3.0);
                                    fdamp = 1.0 / (1.0 + 6.0 * std::pow(sr9 * geomean, alp9));
                                    tmp1 = (rij2 + rjk2 - rik2);
                                    tmp2 = (rij2 + rik2 - rjk2);
                                    tmp3 = (rik2 + rjk2 - rij2);
                                    tmp4 = rij2 * rjk2 * rik2;

                                    ang = (0.375 * tmp1 * tmp2 * tmp3 / tmp4 + 1.0) / std::pow(tmp4, 1.5);

                                    eabc += ang * c9 * fdamp / 2.0;
                                } // end ktau
                    } // end jtauz
                } // end jtauy
            } // end jtaux
        } // end jat
    // end iat

    for (int iat = 0; iat != ucell_.nat; iat++)
    {
        int jat = iat;
        int kat = iat;
        ijvec.set(0, 0, 0);
        ij = lin(iat, iat);
        ik = ij;
        jk = ij;
        ikvec = ijvec;
        jkvec = ikvec;
        c9 = -cc6ab[ij] * cc6ab[ik] * cc6ab[jk];

        r0ij = para_.r0ab()[iz_[iat]][iz_[iat]];
        r0ik = r0ij;
        r0jk = r0ij;

        for (int jtaux = -rep_cn_[0]; jtaux <= rep_cn_[0]; jtaux++)
        {
            repmin[0] = std::max(-rep_cn_[0], jtaux - rep_cn_[0]);
            repmax[0] = std::min(rep_cn_[0], jtaux + rep_cn_[0]);
            for (int jtauy = -rep_cn_[1]; jtauy <= rep_cn_[1]; jtauy++)
            {
                repmin[1] = std::max(-rep_cn_[1], jtauy - rep_cn_[1]);
                repmax[1] = std::min(rep_cn_[1], jtauy + rep_cn_[1]);
                for (int jtauz = -rep_cn_[2]; jtauz <= rep_cn_[2]; jtauz++)
                {
                    repmin[2] = std::max(-rep_cn_[2], jtauz - rep_cn_[2]);
                    repmax[2] = std::min(rep_cn_[2], jtauz + rep_cn_[2]);
                    jtau = static_cast<double>(jtaux) * lat_[0] + static_cast<double>(jtauy) * lat_[1]
                           + static_cast<double>(jtauz) * lat_[2];
                    if (jtaux == 0 && jtauy == 0 && jtauz == 0)
                        continue;
                    rij2 = jtau.norm2();
                    if (rij2 > para_.cn_thr2())
                        continue;
                    rr0ij = std::sqrt(rij2) / r0ij;

                    for (int ktaux = repmin[0]; ktaux <= repmax[0]; ktaux++)
                        for (int ktauy = repmin[1]; ktauy <= repmax[1]; ktauy++)
                            for (int ktauz = repmin[2]; ktauz <= repmax[2]; ktauz++)
                            {
                                if (ktaux == 0 && ktauy == 0 && ktauz == 0)
                                    continue;
                                if (jtaux == ktaux && jtauy == ktauy && jtauz == ktauz)
                                    continue;
                                ktau = static_cast<double>(ktaux) * lat_[0] + static_cast<double>(ktauy) * lat_[1]
                                       + static_cast<double>(ktauz) * lat_[2];
                                rik2 = ktau.norm2();
                                if (rik2 > para_.cn_thr2())
                                    continue;
                                rr0ik = std::sqrt(rik2) / r0ik;

                                rjk2 = (jkvec + ktau - jtau).norm2();
                                if (rjk2 > para_.cn_thr2())
                                    continue;
                                rr0jk = std::sqrt(rjk2) / r0jk;

                                geomean = std::pow(rr0ij * rr0ik * rr0jk, 1.0 / 3.0);
                                fdamp = 1.0 / (1.0 + 6.0 * std::pow(sr9 * geomean, alp9));
                                tmp1 = (rij2 + rjk2 - rik2);
                                tmp2 = (rij2 + rik2 - rjk2);
                                tmp3 = (rik2 + rjk2 - rij2);
                                tmp4 = rij2 * rjk2 * rik2;

                                ang = (0.375 * tmp1 * tmp2 * tmp3 / tmp4 + 1.0) / std::pow(tmp4, 1.5);

                                eabc += ang * c9 * fdamp / 6.0;
                            } // end ktau
                } // end jtauz
            } // end jtauy
        } // end jtaux
    } // end iat
}

void Vdwd3::get_dc6_dcnij(int mxci, int mxcj, double cni, double cnj, int izi, int izj,
                          int iat, int jat, double &c6check, double &dc6i, double &dc6j)
{
    double r_save = 9999.0, c6mem = -1e99, zaehler = 0.0, nenner = 0.0;
    double dzaehler_i = 0.0, dnenner_i = 0.0, dzaehler_j = 0.0, dnenner_j = 0.0;
    double c6ref = 0.0, cn_refi = 0.0, cn_refj = 0.0, r = 0.0, expterm = 0.0, term = 0.0;
    for (size_t a = 0; a != mxci; a++)
        for (size_t b = 0; b != mxcj; b++)
        {
            c6ref = para_.c6ab()[0][b][a][izj][izi];
            if (c6ref > 0)
            {
                cn_refi = para_.c6ab()[1][b][a][izj][izi];
                cn_refj = para_.c6ab()[2][b][a][izj][izi];
                r = (cn_refi - cni) * (cn_refi - cni) + (cn_refj - cnj) * (cn_refj - cnj);
                if (r < r_save)
                {
                    r_save = r;
                    c6mem = c6ref;
                }
                expterm = exp(para_.k3() * r);
                zaehler += c6ref * expterm;
                nenner += expterm;
                expterm *= 2.0 * para_.k3();
                term = expterm * (cni - cn_refi);
                dzaehler_i += c6ref * term;
                dnenner_i += term;

                term = expterm * (cnj - cn_refj);
                dzaehler_j += c6ref * term;
                dnenner_j += term;
            }
        }

    if (nenner > 1e-99)
    {
        c6check = zaehler / nenner;
        dc6i = ((dzaehler_i * nenner) - (dnenner_i * zaehler)) / (nenner * nenner);
        dc6j = ((dzaehler_j * nenner) - (dnenner_j * zaehler)) / (nenner * nenner);
    }
    else
    {
        c6check = c6mem;
        dc6i = 0.0;
        dc6j = 0.0;
    }
}

void Vdwd3::pbc_gdisp(std::vector<ModuleBase::Vector3<double>> &g, ModuleBase::matrix &smearing_sigma)
{
    std::vector<double> c6save(ucell_.nat * (ucell_.nat + 1)), dc6_rest_sum(ucell_.nat * (ucell_.nat + 1) / 2),
        dc6i(ucell_.nat), cn(ucell_.nat);
    pbc_ncoord(cn);
    std::vector<std::vector<double>> dc6ij(ucell_.nat, std::vector<double>(ucell_.nat));
    double c6 = 0.0, dc6iji = 0.0, dc6ijj = 0.0;
    double r = 0.0, r0 = 0.0, r2 = 0.0, r6 = 0.0, r7 = 0.0, r8 = 0.0, r9 = 0.0;
    double r42 = 0.0, rcovij = 0.0, t6 = 0.0, t8 = 0.0, dc6_rest = 0.0;
    int linii = 0, linij = 0;
    ModuleBase::Vector3<double> tau;
    std::vector<std::vector<std::vector<std::vector<double>>>> drij(
        ucell_.nat * (ucell_.nat + 1) / 2,
        std::vector<std::vector<std::vector<double>>>(
            2 * rep_vdw_[0] + 1,
            std::vector<std::vector<double>>(2 * rep_vdw_[1] + 1, std::vector<double>(2 * rep_vdw_[2] + 1))));
    if (para_.version() == "d3_0")
    {
        double damp6 = 0.0, damp8 = 0.0;
        for (int iat = 0; iat != ucell_.nat; iat++)
        {
            get_dc6_dcnij(para_.mxc()[iz_[iat]], para_.mxc()[iz_[iat]], cn[iat], cn[iat],
                          iz_[iat], iz_[iat], iat, iat, c6, dc6iji, dc6ijj);

            linii = lin(iat, iat);
            c6save[linii] = c6;
            dc6ij[iat][iat] = dc6iji;
            r0 = para_.r0ab()[iz_[iat]][iz_[iat]];
            r42 = para_.r2r4()[iz_[iat]] * para_.r2r4()[iz_[iat]];
            rcovij = para_.rcov()[iz_[iat]] + para_.rcov()[iz_[iat]];

            for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                    for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                    {
                        tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                              + static_cast<double>(tauz) * lat_[2];

                        // first dE/d(tau)
                        r2 = tau.norm2();
                        if (r2 > 0.1 && r2 < para_.rthr2())
                        {
                            r = std::sqrt(r2);
                            r6 = std::pow(r2, 3);
                            r7 = r6 * r;
                            r8 = r6 * r2;
                            r9 = r8 * r;

                            t6 = std::pow(r / (para_.rs6() * r0), -para_.alp6());
                            damp6 = 1.0 / (1.0 + 6.0 * t6);
                            t8 = std::pow(r / (para_.rs18() * r0), -para_.alp8());
                            damp8 = 1.0 / (1.0 + 6.0 * t8);

                            // d(r^(-6))/d(tau)
                            drij[linii][taux + rep_vdw_[0]][tauy + rep_vdw_[1]][tauz + rep_vdw_[2]]
                                += (-para_.s6() * (6.0 / (r7)*c6 * damp6)
                                    - para_.s18() * (24.0 / (r9)*c6 * r42 * damp8))
                                   * 0.5;
                            // d(f_dmp)/d(tau)
                            drij[linii][taux + rep_vdw_[0]][tauy + rep_vdw_[1]][tauz + rep_vdw_[2]]
                                += (para_.s6() * c6 / r7 * 6.0 * para_.alp6() * t6 * damp6 * damp6
                                    + para_.s18() * c6 * r42 / r9 * 18.0 * para_.alp8() * t8 * damp8 * damp8)
                                   * 0.5;

                            dc6_rest = (para_.s6() / r6 * damp6 + 3.0 * para_.s18() * r42 / r8 * damp8) * 0.5;
                            dc6i[iat] += dc6_rest * (dc6iji + dc6ijj);
                            dc6_rest_sum[linii] += dc6_rest;
                        }
                    } // end tau
            for (int jat = 0; jat != iat; jat++)
            {
                get_dc6_dcnij(para_.mxc()[iz_[iat]], para_.mxc()[iz_[jat]], cn[iat], cn[jat],
                              iz_[iat], iz_[jat], iat, jat, c6, dc6iji, dc6ijj);

                linij = lin(iat, jat);
                c6save[linij] = c6;
                r0 = para_.r0ab()[iz_[iat]][iz_[jat]];
                r42 = para_.r2r4()[iz_[iat]] * para_.r2r4()[iz_[jat]];
                rcovij = para_.rcov()[iz_[iat]] + para_.rcov()[iz_[jat]];
                dc6ij[jat][iat] = dc6iji;
                dc6ij[iat][jat] = dc6ijj;
                for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                    for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                        for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                        {
                            tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                                  + static_cast<double>(tauz) * lat_[2];
                            r2 = (xyz_[jat] - xyz_[iat] + tau).norm2();
                            if (r2 > para_.rthr2())
                                continue;

                            r = std::sqrt(r2);
                            r6 = std::pow(r2, 3);
                            r7 = r6 * r;
                            r8 = r6 * r2;
                            r9 = r8 * r;

                            t6 = std::pow(r / (para_.rs6() * r0), -para_.alp6());
                            damp6 = 1.0 / (1.0 + 6.0 * t6);
                            t8 = std::pow(r / (para_.rs18() * r0), -para_.alp8());
                            damp8 = 1.0 / (1.0 + 6.0 * t8);

                            // d(r^(-6))/d(r_ij)
                            drij[linij][taux + rep_vdw_[0]][tauy + rep_vdw_[1]][tauz + rep_vdw_[2]]
                                += -para_.s6() * (6.0 / (r7)*c6 * damp6) - para_.s18() * (24.0 / (r9)*c6 * r42 * damp8);
                            // d(f_dmp)/d(r_ij)
                            drij[linij][taux + rep_vdw_[0]][tauy + rep_vdw_[1]][tauz + rep_vdw_[2]]
                                += para_.s6() * c6 / r7 * 6.0 * para_.alp6() * t6 * damp6 * damp6
                                   + para_.s18() * c6 * r42 / r9 * 18.0 * para_.alp8() * t8 * damp8 * damp8;

                            dc6_rest = para_.s6() / r6 * damp6 + 3.0 * para_.s18() * r42 / r8 * damp8;
                            dc6i[iat] += dc6_rest * dc6iji;
                            dc6i[jat] += dc6_rest * dc6ijj;
                            dc6_rest_sum[linij] += dc6_rest;
                        } // end tau
            } // end jat
        } // end iat
    } // end d3_0
    else if (para_.version() == "d3_bj")
    {
        double r4;
        for (int iat = 0; iat != ucell_.nat; iat++)
        {
            get_dc6_dcnij(para_.mxc()[iz_[iat]], para_.mxc()[iz_[iat]], cn[iat], cn[iat],
                          iz_[iat], iz_[iat], iat, iat, c6, dc6iji, dc6ijj);

            linii = lin(iat, iat);
            c6save[linii] = c6;
            dc6ij[iat][iat] = dc6iji;
            r42 = para_.r2r4()[iz_[iat]] * para_.r2r4()[iz_[iat]];
            r0 = para_.rs6() * std::sqrt(3.0 * r42) + para_.rs18();
            rcovij = para_.rcov()[iz_[iat]] + para_.rcov()[iz_[iat]];

            for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                    for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                    {
                        tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                              + static_cast<double>(tauz) * lat_[2];

                        // first dE/d(tau)
                        r2 = tau.norm2();
                        if (r2 > 0.1 && r2 < para_.rthr2())
                        {
                            r = std::sqrt(r2);
                            r4 = r2 * r2;
                            r6 = std::pow(r2, 3);
                            r7 = r6 * r;
                            r8 = r6 * r2;
                            r9 = r8 * r;

                            t6 = r6 + std::pow(r0, 6);
                            t8 = r8 + std::pow(r0, 8);

                            // d(1/r^(-6)+r0^6)/d(r)
                            drij[linii][taux + rep_vdw_[0]][tauy + rep_vdw_[1]][tauz + rep_vdw_[2]]
                                += -para_.s6() * c6 * 6.0 * r4 * r / (t6 * t6) * 0.5
                                   - para_.s18() * c6 * 24.0 * r42 * r7 / (t8 * t8) * 0.5;

                            dc6_rest = (para_.s6() / t6 + 3.0 * para_.s18() * r42 / t8) * 0.5;
                            dc6i[iat] += dc6_rest * (dc6iji + dc6ijj);
                            dc6_rest_sum[linii] += dc6_rest;
                        }
                    } // end tau
            for (int jat = 0; jat != iat; jat++)
            {
                get_dc6_dcnij(para_.mxc()[iz_[iat]], para_.mxc()[iz_[jat]], cn[iat], cn[jat],
                              iz_[iat], iz_[jat], iat, jat, c6, dc6iji, dc6ijj);

                linij = lin(iat, jat);
                c6save[linij] = c6;
                r42 = para_.r2r4()[iz_[iat]] * para_.r2r4()[iz_[jat]];
                r0 = para_.rs6() * std::sqrt(3.0 * r42) + para_.rs18();
                rcovij = para_.rcov()[iz_[iat]] + para_.rcov()[iz_[jat]];
                dc6ij[jat][iat] = dc6iji;
                dc6ij[iat][jat] = dc6ijj;
                for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                    for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                        for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                        {
                            tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                                  + static_cast<double>(tauz) * lat_[2];
                            r2 = (xyz_[jat] - xyz_[iat] + tau).norm2();
                            if (r2 > para_.rthr2())
                                continue;

                            r = std::sqrt(r2);
                            r4 = r2 * r2;
                            r6 = std::pow(r2, 3);
                            r7 = r6 * r;
                            r8 = r6 * r2;
                            r9 = r8 * r;

                            t6 = r6 + std::pow(r0, 6);
                            t8 = r8 + std::pow(r0, 8);

                            drij[linij][taux + rep_vdw_[0]][tauy + rep_vdw_[1]][tauz + rep_vdw_[2]]
                                += -para_.s6() * c6 * 6.0 * r4 * r / (t6 * t6)
                                   - para_.s18() * c6 * 24.0 * r42 * r7 / (t8 * t8);

                            dc6_rest = para_.s6() / t6 + 3.0 * para_.s18() * r42 / t8;
                            dc6i[iat] += dc6_rest * dc6iji;
                            dc6i[jat] += dc6_rest * dc6ijj;
                            dc6_rest_sum[linij] += dc6_rest;
                        } // end tau
            } // end jat
        } // end iat
    } // end d3_bj

    if (para_.abc())
    {
        ModuleBase::Vector3<double> ijvec, ikvec, jkvec, jtau, ktau;
        std::vector<int> repmin(3), repmax(3);
        double sr9 = 0.75, alp9 = -16.0;
        double linik, linjk, rij2, rik2, rjk2, rr0ij, rr0ik, rr0jk, geomean2, geomean, geomean3, r0av, r;
        double c6ij, c6ik, c6jk, c9, damp9, ang, dfdmp, dang, tmp1, dc9;
        for (int iat = 2; iat < ucell_.nat; iat++)
        {
            for (int jat = 1; jat != iat; jat++)
            {
                linij = lin(iat, jat);
                ijvec = xyz_[jat] - xyz_[iat];

                c6ij = c6save[linij];
                for (int kat = 0; kat != jat; kat++)
                {
                    linik = lin(iat, kat);
                    linjk = lin(jat, kat);
                    ikvec = xyz_[kat] - xyz_[iat];
                    jkvec = xyz_[kat] - xyz_[jat];

                    c6ik = c6save[linik];
                    c6jk = c6save[linjk];
                    c9 = -1.0 * std::sqrt(c6ij * c6ik * c6jk);

                    for (int jtaux = -rep_cn_[0]; jtaux <= rep_cn_[0]; jtaux++)
                    {
                        repmin[0] = std::max(-rep_cn_[0], jtaux - rep_cn_[0]);
                        repmax[0] = std::min(rep_cn_[0], jtaux + rep_cn_[0]);
                        for (int jtauy = -rep_cn_[1]; jtauy <= rep_cn_[1]; jtauy++)
                        {
                            repmin[1] = std::max(-rep_cn_[1], jtauy - rep_cn_[1]);
                            repmax[1] = std::min(rep_cn_[1], jtauy + rep_cn_[1]);
                            for (int jtauz = -rep_cn_[2]; jtauz <= rep_cn_[2]; jtauz++)
                            {
                                repmin[2] = std::max(-rep_cn_[2], jtauz - rep_cn_[2]);
                                repmax[2] = std::min(rep_cn_[2], jtauz + rep_cn_[2]);
                                jtau = static_cast<double>(jtaux) * lat_[0] + static_cast<double>(jtauy) * lat_[1]
                                       + static_cast<double>(jtauz) * lat_[2];
                                rij2 = (ijvec + jtau).norm2();
                                if (rij2 > para_.cn_thr2())
                                    continue;
                                rr0ij = std::sqrt(rij2) / para_.r0ab()[iz_[jat]][iz_[iat]];

                                for (int ktaux = repmin[0]; ktaux <= repmax[0]; ktaux++)
                                    for (int ktauy = repmin[1]; ktauy <= repmax[1]; ktauy++)
                                        for (int ktauz = repmin[2]; ktauz <= repmax[2]; ktauz++)
                                        {
                                            ktau = static_cast<double>(ktaux) * lat_[0]
                                                   + static_cast<double>(ktauy) * lat_[1]
                                                   + static_cast<double>(ktauz) * lat_[2];
                                            rik2 = (ikvec + ktau).norm2();
                                            if (rik2 > para_.cn_thr2())
                                                continue;

                                            rjk2 = (jkvec + ktau - jtau).norm2();
                                            if (rjk2 > para_.cn_thr2())
                                                continue;
                                            rr0ik = std::sqrt(rik2) / para_.r0ab()[iz_[kat]][iz_[iat]];
                                            rr0jk = std::sqrt(rjk2) / para_.r0ab()[iz_[kat]][iz_[jat]];

                                            geomean2 = rij2 * rjk2 * rik2;
                                            r0av = std::pow(rr0ij * rr0ik * rr0jk, 1.0 / 3.0);
                                            damp9 = 1.0 / (1.0 + 6.0 * std::pow(sr9 * r0av, alp9));
                                            geomean = std::sqrt(geomean2);
                                            geomean3 = geomean * geomean2;
                                            ang = 0.375 * (rij2 + rjk2 - rik2) * (rij2 - rjk2 + rik2)
                                                      * (-rij2 + rjk2 + rik2) / (geomean3 * geomean2)
                                                  + 1.0 / geomean3;
                                            dc6_rest = ang * damp9;
                                            dfdmp = 2.0 * alp9 * std::pow(0.75 * r0av, alp9) * damp9 * damp9;

                                            r = std::sqrt(rij2);
                                            dang = -0.375
                                                   * (std::pow(rij2, 3) + std::pow(rij2, 2) * (rjk2 + rik2)
                                                      + rij2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rik2
                                                                + 3.0 * std::pow(rik2, 2))
                                                      - 5.0 * std::pow(rjk2 - rik2, 2) * (rjk2 + rik2))
                                                   / (r * geomean3 * geomean2);
                                            tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                            drij[linij][jtaux + rep_vdw_[0]][jtauy + rep_vdw_[1]][jtauz + rep_vdw_[2]]
                                                -= tmp1;

                                            r = std::sqrt(rik2);
                                            dang = -0.375
                                                   * (std::pow(rik2, 3) + std::pow(rik2, 2) * (rjk2 + rij2)
                                                      + rik2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rij2
                                                                + 3.0 * std::pow(rij2, 2))
                                                      - 5.0 * std::pow(rjk2 - rij2, 2) * (rjk2 + rij2))
                                                   / (r * geomean3 * geomean2);
                                            tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                            drij[linik][ktaux + rep_vdw_[0]][ktauy + rep_vdw_[1]][ktauz + rep_vdw_[2]]
                                                -= tmp1;

                                            r = std::sqrt(rjk2);
                                            dang = -0.375
                                                   * (std::pow(rjk2, 3) + std::pow(rjk2, 2) * (rik2 + rij2)
                                                      + rjk2 * (3.0 * std::pow(rik2, 2) + 2.0 * rik2 * rij2
                                                                + 3.0 * std::pow(rij2, 2))
                                                      - 5.0 * std::pow(rik2 - rij2, 2) * (rik2 + rij2))
                                                   / (r * geomean3 * geomean2);
                                            tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                            drij[linjk][ktaux - jtaux + rep_vdw_[0]][ktauy - jtauy + rep_vdw_[1]]
                                                [ktauz - jtauz + rep_vdw_[2]]
                                                -= tmp1;

                                            dc9 = (dc6ij[jat][iat] / c6ij + dc6ij[kat][iat] / c6ik) * c9 * 0.5;
                                            dc6i[iat] += dc6_rest * dc9;

                                            dc9 = (dc6ij[iat][jat] / c6ij + dc6ij[kat][jat] / c6jk) * c9 * 0.5;
                                            dc6i[jat] += dc6_rest * dc9;

                                            dc9 = (dc6ij[iat][kat] / c6ik + dc6ij[jat][kat] / c6jk) * c9 * 0.5;
                                            dc6i[kat] += dc6_rest * dc9;
                                        } // end ktau
                            } // end jtauz
                        } // end jtauy
                    } // end jtaux
                } // end kat
            } // end jat
        }
        for (int iat = 1; iat != ucell_.nat; iat++)
        {
            int jat = iat;
            linij = lin(iat, jat);
            ijvec.set(0, 0, 0);

            c6ij = c6save[linij];
            for (int kat = 0; kat != iat; kat++)
            {
                linjk = lin(jat, kat);
                linik = linjk;
                ikvec = xyz_[kat] - xyz_[iat];
                jkvec = ikvec;

                c6ik = c6save[linik];
                c6jk = c6ik;
                c9 = -1.0 * std::sqrt(c6ij * c6ik * c6jk);

                for (int jtaux = -rep_cn_[0]; jtaux <= rep_cn_[0]; jtaux++)
                {
                    repmin[0] = std::max(-rep_cn_[0], jtaux - rep_cn_[0]);
                    repmax[0] = std::min(rep_cn_[0], jtaux + rep_cn_[0]);
                    for (int jtauy = -rep_cn_[1]; jtauy <= rep_cn_[1]; jtauy++)
                    {
                        repmin[1] = std::max(-rep_cn_[1], jtauy - rep_cn_[1]);
                        repmax[1] = std::min(rep_cn_[1], jtauy + rep_cn_[1]);
                        for (int jtauz = -rep_cn_[2]; jtauz <= rep_cn_[2]; jtauz++)
                        {
                            repmin[2] = std::max(-rep_cn_[2], jtauz - rep_cn_[2]);
                            repmax[2] = std::min(rep_cn_[2], jtauz + rep_cn_[2]);
                            if (jtaux == 0 && jtauy == 0 && jtauz == 0)
                                continue;
                            jtau = static_cast<double>(jtaux) * lat_[0] + static_cast<double>(jtauy) * lat_[1]
                                   + static_cast<double>(jtauz) * lat_[2];
                            rij2 = jtau.norm2();
                            if (rij2 > para_.cn_thr2())
                                continue;
                            rr0ij = std::sqrt(rij2) / para_.r0ab()[iz_[jat]][iz_[iat]];

                            for (int ktaux = repmin[0]; ktaux <= repmax[0]; ktaux++)
                                for (int ktauy = repmin[1]; ktauy <= repmax[1]; ktauy++)
                                    for (int ktauz = repmin[2]; ktauz <= repmax[2]; ktauz++)
                                    {
                                        ktau = static_cast<double>(ktaux) * lat_[0]
                                               + static_cast<double>(ktauy) * lat_[1]
                                               + static_cast<double>(ktauz) * lat_[2];
                                        rik2 = (ikvec + ktau).norm2();
                                        if (rik2 > para_.cn_thr2())
                                            continue;

                                        rjk2 = (jkvec + ktau - jtau).norm2();
                                        if (rjk2 > para_.cn_thr2())
                                            continue;
                                        rr0ik = std::sqrt(rik2) / para_.r0ab()[iz_[kat]][iz_[iat]];
                                        rr0jk = std::sqrt(rjk2) / para_.r0ab()[iz_[kat]][iz_[jat]];

                                        geomean2 = rij2 * rjk2 * rik2;
                                        r0av = std::pow(rr0ij * rr0ik * rr0jk, 1.0 / 3.0);
                                        damp9 = 1.0 / (1.0 + 6.0 * std::pow(sr9 * r0av, alp9));
                                        geomean = std::sqrt(geomean2);
                                        geomean3 = geomean * geomean2;
                                        ang = 0.375 * (rij2 + rjk2 - rik2) * (rij2 - rjk2 + rik2)
                                                  * (-rij2 + rjk2 + rik2) / (geomean3 * geomean2)
                                              + 1.0 / geomean3;
                                        dc6_rest = ang * damp9 / 2.0;
                                        dfdmp = 2.0 * alp9 * std::pow(0.75 * r0av, alp9) * damp9 * damp9;

                                        r = std::sqrt(rij2);
                                        dang = -0.375
                                               * (std::pow(rij2, 3) + std::pow(rij2, 2) * (rjk2 + rik2)
                                                  + rij2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rik2
                                                            + 3.0 * std::pow(rik2, 2))
                                                  - 5.0 * std::pow(rjk2 - rik2, 2) * (rjk2 + rik2))
                                               / (r * geomean3 * geomean2);
                                        tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                        drij[linij][jtaux + rep_vdw_[0]][jtauy + rep_vdw_[1]][jtauz + rep_vdw_[2]]
                                            -= tmp1 / 2.0;

                                        r = std::sqrt(rik2);
                                        dang = -0.375
                                               * (std::pow(rik2, 3) + std::pow(rik2, 2) * (rjk2 + rij2)
                                                  + rik2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rij2
                                                            + 3.0 * std::pow(rij2, 2))
                                                  - 5.0 * std::pow(rjk2 - rij2, 2) * (rjk2 + rij2))
                                               / (r * geomean3 * geomean2);
                                        tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                        drij[linik][ktaux + rep_vdw_[0]][ktauy + rep_vdw_[1]][ktauz + rep_vdw_[2]]
                                            -= tmp1 / 2.0;

                                        r = std::sqrt(rjk2);
                                        dang = -0.375
                                               * (std::pow(rjk2, 3) + std::pow(rjk2, 2) * (rik2 + rij2)
                                                  + rjk2 * (3.0 * std::pow(rik2, 2) + 2.0 * rik2 * rij2
                                                            + 3.0 * std::pow(rij2, 2))
                                                  - 5.0 * std::pow(rik2 - rij2, 2) * (rik2 + rij2))
                                               / (r * geomean3 * geomean2);
                                        tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                        drij[linjk][ktaux - jtaux + rep_vdw_[0]][ktauy - jtauy + rep_vdw_[1]]
                                            [ktauz - jtauz + rep_vdw_[2]]
                                            -= tmp1 / 2.0;

                                        dc9 = (dc6ij[jat][iat] / c6ij + dc6ij[kat][iat] / c6ik) * c9 * 0.5;
                                        dc6i[iat] += dc6_rest * dc9;

                                        dc9 = (dc6ij[iat][jat] / c6ij + dc6ij[kat][jat] / c6jk) * c9 * 0.5;
                                        dc6i[jat] += dc6_rest * dc9;

                                        dc9 = (dc6ij[iat][kat] / c6ik + dc6ij[jat][kat] / c6jk) * c9 * 0.5;
                                        dc6i[kat] += dc6_rest * dc9;
                                    } // end ktau
                        } // end jtauz
                    } // end jtauy
                } // end jtaux
            } // end kat
        } // end iat

        for (int iat = 1; iat != ucell_.nat; iat++)
            for (int jat = 0; jat != iat; jat++)
            {
                int kat = jat;
                linij = lin(iat, jat);
                linjk = lin(jat, kat);
                linik = linij;
                ikvec = xyz_[kat] - xyz_[iat];
                ijvec = ikvec;
                jkvec.set(0, 0, 0);

                c6ij = c6save[linij];
                c6ik = c6ij;
                c6jk = c6save[linjk];
                c9 = -1.0 * std::sqrt(c6ij * c6ik * c6jk);

                for (int jtaux = -rep_cn_[0]; jtaux <= rep_cn_[0]; jtaux++)
                {
                    repmin[0] = std::max(-rep_cn_[0], jtaux - rep_cn_[0]);
                    repmax[0] = std::min(rep_cn_[0], jtaux + rep_cn_[0]);
                    for (int jtauy = -rep_cn_[1]; jtauy <= rep_cn_[1]; jtauy++)
                    {
                        repmin[1] = std::max(-rep_cn_[1], jtauy - rep_cn_[1]);
                        repmax[1] = std::min(rep_cn_[1], jtauy + rep_cn_[1]);
                        for (int jtauz = -rep_cn_[2]; jtauz <= rep_cn_[2]; jtauz++)
                        {
                            repmin[2] = std::max(-rep_cn_[2], jtauz - rep_cn_[2]);
                            repmax[2] = std::min(rep_cn_[2], jtauz + rep_cn_[2]);
                            jtau = static_cast<double>(jtaux) * lat_[0] + static_cast<double>(jtauy) * lat_[1]
                                   + static_cast<double>(jtauz) * lat_[2];
                            rij2 = (ijvec + jtau).norm2();
                            if (rij2 > para_.cn_thr2())
                                continue;
                            rr0ij = std::sqrt(rij2) / para_.r0ab()[iz_[jat]][iz_[iat]];

                            for (int ktaux = repmin[0]; ktaux <= repmax[0]; ktaux++)
                                for (int ktauy = repmin[1]; ktauy <= repmax[1]; ktauy++)
                                    for (int ktauz = repmin[2]; ktauz <= repmax[2]; ktauz++)
                                    {
                                        if (jtaux == ktaux && jtauy == ktauy && jtauz == ktauz)
                                            continue;
                                        ktau = static_cast<double>(ktaux) * lat_[0]
                                               + static_cast<double>(ktauy) * lat_[1]
                                               + static_cast<double>(ktauz) * lat_[2];
                                        rik2 = (ikvec + ktau).norm2();
                                        if (rik2 > para_.cn_thr2())
                                            continue;
                                        rr0ik = std::sqrt(rik2) / para_.r0ab()[iz_[kat]][iz_[iat]];

                                        rjk2 = (jkvec + ktau - jtau).norm2();
                                        if (rjk2 > para_.cn_thr2())
                                            continue;
                                        rr0jk = std::sqrt(rjk2) / para_.r0ab()[iz_[kat]][iz_[jat]];

                                        geomean2 = rij2 * rjk2 * rik2;
                                        r0av = std::pow(rr0ij * rr0ik * rr0jk, 1.0 / 3.0);
                                        damp9 = 1.0 / (1.0 + 6.0 * std::pow(sr9 * r0av, alp9));
                                        geomean = std::sqrt(geomean2);
                                        geomean3 = geomean * geomean2;
                                        ang = 0.375 * (rij2 + rjk2 - rik2) * (rij2 - rjk2 + rik2)
                                                  * (-rij2 + rjk2 + rik2) / (geomean3 * geomean2)
                                              + 1.0 / geomean3;
                                        dc6_rest = ang * damp9 / 2.0;
                                        dfdmp = 2.0 * alp9 * std::pow(0.75 * r0av, alp9) * damp9 * damp9;

                                        r = std::sqrt(rij2);
                                        dang = -0.375
                                               * (std::pow(rij2, 3) + std::pow(rij2, 2) * (rjk2 + rik2)
                                                  + rij2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rik2
                                                            + 3.0 * std::pow(rik2, 2))
                                                  - 5.0 * std::pow(rjk2 - rik2, 2) * (rjk2 + rik2))
                                               / (r * geomean3 * geomean2);
                                        tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                        drij[linij][jtaux + rep_vdw_[0]][jtauy + rep_vdw_[1]][jtauz + rep_vdw_[2]]
                                            -= tmp1 / 2.0;

                                        r = std::sqrt(rik2);
                                        dang = -0.375
                                               * (std::pow(rik2, 3) + std::pow(rik2, 2) * (rjk2 + rij2)
                                                  + rik2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rij2
                                                            + 3.0 * std::pow(rij2, 2))
                                                  - 5.0 * std::pow(rjk2 - rij2, 2) * (rjk2 + rij2))
                                               / (r * geomean3 * geomean2);
                                        tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                        drij[linik][ktaux + rep_vdw_[0]][ktauy + rep_vdw_[1]][ktauz + rep_vdw_[2]]
                                            -= tmp1 / 2.0;

                                        r = std::sqrt(rjk2);
                                        dang = -0.375
                                               * (std::pow(rjk2, 3) + std::pow(rjk2, 2) * (rik2 + rij2)
                                                  + rjk2 * (3.0 * std::pow(rik2, 2) + 2.0 * rik2 * rij2
                                                            + 3.0 * std::pow(rij2, 2))
                                                  - 5.0 * std::pow(rik2 - rij2, 2) * (rik2 + rij2))
                                               / (r * geomean3 * geomean2);
                                        tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                        drij[linjk][ktaux - jtaux + rep_vdw_[0]][ktauy - jtauy + rep_vdw_[1]]
                                            [ktauz - jtauz + rep_vdw_[2]]
                                            -= tmp1 / 2.0;

                                        dc9 = (dc6ij[jat][iat] / c6ij + dc6ij[kat][iat] / c6ik) * c9 * 0.5;
                                        dc6i[iat] += dc6_rest * dc9;

                                        dc9 = (dc6ij[iat][jat] / c6ij + dc6ij[kat][jat] / c6jk) * c9 * 0.5;
                                        dc6i[jat] += dc6_rest * dc9;

                                        dc9 = (dc6ij[iat][kat] / c6ik + dc6ij[jat][kat] / c6jk) * c9 * 0.5;
                                        dc6i[kat] += dc6_rest * dc9;
                                    } // end ktau
                        } // end jtauz
                    } // end jtauy
                } // end jtaux
            } // end jat
        // end iat

        for (int iat = 0; iat != ucell_.nat; iat++)
        {
            int jat = iat;
            int kat = iat;
            ijvec.set(0, 0, 0);
            linij = lin(iat, jat);
            linik = lin(iat, kat);
            linjk = lin(jat, kat);
            ikvec = ijvec;
            jkvec = ikvec;
            c6ij = c6save[linij];
            c6ik = c6ij;
            c6jk = c6ij;
            c9 = -1.0 * std::sqrt(c6ij * c6ik * c6jk);

            for (int jtaux = -rep_cn_[0]; jtaux <= rep_cn_[0]; jtaux++)
            {
                repmin[0] = std::max(-rep_cn_[0], jtaux - rep_cn_[0]);
                repmax[0] = std::min(rep_cn_[0], jtaux + rep_cn_[0]);
                for (int jtauy = -rep_cn_[1]; jtauy <= rep_cn_[1]; jtauy++)
                {
                    repmin[1] = std::max(-rep_cn_[1], jtauy - rep_cn_[1]);
                    repmax[1] = std::min(rep_cn_[1], jtauy + rep_cn_[1]);
                    for (int jtauz = -rep_cn_[2]; jtauz <= rep_cn_[2]; jtauz++)
                    {
                        repmin[2] = std::max(-rep_cn_[2], jtauz - rep_cn_[2]);
                        repmax[2] = std::min(rep_cn_[2], jtauz + rep_cn_[2]);
                        if (jtaux == 0 && jtauy == 0 && jtauz == 0)
                            continue;
                        jtau = static_cast<double>(jtaux) * lat_[0] + static_cast<double>(jtauy) * lat_[1]
                               + static_cast<double>(jtauz) * lat_[2];
                        rij2 = jtau.norm2();
                        if (rij2 > para_.cn_thr2())
                            continue;
                        rr0ij = std::sqrt(rij2) / para_.r0ab()[iz_[jat]][iz_[iat]];

                        for (int ktaux = repmin[0]; ktaux <= repmax[0]; ktaux++)
                            for (int ktauy = repmin[1]; ktauy <= repmax[1]; ktauy++)
                                for (int ktauz = repmin[2]; ktauz <= repmax[2]; ktauz++)
                                {
                                    if (ktaux == 0 && ktauy == 0 && ktauz == 0)
                                        continue;
                                    if (jtaux == ktaux && jtauy == ktauy && jtauz == ktauz)
                                        continue;
                                    ktau = static_cast<double>(ktaux) * lat_[0] + static_cast<double>(ktauy) * lat_[1]
                                           + static_cast<double>(ktauz) * lat_[2];
                                    rik2 = ktau.norm2();
                                    if (rik2 > para_.cn_thr2())
                                        continue;
                                    rr0ik = std::sqrt(rik2) / para_.r0ab()[iz_[kat]][iz_[iat]];

                                    rjk2 = (jkvec + ktau - jtau).norm2();
                                    if (rjk2 > para_.cn_thr2())
                                        continue;
                                    rr0jk = std::sqrt(rjk2) / para_.r0ab()[iz_[kat]][iz_[jat]];

                                    geomean2 = rij2 * rjk2 * rik2;
                                    r0av = std::pow(rr0ij * rr0ik * rr0jk, 1.0 / 3.0);
                                    damp9 = 1.0 / (1.0 + 6.0 * std::pow(sr9 * r0av, alp9));
                                    geomean = std::sqrt(geomean2);
                                    geomean3 = geomean * geomean2;
                                    ang = 0.375 * (rij2 + rjk2 - rik2) * (rij2 - rjk2 + rik2) * (-rij2 + rjk2 + rik2)
                                              / (geomean3 * geomean2)
                                          + 1.0 / geomean3;
                                    dc6_rest = ang * damp9 / 6.0;
                                    dfdmp = 2.0 * alp9 * std::pow(0.75 * r0av, alp9) * damp9 * damp9;

                                    r = std::sqrt(rij2);
                                    dang = -0.375
                                           * (std::pow(rij2, 3) + std::pow(rij2, 2) * (rjk2 + rik2)
                                              + rij2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rik2
                                                        + 3.0 * std::pow(rik2, 2))
                                              - 5.0 * std::pow(rjk2 - rik2, 2) * (rjk2 + rik2))
                                           / (r * geomean3 * geomean2);
                                    tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                    drij[linij][jtaux + rep_vdw_[0]][jtauy + rep_vdw_[1]][jtauz + rep_vdw_[2]]
                                        -= tmp1 / 6.0;

                                    r = std::sqrt(rik2);
                                    dang = -0.375
                                           * (std::pow(rik2, 3) + std::pow(rik2, 2) * (rjk2 + rij2)
                                              + rik2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rij2
                                                        + 3.0 * std::pow(rij2, 2))
                                              - 5.0 * std::pow(rjk2 - rij2, 2) * (rjk2 + rij2))
                                           / (r * geomean3 * geomean2);
                                    tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                    drij[linik][ktaux + rep_vdw_[0]][ktauy + rep_vdw_[1]][ktauz + rep_vdw_[2]]
                                        -= tmp1 / 6.0;

                                    r = std::sqrt(rjk2);
                                    dang = -0.375
                                           * (std::pow(rjk2, 3) + std::pow(rjk2, 2) * (rik2 + rij2)
                                              + rjk2 * (3.0 * std::pow(rik2, 2) + 2.0 * rik2 * rij2
                                                        + 3.0 * std::pow(rij2, 2))
                                              - 5.0 * std::pow(rik2 - rij2, 2) * (rik2 + rij2))
                                           / (r * geomean3 * geomean2);
                                    tmp1 = -dang * c9 * damp9 + dfdmp / r * c9 * ang;
                                    drij[linjk][ktaux - jtaux + rep_vdw_[0]][ktauy - jtauy + rep_vdw_[1]]
                                        [ktauz - jtauz + rep_vdw_[2]]
                                        -= tmp1 / 6.0;

                                    dc9 = (dc6ij[jat][iat] / c6ij + dc6ij[kat][iat] / c6ik) * c9 * 0.5;
                                    dc6i[iat] += dc6_rest * dc9;

                                    dc9 = (dc6ij[iat][jat] / c6ij + dc6ij[kat][jat] / c6jk) * c9 * 0.5;
                                    dc6i[jat] += dc6_rest * dc9;

                                    dc9 = (dc6ij[iat][kat] / c6ik + dc6ij[jat][kat] / c6jk) * c9 * 0.5;
                                    dc6i[kat] += dc6_rest * dc9;
                                } // end ktau
                    } // end jtauz
                } // end jtauy
            } // jtaux
        } // end iat
    }

    // dE/dr_ij * dr_ij/dxyz_i
    double expterm, dcnn, x1;
    ModuleBase::Vector3<double> rij, vec3;
    for (int iat = 1; iat != ucell_.nat; iat++)
        for (int jat = 0; jat != iat; jat++)
        {
            linij = lin(iat, jat);
            rcovij = para_.rcov()[iz_[iat]] + para_.rcov()[iz_[jat]];
            for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
                for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                    for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                    {
                        tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                              + static_cast<double>(tauz) * lat_[2];
                        rij = xyz_[jat] - xyz_[iat] + tau;
                        r2 = rij.norm2();
                        if (r2 > para_.rthr2() || r2 < 0.5)
                            continue;
                        r = std::sqrt(r2);
                        if (r2 < para_.cn_thr2())
                        {
                            expterm = exp(-para_.k1() * (rcovij / r - 1.0));
                            dcnn = -para_.k1() * rcovij * expterm / (r2 * (expterm + 1.0) * (expterm + 1.0));
                        }
                        else
                            dcnn = 0.0;
                        x1 = drij[linij][taux + rep_vdw_[0]][tauy + rep_vdw_[1]][tauz + rep_vdw_[2]]
                             + dcnn * (dc6i[iat] + dc6i[jat]);
                        vec3 = x1 * rij / r;
                        g[iat] += vec3;
                        g[jat] -= vec3;

                        std::vector<double> vec = {vec3.x, vec3.y, vec3.z};
                        std::vector<double> rij_vec = {rij.x, rij.y, rij.z};
                        for (size_t i = 0; i != 3; i++)
                            for (size_t j = 0; j != 3; j++)
                            {
                                smearing_sigma(i, j) += vec[j] * rij_vec[i];
                            }
                    } // end tau
        } // end iat, jat
    for (int iat = 0; iat != ucell_.nat; iat++)
    {
        linii = lin(iat, iat);
        rcovij = para_.rcov()[iz_[iat]] + para_.rcov()[iz_[iat]];
        for (int taux = -rep_vdw_[0]; taux <= rep_vdw_[0]; taux++)
            for (int tauy = -rep_vdw_[1]; tauy <= rep_vdw_[1]; tauy++)
                for (int tauz = -rep_vdw_[2]; tauz <= rep_vdw_[2]; tauz++)
                {
                    if (taux == 0 && tauy == 0 && tauz == 0)
                        continue;
                    tau = static_cast<double>(taux) * lat_[0] + static_cast<double>(tauy) * lat_[1]
                          + static_cast<double>(tauz) * lat_[2];
                    r2 = tau.norm2();
                    r = std::sqrt(r2);
                    if (r2 < para_.cn_thr2())
                    {
                        expterm = exp(-para_.k1() * (rcovij / r - 1.0));
                        dcnn = -para_.k1() * rcovij * expterm / (r2 * (expterm + 1.0) * (expterm + 1.0));
                    }
                    else
                        dcnn = 0.0;
                    x1 = drij[linii][taux + rep_vdw_[0]][tauy + rep_vdw_[1]][tauz + rep_vdw_[2]] + dcnn * dc6i[iat];

                    vec3 = x1 * tau / r;
                    std::vector<double> vec = {vec3.x, vec3.y, vec3.z};
                    std::vector<double> tau_vec = {tau.x, tau.y, tau.z};
                    for (size_t i = 0; i != 3; i++)
                        for (size_t j = 0; j != 3; j++)
                        {
                            smearing_sigma(i, j) += vec[j] * tau_vec[i];
                        }
                } // end tau
    } // end iat
}

} // namespace vdw
