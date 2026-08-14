/**
 * @file reciprocal_grid.cpp
 * @brief Implementation of the ModuleCell::ReciprocalGrid base class.
 * @note Spin-free logic migrated from K_Vectors (klist.cpp) and
 *       KVectorUtils (k_vector_utils.cpp) on 2026-08-14.
 */
#include "reciprocal_grid.h"

#include "source_cell/unitcell.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_base/global_function.h"
#include "source_base/formatter.h"
#include "source_base/global_variable.h"
#include "source_base/matrix3.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"

namespace ModuleCell
{

void ReciprocalGrid::renew(const int& kpoint_number)
{
    kvec_c.resize(kpoint_number);
    kvec_d.resize(kpoint_number);
    kvec_c_full.resize(kpoint_number);
    wk.resize(kpoint_number);
    ngk.resize(kpoint_number);

    return;
}

double ReciprocalGrid::Monkhorst_Pack_formula(const int& k_type, const double& offset, const int& n, const int& dim)
{
    double coordinate = 0.0;
    if (k_type == 1)
    {
        coordinate = (offset + 2.0 * (double)n - (double)dim - 1.0) / (2.0 * (double)dim);
    }
    else
    {
        coordinate = (offset + (double)n - 1.0) / (double)dim;
    }
    return coordinate;
}

void ReciprocalGrid::Monkhorst_Pack(const int* nmp_in, const double* koffset_in, const int k_type)
{
    const int mpnx = nmp_in[0];
    const int mpny = nmp_in[1];
    const int mpnz = nmp_in[2];

    this->nkstot = mpnx * mpny * mpnz;
    // only can renew after nkstot is estimated.
    this->renew(nkstot * spin_factor());

    for (int x = 1; x <= mpnx; x++)
    {
        double v1 = Monkhorst_Pack_formula(k_type, koffset_in[0], x, mpnx);
        if (std::abs(v1) < 1.0e-10) {
            v1 = 0.0; // mohan update 2012-06-10
        }
        for (int y = 1; y <= mpny; y++)
        {
            double v2 = Monkhorst_Pack_formula(k_type, koffset_in[1], y, mpny);
            if (std::abs(v2) < 1.0e-10) {
                v2 = 0.0;
            }
            for (int z = 1; z <= mpnz; z++)
            {
                double v3 = Monkhorst_Pack_formula(k_type, koffset_in[2], z, mpnz);
                if (std::abs(v3) < 1.0e-10) {
                    v3 = 0.0;
                }
                // index of nks kpoint
                const int i = mpnx * mpny * (z - 1) + mpnx * (y - 1) + (x - 1);
                kvec_d[i].set(v1, v2, v3);
            }
        }
    }

    const double weight = 1.0 / static_cast<double>(nkstot);
    for (int ik = 0; ik < nkstot; ik++)
    {
        wk[ik] = weight;
    }
    this->kd_done = true;

    return;
}

void ReciprocalGrid::kvec_d2c(const ModuleBase::Matrix3& reciprocal_vec)
{
    if (this->kvec_d.size() != this->kvec_c.size())
    {
        this->kvec_c.resize(this->kvec_d.size());
    }
    int nks = this->kvec_d.size(); // always convert all k vectors

    for (int i = 0; i < nks; i++)
    {
        // mohan fixed bug 2010-1-10
        if (std::abs(this->kvec_d[i].x) < 1.0e-10)
        {
            this->kvec_d[i].x = 0.0;
        }
        if (std::abs(this->kvec_d[i].y) < 1.0e-10)
        {
            this->kvec_d[i].y = 0.0;
        }
        if (std::abs(this->kvec_d[i].z) < 1.0e-10)
        {
            this->kvec_d[i].z = 0.0;
        }

        this->kvec_c[i] = this->kvec_d[i] * reciprocal_vec;

        // mohan add2012-06-10
        if (std::abs(this->kvec_c[i].x) < 1.0e-10)
        {
            this->kvec_c[i].x = 0.0;
        }
        if (std::abs(this->kvec_c[i].y) < 1.0e-10)
        {
            this->kvec_c[i].y = 0.0;
        }
        if (std::abs(this->kvec_c[i].z) < 1.0e-10)
        {
            this->kvec_c[i].z = 0.0;
        }
    }
}

void ReciprocalGrid::kvec_c2d(const ModuleBase::Matrix3& latvec)
{
    if (this->kvec_d.size() != this->kvec_c.size())
    {
        this->kvec_d.resize(this->kvec_c.size());
    }
    int nks = this->kvec_d.size(); // always convert all k vectors

    ModuleBase::Matrix3 RT = latvec.Transpose();
    for (int i = 0; i < nks; i++)
    {
        // mohan fixed bug 2011-03-07
        this->kvec_d[i] = this->kvec_c[i] * RT;
    }
}

void ReciprocalGrid::set_both_kvec(const ModuleBase::Matrix3& G, const ModuleBase::Matrix3& R, std::string& skpt)
{
    if (true) // Originally GlobalV::FINAL_SCF
    {
        if (this->k_nkstot == 0)
        {
            this->kd_done = true;
            this->kc_done = false;
        }
        else
        {
            if (this->k_kword == "Cartesian" || this->k_kword == "C")
            {
                this->kc_done = true;
                this->kd_done = false;
            }
            else if (this->k_kword == "Direct" || this->k_kword == "D")
            {
                this->kd_done = true;
                this->kc_done = false;
            }
            else
            {
                GlobalV::ofs_warning << " Error : neither Cartesian nor Direct kpoint." << std::endl;
            }
        }
    }

    // set cartesian k vectors.
    if (!this->kc_done && this->kd_done)
    {
        this->kvec_d2c(G);
        this->kc_done = true;
    }

    // set direct k vectors
    else if (this->kc_done && !this->kd_done)
    {
        this->kvec_c2d(R);
        this->kd_done = true;
    }
    std::string table;
    table += " K-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < this->nkstot; i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 this->kvec_d[i].x,
                                 this->kvec_d[i].y,
                                 this->kvec_d[i].z,
                                 this->wk[i]);
    }
    GlobalV::ofs_running << table << std::endl;
    if (GlobalV::MY_RANK == 0)
    {
        std::stringstream ss;
        ss << " " << std::setw(40) << "nkstot now"
           << " = " << this->nkstot << std::endl;
        ss << table << std::endl;
        skpt = ss.str();
    }
    return;
}

void ReciprocalGrid::normalize_wk(const int& degspin)
{
    if (GlobalV::MY_RANK != 0) {
        return;
    }
    double sum = 0.0;

    for (int ik = 0; ik < nkstot; ik++)
    {
        sum += this->wk[ik];
    }

    // If sum of weights is zero or very small, set equal weights
    if (sum < 1e-10)
    {
        ModuleBase::WARNING("ReciprocalGrid::normalize_wk",
                            "Sum of k-point weights is zero or very small. "
                            "Setting equal weights for all k-points.");
        for (int ik = 0; ik < nkstot; ik++)
        {
            this->wk[ik] = 1.0 / double(nkstot);
        }
        sum = 1.0;
    }

    for (int ik = 0; ik < nkstot; ik++)
    {
        this->wk[ik] /= sum;
    }

    for (int ik = 0; ik < nkstot; ik++)
    {
        this->wk[ik] *= degspin;
    }

    return;
}

void ReciprocalGrid::print_klists(std::ofstream& ofs) const
{
    ModuleBase::TITLE("ReciprocalGrid", "print_klists");
    int nks = this->nks;
    int nkstot = this->nkstot;

    if (nkstot < nks)
    {
        std::cout << "\n nkstot=" << nkstot;
        std::cout << "\n nks=" << nks;
        ModuleBase::WARNING_QUIT("print_klists", "nkstot < nks");
    }
    std::string table;
    table += " K-POINTS CARTESIAN COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "CARTESIAN_X", "CARTESIAN_Y", "CARTESIAN_Z", "WEIGHT");
    for (int i = 0; i < nks; i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 this->kvec_c[i].x,
                                 this->kvec_c[i].y,
                                 this->kvec_c[i].z,
                                 this->wk[i]);
    }
    GlobalV::ofs_running << "\n" << table << std::endl;

    table.clear();
    table += " K-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < nks; i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 this->kvec_d[i].x,
                                 this->kvec_d[i].y,
                                 this->kvec_d[i].z,
                                 this->wk[i]);
    }
    GlobalV::ofs_running << "\n" << table << std::endl;
    return;
}

void ReciprocalGrid::reduce_ibz(const ModuleBase::Matrix3* rot_ops,
                                int nrotkm,
                                const ModuleBase::Matrix3& G,
                                const ModuleBase::Matrix3& k_lattice,
                                const ModuleBase::Matrix3* kkmatrix,
                                double epsilon,
                                std::vector<ModuleBase::Vector3<double>>& vec_ibz,
                                std::vector<double>& wk_ibz,
                                std::vector<int>& ibz_index,
                                std::vector<int>& ibz2bz)
{
    auto equal = [epsilon](double m, double n) { return fabs(m - n) < epsilon; };
    // restrict a vector to (-0.5, 0.5]
    auto restrict_kpt = [epsilon](ModuleBase::Vector3<double>& kvec) {
        kvec.x = fmod(kvec.x + 100.5 - 0.5 * epsilon, 1) - 0.5 + 0.5 * epsilon;
        kvec.y = fmod(kvec.y + 100.5 - 0.5 * epsilon, 1) - 0.5 + 0.5 * epsilon;
        kvec.z = fmod(kvec.z + 100.5 - 0.5 * epsilon, 1) - 0.5 + 0.5 * epsilon;
        if (std::abs(kvec.x) < epsilon)
        {
            kvec.x = 0.0;
        }
        if (std::abs(kvec.y) < epsilon)
        {
            kvec.y = 0.0;
        }
        if (std::abs(kvec.z) < epsilon)
        {
            kvec.z = 0.0;
        }
        return;
    };

    // direct coordinates of points in the k-lattice
    std::vector<ModuleBase::Vector3<double>> kvec_d_k(this->nkstot);
    if (this->is_mp)
    {
        for (int i = 0; i < this->nkstot; ++i)
        {
            kvec_d_k[i] = this->kvec_d[i] * G * k_lattice.Inverse();
        }
    }

    int nkstot_ibz = 0;

    assert(this->nkstot > 0);
    std::vector<ModuleBase::Vector3<double>> kvec_d_ibz(this->nkstot);
    std::vector<double> wk_ibz_tmp(this->nkstot); // ibz point weight
    ibz2bz.resize(this->nkstot);

    // nkstot is the total input points number.
    double weight = 1.0 / static_cast<double>(this->nkstot);

    ModuleBase::Vector3<double> kvec_rot;
    ModuleBase::Vector3<double> kvec_rot_k;

    // update map k -> irreducible k
    ibz_index.assign(this->nkstot_full, -1); // -1 means not in ibz list
    // search in all k-points.
    for (int i = 0; i < this->nkstot; ++i)
    {
        if (!this->is_mp) { weight = this->wk[i]; } // use the input weight, instead of 1/nkstot

        // restrict to (-0.5, 0.5]
        restrict_kpt(this->kvec_d[i]);

        bool already_exist = false;
        int exist_number = -1;
        // search over all symmetry operations
        for (int j = 0; j < nrotkm; ++j)
        {
            if (!already_exist)
            {
                kvec_rot = this->kvec_d[i] * rot_ops[j]; // wrong for total energy, but correct for nonlocal force.
                restrict_kpt(kvec_rot);
                if (this->is_mp)
                {
                    kvec_rot_k = kvec_d_k[i] * kkmatrix[j];              // k-lattice rotation
                    kvec_rot_k = kvec_rot_k * k_lattice * G.Inverse();   // convert to recip lattice
                    restrict_kpt(kvec_rot_k);

                    assert(equal(kvec_rot.x, kvec_rot_k.x));
                    assert(equal(kvec_rot.y, kvec_rot_k.y));
                    assert(equal(kvec_rot.z, kvec_rot_k.z));
                    kvec_rot_k = kvec_rot_k * G * k_lattice.Inverse(); // convert back to k-lattice
                }
                for (int k = 0; k < nkstot_ibz; ++k)
                {
                    if (equal(kvec_rot.x, kvec_d_ibz[k].x) && equal(kvec_rot.y, kvec_d_ibz[k].y)
                        && equal(kvec_rot.z, kvec_d_ibz[k].z))
                    {
                        already_exist = true;
                        // find another ibz point,
                        // but is already in the ibz list.
                        // so the weight need to +1;
                        wk_ibz_tmp[k] += weight;
                        exist_number = k;
                        break;
                    }
                }
            } // end !already_exist
        }
        // if really there is no equivalent point in the list, then add it.
        if (!already_exist)
        {
            kvec_d_ibz[nkstot_ibz] = this->kvec_d[i];
            ibz_index[i] = nkstot_ibz;

            // the weight should be averaged point weight.
            wk_ibz_tmp[nkstot_ibz] = weight;

            // ibz2bz records the index of origin points.
            ibz2bz[nkstot_ibz] = i;
            ++nkstot_ibz;
        }
        else
        {
            double kmol_new = this->kvec_d[i].norm2();
            double kmol_old = kvec_d_ibz[exist_number].norm2();

            ibz_index[i] = exist_number;

            // why we need this step?
            // because in pw_basis.cpp, while calculate ggwfc2,
            // if we want to keep the result of symmetry operation is right.
            // we need to fix the number of plane wave.
            // and the number of plane wave is depending on the |K+G|,
            // so we need to |K|max to be the same as 'no symmetry'.
            // mohan 2010-01-30
            if (kmol_new > kmol_old)
            {
                kvec_d_ibz[exist_number] = this->kvec_d[i];
            }
        }
    }

    vec_ibz.resize(nkstot_ibz);
    wk_ibz.resize(nkstot_ibz);
    ibz2bz.resize(nkstot_ibz);
    for (int i = 0; i < nkstot_ibz; ++i)
    {
        vec_ibz[i] = kvec_d_ibz[i];
        wk_ibz[i] = wk_ibz_tmp[i];
    }

    return;
}

bool ReciprocalGrid::build_star_ops(const UnitCell& ucell,
                                    const ModuleSymmetry::Symmetry& symm,
                                    bool use_symm,
                                    ModuleBase::Matrix3& k_vec,
                                    std::vector<ModuleBase::Matrix3>& kgmatrix,
                                    int& nrotkm) const
{
    // k-lattice: "pricell" of reciprocal space
    // CAUTION: should fit into all k-input method, not only MP  !!!
    // the basis vector of reciprocal lattice: recip_vec1, recip_vec2, recip_vec3
    ModuleBase::Vector3<double> recip_vec1(ucell.G.e11, ucell.G.e12, ucell.G.e13);
    ModuleBase::Vector3<double> recip_vec2(ucell.G.e21, ucell.G.e22, ucell.G.e23);
    ModuleBase::Vector3<double> recip_vec3(ucell.G.e31, ucell.G.e32, ucell.G.e33);
    ModuleBase::Vector3<double> k_vec1, k_vec2, k_vec3;
    if (this->is_mp)
    {
        k_vec1 = ModuleBase::Vector3<double>(recip_vec1.x / this->nmp[0], recip_vec1.y / this->nmp[0], recip_vec1.z / this->nmp[0]);
        k_vec2 = ModuleBase::Vector3<double>(recip_vec2.x / this->nmp[1], recip_vec2.y / this->nmp[1], recip_vec2.z / this->nmp[1]);
        k_vec3 = ModuleBase::Vector3<double>(recip_vec3.x / this->nmp[2], recip_vec3.y / this->nmp[2], recip_vec3.z / this->nmp[2]);
        k_vec = ModuleBase::Matrix3(k_vec1.x,
                                    k_vec1.y,
                                    k_vec1.z,
                                    k_vec2.x,
                                    k_vec2.y,
                                    k_vec2.z,
                                    k_vec3.x,
                                    k_vec3.y,
                                    k_vec3.z);
    }

    ModuleBase::Matrix3 inv(-1, 0, 0, 0, -1, 0, 0, 0, -1);
    ModuleBase::Matrix3 ind(1, 0, 0, 0, 1, 0, 0, 0, 1);

    nrotkm = 0;
    if (use_symm)
    {
        // bravais type of reciprocal lattice and k-lattice

        double recip_vec_const[6];
        double recip_vec0_const[6];
        double k_vec_const[6];
        double k_vec0_const[6];
        int recip_brav_type = 15;
        int k_brav_type = 15;
        std::string recip_brav_name;
        std::string k_brav_name;
        ModuleBase::Vector3<double> k_vec01 = k_vec1, k_vec02 = k_vec2, k_vec03 = k_vec3;

        // determine the Bravais type and related parameters of the lattice
        symm.lattice_type(recip_vec1,
                          recip_vec2,
                          recip_vec3,
                          recip_vec1,
                          recip_vec2,
                          recip_vec3,
                          recip_vec_const,
                          recip_vec0_const,
                          recip_brav_type,
                          recip_brav_name,
                          ucell.atoms,
                          false,
                          nullptr,
                          1e-6);
        GlobalV::ofs_running << "\n For reciprocal-space lattice" << std::endl;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Bravais lattice type", recip_brav_type);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Bravais lattice name", recip_brav_name);

        // the map of bravis lattice from real to reciprocal space
        // for example, 3(fcc) in real space matches 2(bcc) in reciprocal space
        std::vector<int> ibrav_a2b{1, 3, 2, 4, 5, 6, 7, 8, 10, 9, 11, 12, 13, 14};
        // check if the reciprocal lattice is compatible with the real space lattice
        auto ibrav_match = [&](int ibrav_b) -> bool {
            const int& ibrav_a = symm.real_brav;
            if (ibrav_a < 1 || ibrav_a > 14)
            {
                return false;
            }
            return (ibrav_b == ibrav_a2b[ibrav_a - 1]);
        };
        if (!ibrav_match(recip_brav_type)) // if not match, exit and return
        {
            GlobalV::ofs_running << "Error: Bravais lattice type of reciprocal lattice is not compatible with that of "
                                    "real space lattice:"
                                 << std::endl;
            GlobalV::ofs_running << "ibrav of real space lattice: " << symm.ilattname << std::endl;
            GlobalV::ofs_running << "ibrav of reciprocal lattice: " << recip_brav_name << std::endl;
            GlobalV::ofs_running << "(which should be " << ibrav_a2b[symm.real_brav - 1] << ")." << std::endl;
            return false;
        }

        // if match, continue
        if (this->is_mp)
        {
            symm.lattice_type(k_vec1,
                              k_vec2,
                              k_vec3,
                              k_vec01,
                              k_vec02,
                              k_vec03,
                              k_vec_const,
                              k_vec0_const,
                              k_brav_type,
                              k_brav_name,
                              ucell.atoms,
                              false,
                              nullptr,
                              1e-6);
            GlobalV::ofs_running << "\n For k-vectors" << std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Bravais lattice type", k_brav_type);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Bravais lattice name", k_brav_name);
        }
        // point-group analysis of reciprocal lattice
        ModuleBase::Matrix3 bsymop[48];
        int bnop = 0;
        // search again
        symm.lattice_type(recip_vec1,
                          recip_vec2,
                          recip_vec3,
                          recip_vec1,
                          recip_vec2,
                          recip_vec3,
                          recip_vec_const,
                          recip_vec0_const,
                          recip_brav_type,
                          recip_brav_name,
                          ucell.atoms,
                          false,
                          nullptr,
                          1e-6);
        ModuleBase::Matrix3 b_optlat_new(recip_vec1.x, recip_vec1.y, recip_vec1.z,
                                         recip_vec2.x, recip_vec2.y, recip_vec2.z,
                                         recip_vec3.x, recip_vec3.y, recip_vec3.z);
        // set the crystal point-group symmetry operation
        const int cal_symm_repr[2] = {0, 6};
        symm.setgroup(bsymop, bnop, recip_brav_type, cal_symm_repr);
        // transform the above symmetric operation matrices between different coordinate
        symm.gmatrix_convert(bsymop, bsymop, bnop, b_optlat_new, ucell.G);

        // check if all the kgmatrix are in bsymop
        auto matequal = [&symm](ModuleBase::Matrix3 a, ModuleBase::Matrix3 b) {
            return (symm.equal(a.e11, b.e11) && symm.equal(a.e12, b.e12) && symm.equal(a.e13, b.e13)
                    && symm.equal(a.e21, b.e21) && symm.equal(a.e22, b.e22) && symm.equal(a.e23, b.e23)
                    && symm.equal(a.e31, b.e31) && symm.equal(a.e32, b.e32) && symm.equal(a.e33, b.e33));
        };
        for (int i = 0; i < symm.nrotk; ++i)
        {
            bool found = false;
            for (int j = 0; j < bnop; ++j)
            {
                if (matequal(symm.kgmatrix[i], bsymop[j]))
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                return false;
            }
        }
        nrotkm = symm.nrotk;
        for (int i = 0; i < nrotkm; ++i)
        {
            kgmatrix[i] = symm.kgmatrix[i];
        }
    }
    else if (this->is_mp) // only include for Monkhorst-Pack grid
    {
        nrotkm = 2;
        kgmatrix[0] = ind;
        kgmatrix[1] = inv;
    }

    return true;
}

} // namespace ModuleCell
