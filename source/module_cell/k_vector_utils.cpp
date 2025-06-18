//
// Created by rhx on 25-6-3.
//
#include "k_vector_utils.h"

#include "klist.h"
#include "source_base/global_variable.h"
#include "source_base/matrix3.h"

#include "source_base/formatter.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"
#include "module_parameter/parameter.h"

namespace KVectorUtils
{
void kvec_d2c(K_Vectors& kv, const ModuleBase::Matrix3& reciprocal_vec)
{
    //    throw std::runtime_error("k_vec_d2c: This function is not implemented in the new codebase. Please use the new
    //    implementation.");
    if (kv.kvec_d.size() != kv.kvec_c.size())
    {
        //        ModuleBase::WARNING_QUIT("k_vec_d2c", "Size of Cartesian and Direct K vectors mismatch. ");
        kv.kvec_c.resize(kv.kvec_d.size());
    }
    int nks = kv.kvec_d.size(); // always convert all k vectors

    for (int i = 0; i < nks; i++)
    {
        // wrong!!   kvec_c[i] = G * kvec_d[i];
        //  mohan fixed bug 2010-1-10
        if (std::abs(kv.kvec_d[i].x) < 1.0e-10)
        {
            kv.kvec_d[i].x = 0.0;
        }
        if (std::abs(kv.kvec_d[i].y) < 1.0e-10)
        {
            kv.kvec_d[i].y = 0.0;
        }
        if (std::abs(kv.kvec_d[i].z) < 1.0e-10)
        {
            kv.kvec_d[i].z = 0.0;
        }

        kv.kvec_c[i] = kv.kvec_d[i] * reciprocal_vec;

        // mohan add2012-06-10
        if (std::abs(kv.kvec_c[i].x) < 1.0e-10)
        {
            kv.kvec_c[i].x = 0.0;
        }
        if (std::abs(kv.kvec_c[i].y) < 1.0e-10)
        {
            kv.kvec_c[i].y = 0.0;
        }
        if (std::abs(kv.kvec_c[i].z) < 1.0e-10)
        {
            kv.kvec_c[i].z = 0.0;
        }
    }
}
void kvec_c2d(K_Vectors& kv, const ModuleBase::Matrix3& latvec)
{
    if (kv.kvec_d.size() != kv.kvec_c.size())
    {
        kv.kvec_d.resize(kv.kvec_c.size());
    }
    int nks = kv.kvec_d.size(); // always convert all k vectors

    ModuleBase::Matrix3 RT = latvec.Transpose();
    for (int i = 0; i < nks; i++)
    {
        //			std::cout << " ik=" << i
        //				<< " kvec.x=" << kvec_c[i].x
        //				<< " kvec.y=" << kvec_c[i].y
        //				<< " kvec.z=" << kvec_c[i].z << std::endl;
        // wrong!            kvec_d[i] = RT * kvec_c[i];
        // mohan fixed bug 2011-03-07
        kv.kvec_d[i] = kv.kvec_c[i] * RT;
    }
}

void set_both_kvec(K_Vectors& kv, const ModuleBase::Matrix3& G, const ModuleBase::Matrix3& R, std::string& skpt)
{
    if (true) // Originally GlobalV::FINAL_SCF, but we don't have this variable in the new code.
    {
        if (kv.get_k_nkstot() == 0)
        {
            kv.kd_done = true;
            kv.kc_done = false;
        }
        else
        {
            if (kv.get_k_kword() == "Cartesian" || kv.get_k_kword() == "C")
            {
                kv.kc_done = true;
                kv.kd_done = false;
            }
            else if (kv.get_k_kword() == "Direct" || kv.get_k_kword() == "D")
            {
                kv.kd_done = true;
                kv.kc_done = false;
            }
            else
            {
                GlobalV::ofs_warning << " Error : neither Cartesian nor Direct kpoint." << std::endl;
            }
        }
    }

    // set cartesian k vectors.
    if (!kv.kc_done && kv.kd_done)
    {
        KVectorUtils::kvec_d2c(kv, G);
        kv.kc_done = true;
    }

    // set direct k vectors
    else if (kv.kc_done && !kv.kd_done)
    {
        KVectorUtils::kvec_c2d(kv, R);
        kv.kd_done = true;
    }
    std::string table;
    table += " K-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < kv.get_nkstot(); i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 kv.kvec_d[i].x,
                                 kv.kvec_d[i].y,
                                 kv.kvec_d[i].z,
                                 kv.wk[i]);
    }
    GlobalV::ofs_running << table << std::endl;
    if (GlobalV::MY_RANK == 0)
    {
        std::stringstream ss;
        ss << " " << std::setw(40) << "nkstot now"
           << " = " << kv.get_nkstot() << std::endl;
        ss << table << std::endl;
        skpt = ss.str();
    }
    return;
}

void set_after_vc(K_Vectors& kv, const int& nspin_in, const ModuleBase::Matrix3& reciprocal_vec)
{
    GlobalV::ofs_running << "\n SETUP K-POINTS" << std::endl;
    //    kv.nspin = nspin_in;
    kv.set_nspin(nspin_in);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nspin", kv.get_nspin());

    // set cartesian k vectors.
    KVectorUtils::kvec_d2c(kv, reciprocal_vec);

    std::string table;
    table += "K-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < kv.get_nks(); i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 kv.kvec_d[i].x,
                                 kv.kvec_d[i].y,
                                 kv.kvec_d[i].z,
                                 kv.wk[i]);
    }
    GlobalV::ofs_running << table << std::endl;

    kv.kd_done = true;
    kv.kc_done = true;

    print_klists(kv, GlobalV::ofs_running);
}

void print_klists(const K_Vectors& kv, std::ofstream& ofs)
{
    ModuleBase::TITLE("KVectorUtils", "print_klists");
    int nks = kv.get_nks();
    int nkstot = kv.get_nkstot();

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
                                 kv.kvec_c[i].x,
                                 kv.kvec_c[i].y,
                                 kv.kvec_c[i].z,
                                 kv.wk[i]);
    }
    GlobalV::ofs_running << "\n" << table << std::endl;

    table.clear();
    table += " K-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < nks; i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 kv.kvec_d[i].x,
                                 kv.kvec_d[i].y,
                                 kv.kvec_d[i].z,
                                 kv.wk[i]);
    }
    GlobalV::ofs_running << "\n" << table << std::endl;
    return;
}

#ifdef __MPI
void kvec_mpi_k(K_Vectors& kv)
{
    ModuleBase::TITLE("KVectorUtils", "kvec_mpi_k");

    Parallel_Common::bcast_bool(kv.kc_done);

    Parallel_Common::bcast_bool(kv.kd_done);

    Parallel_Common::bcast_int(kv.nspin);

    Parallel_Common::bcast_int(kv.nkstot);

    Parallel_Common::bcast_int(kv.nkstot_full);

    Parallel_Common::bcast_int(kv.nmp, 3);

    kv.kl_segids.resize(kv.nkstot);
    Parallel_Common::bcast_int(kv.kl_segids.data(), kv.nkstot);

    Parallel_Common::bcast_double(kv.koffset, 3);

    kv.nks = kv.para_k.nks_pool[GlobalV::MY_POOL];

    GlobalV::ofs_running << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "k-point number in this process", kv.nks);
    int nks_minimum = kv.nks;

    Parallel_Reduce::gather_min_int_all(GlobalV::NPROC, nks_minimum);

    if (nks_minimum == 0)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::mpi_k()", " nks == 0, some processor have no k point!");
    }
    else
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "minimum distributed K point number", nks_minimum);
    }

    std::vector<int> isk_aux(kv.nkstot);
    std::vector<double> wk_aux(kv.nkstot);
    std::vector<double> kvec_c_aux(kv.nkstot * 3);
    std::vector<double> kvec_d_aux(kv.nkstot * 3);

    // collect and process in rank 0
    if (GlobalV::MY_RANK == 0)
    {
        for (int ik = 0; ik < kv.nkstot; ik++)
        {
            isk_aux[ik] = kv.isk[ik];
            wk_aux[ik] = kv.wk[ik];
            kvec_c_aux[3 * ik] = kv.kvec_c[ik].x;
            kvec_c_aux[3 * ik + 1] = kv.kvec_c[ik].y;
            kvec_c_aux[3 * ik + 2] = kv.kvec_c[ik].z;
            kvec_d_aux[3 * ik] = kv.kvec_d[ik].x;
            kvec_d_aux[3 * ik + 1] = kv.kvec_d[ik].y;
            kvec_d_aux[3 * ik + 2] = kv.kvec_d[ik].z;
        }
    }

    // broadcast k point data to all processors
    Parallel_Common::bcast_int(isk_aux.data(), kv.nkstot);

    Parallel_Common::bcast_double(wk_aux.data(), kv.nkstot);
    Parallel_Common::bcast_double(kvec_c_aux.data(), kv.nkstot * 3);
    Parallel_Common::bcast_double(kvec_d_aux.data(), kv.nkstot * 3);

    // process k point data in each processor
    kv.renew(kv.nks * kv.nspin);

    // distribute
    int k_index = 0;

    for (int i = 0; i < kv.nks; i++)
    {
        // 3 is because each k point has three value:kx, ky, kz
        k_index = i + kv.para_k.startk_pool[GlobalV::MY_POOL];
        kv.kvec_c[i].x = kvec_c_aux[k_index * 3];
        kv.kvec_c[i].y = kvec_c_aux[k_index * 3 + 1];
        kv.kvec_c[i].z = kvec_c_aux[k_index * 3 + 2];
        kv.kvec_d[i].x = kvec_d_aux[k_index * 3];
        kv.kvec_d[i].y = kvec_d_aux[k_index * 3 + 1];
        kv.kvec_d[i].z = kvec_d_aux[k_index * 3 + 2];
        kv.wk[i] = wk_aux[k_index];
        kv.isk[i] = isk_aux[k_index];
    }

#ifdef __EXX
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    { // bcast kstars
        kv.kstars.resize(kv.nkstot);
        for (int ikibz = 0; ikibz < kv.nkstot; ++ikibz)
        {
            int starsize = kv.kstars[ikibz].size();
            Parallel_Common::bcast_int(starsize);
            GlobalV::ofs_running << "starsize: " << starsize << std::endl;
            auto ks = kv.kstars[ikibz].begin();
            for (int ik = 0; ik < starsize; ++ik)
            {
                int isym = 0;
                ModuleBase::Vector3<double> ks_vec(0, 0, 0);
                if (GlobalV::MY_RANK == 0)
                {
                    isym = ks->first;
                    ks_vec = ks->second;
                    ++ks;
                }
                Parallel_Common::bcast_int(isym);
                Parallel_Common::bcast_double(ks_vec.x);
                Parallel_Common::bcast_double(ks_vec.y);
                Parallel_Common::bcast_double(ks_vec.z);
                GlobalV::ofs_running << "isym: " << isym << " ks_vec: " << ks_vec.x << " " << ks_vec.y << " "
                                     << ks_vec.z << std::endl;
                if (GlobalV::MY_RANK != 0)
                {
                    kv.kstars[ikibz].insert(std::make_pair(isym, ks_vec));
                }
            }
        }
    }
#endif
} // END SUBROUTINE
#endif


void kvec_ibz_kpoint(K_Vectors& kv,
                     const ModuleSymmetry::Symmetry& symm,
                     bool use_symm,
                     std::string& skpt,
                     const UnitCell& ucell,
                     bool& match)
{
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    ModuleBase::TITLE("K_Vectors", "ibz_kpoint");

    // k-lattice: "pricell" of reciprocal space
    // CAUTION: should fit into all k-input method, not only MP  !!!
    // the basis vector of reciprocal lattice: recip_vec1, recip_vec2, recip_vec3
    ModuleBase::Vector3<double> recip_vec1(ucell.G.e11, ucell.G.e12, ucell.G.e13);
    ModuleBase::Vector3<double> recip_vec2(ucell.G.e21, ucell.G.e22, ucell.G.e23);
    ModuleBase::Vector3<double> recip_vec3(ucell.G.e31, ucell.G.e32, ucell.G.e33);
    ModuleBase::Vector3<double> k_vec1, k_vec2, k_vec3;
    ModuleBase::Matrix3 k_vec;
    if (kv.get_is_mp())
    {
        k_vec1 = ModuleBase::Vector3<double>(recip_vec1.x / kv.nmp[0], recip_vec1.y / kv.nmp[0], recip_vec1.z / kv.nmp[0]);
        k_vec2 = ModuleBase::Vector3<double>(recip_vec2.x / kv.nmp[1], recip_vec2.y / kv.nmp[1], recip_vec2.z / kv.nmp[1]);
        k_vec3 = ModuleBase::Vector3<double>(recip_vec3.x / kv.nmp[2], recip_vec3.y / kv.nmp[2], recip_vec3.z / kv.nmp[2]);
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

    //===============================================
    // search in all space group operations
    // if the operations does not already included
    // inverse operation, double it.
    //===============================================
    bool include_inv = false;
    std::vector<ModuleBase::Matrix3> kgmatrix(48 * 2);
    ModuleBase::Matrix3 inv(-1, 0, 0, 0, -1, 0, 0, 0, -1);
    ModuleBase::Matrix3 ind(1, 0, 0, 0, 1, 0, 0, 0, 1);

    int nrotkm = 0;
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

        // it's not necessary to calculate gb01, gb02, gb03,
        // because they are only used as a vector, no need to be assigned values

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
                          nullptr);
        GlobalV::ofs_running << "\n For reciprocal-space lattice:" << std::endl;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "BRAVAIS TYPE", recip_brav_type);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "BRAVAIS LATTICE NAME", recip_brav_name);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ibrav", recip_brav_type);

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
            match = false;
            return;
        }

        // if match, continue
        if (kv.get_is_mp())
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
                              nullptr);
            GlobalV::ofs_running << "\n For k vectors:" << std::endl;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "BRAVAIS TYPE", k_brav_type);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "BRAVAIS LATTICE NAME", k_brav_name);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ibrav", k_brav_type);
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
                          nullptr);
        ModuleBase::Matrix3 b_optlat_new(recip_vec1.x,
                                         recip_vec1.y,
                                         recip_vec1.z,
                                         recip_vec2.x,
                                         recip_vec2.y,
                                         recip_vec2.z,
                                         recip_vec3.x,
                                         recip_vec3.y,
                                         recip_vec3.z);
        // set the crystal point-group symmetry operation
        symm.setgroup(bsymop, bnop, recip_brav_type);
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
            match = false;
            for (int j = 0; j < bnop; ++j)
            {
                if (matequal(symm.kgmatrix[i], bsymop[j]))
                {
                    match = true;
                    break;
                }
            }
            if (!match)
            {
                return;
            }
        }
        nrotkm = symm.nrotk; // change if inv not included
        for (int i = 0; i < nrotkm; ++i)
        {
            if (symm.kgmatrix[i] == inv)
            {
                include_inv = true;
            }
            kgmatrix[i] = symm.kgmatrix[i];
        }

        if (!include_inv)
        {
            for (int i = 0; i < symm.nrotk; ++i)
            {
                kgmatrix[i + symm.nrotk] = inv * symm.kgmatrix[i];
            }
            nrotkm = 2 * symm.nrotk;
        }
    }
    else if (kv.get_is_mp()) // only include for mp grid
    {
        nrotkm = 2;
        kgmatrix[0] = ind;
        kgmatrix[1] = inv;
    }
    else
    {
        return;
    }

    // convert kgmatrix to k-lattice
    ModuleBase::Matrix3* kkmatrix = new ModuleBase::Matrix3[nrotkm];
    if (kv.get_is_mp())
    {
        symm.gmatrix_convert(kgmatrix.data(), kkmatrix, nrotkm, ucell.G, k_vec);
    }
    // direct coordinates of k-points in k-lattice
    std::vector<ModuleBase::Vector3<double>> kvec_d_k(kv.get_nkstot());
    if (kv.get_is_mp())
    {
        for (int i = 0; i < kv.get_nkstot(); ++i)
        {
            kvec_d_k[i] = kv.kvec_d[i] * ucell.G * k_vec.Inverse();
        }
    }

    // use operation : kgmatrix to find
    // the new set kvec_d : ir_kpt
    int nkstot_ibz = 0;

    assert(kv.get_nkstot() > 0);
    std::vector<ModuleBase::Vector3<double>> kvec_d_ibz(kv.get_nkstot());
    std::vector<double> wk_ibz(kv.get_nkstot()); // ibz kpoint wk ,weight of k points
    std::vector<int> ibz2bz(kv.get_nkstot());

    // nkstot is the total input k-points number.
    const double weight = 1.0 / static_cast<double>(kv.get_nkstot());

    ModuleBase::Vector3<double> kvec_rot;
    ModuleBase::Vector3<double> kvec_rot_k;

    //	for(int i=0; i<nrotkm; i++)
    //	{
    //		out.printM3("rot matrix",kgmatrix[i]);
    //	}
    auto restrict_kpt = [&symm](ModuleBase::Vector3<double>& kvec) {
        // in (-0.5, 0.5]
        kvec.x = fmod(kvec.x + 100.5 - 0.5 * symm.epsilon, 1) - 0.5 + 0.5 * symm.epsilon;
        kvec.y = fmod(kvec.y + 100.5 - 0.5 * symm.epsilon, 1) - 0.5 + 0.5 * symm.epsilon;
        kvec.z = fmod(kvec.z + 100.5 - 0.5 * symm.epsilon, 1) - 0.5 + 0.5 * symm.epsilon;
        // in [0, 1)
        // kvec.x = fmod(kvec.x + 100 + symm.epsilon, 1) - symm.epsilon;
        // kvec.y = fmod(kvec.y + 100 + symm.epsilon, 1) - symm.epsilon;
        // kvec.z = fmod(kvec.z + 100 + symm.epsilon, 1) - symm.epsilon;
        if (std::abs(kvec.x) < symm.epsilon)
        {
            kvec.x = 0.0;
        }
        if (std::abs(kvec.y) < symm.epsilon)
        {
            kvec.y = 0.0;
        }
        if (std::abs(kvec.z) < symm.epsilon)
        {
            kvec.z = 0.0;
        }
        return;
    };
    // for output in kpoints file
    int ibz_index[kv.get_nkstot()];
    // search in all k-poins.
    for (int i = 0; i < kv.get_nkstot(); ++i)
    {
        // restrict to [0, 1)
        restrict_kpt(kv.kvec_d[i]);

        // std::cout << "\n kpoint = " << i << std::endl;
        // std::cout << "\n kvec_d = " << kvec_d[i].x << " " << kvec_d[i].y << " " << kvec_d[i].z;
        bool already_exist = false;
        int exist_number = -1;
        // search over all symmetry operations
        for (int j = 0; j < nrotkm; ++j)
        {
            if (!already_exist)
            {
                // rotate the kvec_d within all operations.
                // here use direct coordinates.
                //                kvec_rot = kgmatrix[j] * kvec_d[i];
                // mohan modify 2010-01-30.
                // mohan modify again 2010-01-31
                // fix the bug like kvec_d * G; is wrong
                kvec_rot = kv.kvec_d[i] * kgmatrix[j]; // wrong for total energy, but correct for nonlocal force.
                // kvec_rot = kgmatrix[j] * kvec_d[i]; //correct for total energy, but wrong for nonlocal force.
                restrict_kpt(kvec_rot);
                if (kv.get_is_mp())
                {
                    kvec_rot_k = kvec_d_k[i] * kkmatrix[j];              // k-lattice rotation
                    kvec_rot_k = kvec_rot_k * k_vec * ucell.G.Inverse(); // convert to recip lattice
                    restrict_kpt(kvec_rot_k);

                    assert(symm.equal(kvec_rot.x, kvec_rot_k.x));
                    assert(symm.equal(kvec_rot.y, kvec_rot_k.y));
                    assert(symm.equal(kvec_rot.z, kvec_rot_k.z));
                    // std::cout << "\n kvec_rot (in recip) = " << kvec_rot.x << " " << kvec_rot.y << " " << kvec_rot.z;
                    // std::cout << "\n kvec_rot(k to recip)= " << kvec_rot_k.x << " " << kvec_rot_k.y << " " <<
                    // kvec_rot_k.z;
                    kvec_rot_k = kvec_rot_k * ucell.G * k_vec.Inverse(); // convert back to k-latice
                }
                for (int k = 0; k < nkstot_ibz; ++k)
                {
                    if (symm.equal(kvec_rot.x, kvec_d_ibz[k].x) && symm.equal(kvec_rot.y, kvec_d_ibz[k].y)
                        && symm.equal(kvec_rot.z, kvec_d_ibz[k].z))
                    {
                        already_exist = true;
                        // find another ibz k point,
                        // but is already in the ibz_kpoint list.
                        // so the weight need to +1;
                        wk_ibz[k] += weight;
                        exist_number = k;
                        break;
                    }
                }
            } // end !already_exist
        }
        // if really there is no equivalent k point in the list, then add it.
        if (!already_exist)
        {
            // if it's a new ibz kpoint.
            // nkstot_ibz indicate the index of ibz kpoint.
            kvec_d_ibz[nkstot_ibz] = kv.kvec_d[i];
            // output in kpoints file
            ibz_index[i] = nkstot_ibz;

            // the weight should be averged k-point weight.
            wk_ibz[nkstot_ibz] = weight;

            // ibz2bz records the index of origin k points.
            ibz2bz[nkstot_ibz] = i;
            ++nkstot_ibz;
        }
        else // mohan fix bug 2010-1-30
        {
            //			std::cout << "\n\n already exist ! ";

            //			std::cout << "\n kvec_rot = " << kvec_rot.x << " " << kvec_rot.y << " " << kvec_rot.z;
            //			std::cout << "\n kvec_d_ibz = " << kvec_d_ibz[exist_number].x
            //			<< " " << kvec_d_ibz[exist_number].y
            //			<< " " << kvec_d_ibz[exist_number].z;

            double kmol_new = kv.kvec_d[i].norm2();
            double kmol_old = kvec_d_ibz[exist_number].norm2();

            ibz_index[i] = exist_number;

            //			std::cout << "\n kmol_new = " << kmol_new;
            //			std::cout << "\n kmol_old = " << kmol_old;

            // why we need this step?
            // because in pw_basis.cpp, while calculate ggwfc2,
            // if we want to keep the result of symmetry operation is right.
            // we need to fix the number of plane wave.
            // and the number of plane wave is depending on the |K+G|,
            // so we need to |K|max to be the same as 'no symmetry'.
            // mohan 2010-01-30
            if (kmol_new > kmol_old)
            {
                kvec_d_ibz[exist_number] = kv.kvec_d[i];
            }
        }
        //		BLOCK_HERE("check k point");
    }

    delete[] kkmatrix;

#ifdef __EXX
    // setup kstars according to the final (max-norm) kvec_d_ibz
    kv.kstars.resize(nkstot_ibz);
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        for (int i = 0; i < kv.get_nkstot(); ++i)
        {
            int exist_number = -1;
            int isym = 0;
            for (int j = 0; j < nrotkm; ++j)
            {
                kvec_rot = kv.kvec_d[i] * kgmatrix[j];
                restrict_kpt(kvec_rot);
                for (int k = 0; k < nkstot_ibz; ++k)
                {
                    if (symm.equal(kvec_rot.x, kvec_d_ibz[k].x) && symm.equal(kvec_rot.y, kvec_d_ibz[k].y)
                        && symm.equal(kvec_rot.z, kvec_d_ibz[k].z))
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
            kv.kstars[exist_number].insert(std::make_pair(isym, kv.kvec_d[i]));
        }
    }
#endif

    // output in kpoints file
    std::stringstream ss;
    ss << " " << std::setw(40) << "nkstot"
       << " = " << kv.get_nkstot() << std::setw(66) << "ibzkpt" << std::endl;
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
    for (int i = 0; i < kv.get_nkstot(); ++i)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8d%12.8f%12.8f%12.8f\n",
                                 i + 1,
                                 kv.kvec_d[i].x,
                                 kv.kvec_d[i].y,
                                 kv.kvec_d[i].z,
                                 ibz_index[i] + 1,
                                 kvec_d_ibz[ibz_index[i]].x,
                                 kvec_d_ibz[ibz_index[i]].y,
                                 kvec_d_ibz[ibz_index[i]].z);
    }
    ss << table << std::endl;
    skpt = ss.str();
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nkstot_ibz", nkstot_ibz);

    table.clear();
    table += "\n K-POINTS REDUCTION ACCORDING TO SYMMETRY:\n";
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
    GlobalV::ofs_running << table << std::endl;

    // resize the kpoint container according to nkstot_ibz
    if (use_symm || kv.get_is_mp())
    {
        kv.update_use_ibz(nkstot_ibz, kvec_d_ibz, wk_ibz);
    }

    return;
}
} // namespace KVectorUtils
