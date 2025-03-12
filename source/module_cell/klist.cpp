#include "klist.h"

#include "module_base/formatter.h"
#include "module_base/memory.h"
#include "module_base/parallel_common.h"
#include "module_base/parallel_global.h"
#include "module_base/parallel_reduce.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/berryphase.h"
#include "module_parameter/parameter.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

void K_Vectors::cal_ik_global()
{
    const int my_pool = this->para_k.my_pool;
    this->ik2iktot.resize(this->nks);
#ifdef __MPI
    if(this->nspin == 2)
    {
        for (int ik = 0; ik < this->nks / 2; ++ik)
        {
            this->ik2iktot[ik] = this->para_k.startk_pool[my_pool] + ik;
            this->ik2iktot[ik + this->nks / 2] = this->nkstot / 2 + this->para_k.startk_pool[my_pool] + ik;
        }
    }
    else
    {
        for (int ik = 0; ik < this->nks; ++ik)
        {
            this->ik2iktot[ik] = this->para_k.startk_pool[my_pool] + ik;
        }
    }
#else
    for (int ik = 0; ik < this->nks; ++ik)
    {
        this->ik2iktot[ik] = ik;
    }
#endif

}

void K_Vectors::set(const UnitCell& ucell,
                    const ModuleSymmetry::Symmetry& symm,
                    const std::string& k_file_name,
                    const int& nspin_in,
                    const ModuleBase::Matrix3& reciprocal_vec,
                    const ModuleBase::Matrix3& latvec,
                    std::ofstream& ofs)
{
    ModuleBase::TITLE("K_Vectors", "set");

    ofs << "\n\n\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " | Setup K-points                                                     |" << std::endl;
    ofs << " | We setup the k-points according to input parameters.               |" << std::endl;
    ofs << " | The reduced k-points are set according to symmetry operations.     |" << std::endl;
    ofs << " | We treat the spin as another set of k-points.                      |" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    ofs << "\n\n\n\n";

    ofs << "\n SETUP K-POINTS" << std::endl;

    // (1) set nspin, read kpoints.
    this->nspin = nspin_in;
    ModuleBase::GlobalFunc::OUT(ofs, "nspin", nspin);

    if (this->nspin != 1 && this->nspin != 2 && this->nspin != 4)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::set", "Only available for nspin = 1 or 2 or 4");
    }

    this->nspin = (this->nspin == 4) ? 1 : this->nspin;

    // read KPT file and generate K-point grid
    bool read_succesfully = this->read_kpoints(ucell,k_file_name);
#ifdef __MPI
    Parallel_Common::bcast_bool(read_succesfully);
#endif
    if (!read_succesfully)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::set", "Something wrong while reading KPOINTS.");
    }

    // output kpoints file
    std::string skpt1;
    std::string skpt2;

    // (2)
    // only berry phase need all kpoints including time-reversal symmetry!
    // if symm_flag is not set, only time-reversal symmetry would be considered.
    if (!berryphase::berry_phase_flag && ModuleSymmetry::Symmetry::symm_flag != -1)
    {
        bool match = true;
        // calculate kpoints in IBZ and reduce kpoints according to symmetry
        this->ibz_kpoint(symm, ModuleSymmetry::Symmetry::symm_flag, skpt1, ucell, match);
#ifdef __MPI
        Parallel_Common::bcast_bool(match);
#endif
        if (!match)
        {
            std::cout << "Optimized lattice type of reciprocal lattice cannot match the optimized real lattice. "
                      << std::endl;
            std::cout << "It is often because the inaccuracy of lattice parameters in STRU." << std::endl;
            if (ModuleSymmetry::Symmetry::symm_autoclose)
            {
                ModuleBase::WARNING("K_Vectors::ibz_kpoint", "Automatically set symmetry to 0 and continue ...");
                std::cout << "Automatically set symmetry to 0 and continue ..." << std::endl;
                ModuleSymmetry::Symmetry::symm_flag = 0;
                match = true;
                this->ibz_kpoint(symm, ModuleSymmetry::Symmetry::symm_flag, skpt1, ucell, match);
            } else {
                ModuleBase::WARNING_QUIT("K_Vectors::ibz_kpoint",
                                         "Possible solutions: \n \
1. Refine the lattice parameters in STRU;\n \
2. Use a different`symmetry_prec`.  \n \
3. Close symemtry: set `symmetry` to 0 in INPUT. \n \
4. Set `symmetry_autoclose` to 1 in INPUT to automatically close symmetry when this error occurs.");
            }
        }
    }

    // (3)
    // Improve k point information

    // Complement the coordinates of k point
    this->set_both_kvec(reciprocal_vec, latvec, skpt2);

    if (GlobalV::MY_RANK == 0)
    {
        // output kpoints file
        std::stringstream skpt;
        skpt << PARAM.globalv.global_readin_dir << "kpoints";
        std::ofstream ofkpt(skpt.str().c_str()); // clear kpoints
        ofkpt << skpt2 << skpt1;
        ofkpt.close();
    }

    int deg = (nspin_in == 1) ? 2 : 1;
    // normalize k points weights according to nspin
    this->normalize_wk(deg);

    // It's very important in parallel case,
    // firstly do the mpi_k() and then
    // do set_kup_and_kdw()
    this->para_k.kinfo(nkstot,
                       GlobalV::KPAR,
                       GlobalV::MY_POOL,
                       GlobalV::RANK_IN_POOL,
                       GlobalV::NPROC,
                       nspin_in); // assign k points to several process pools
#ifdef __MPI
    // distribute K point data to the corresponding process
    this->mpi_k(); // 2008-4-29
#endif

    // set the k vectors for the up and down spin
    this->set_kup_and_kdw();

    // get ik2iktot
    this->cal_ik_global();

    this->print_klists(ofs);

    // std::cout << " NUMBER OF K-POINTS   : " << nkstot << std::endl;

#ifdef USE_PAW
    GlobalC::paw_cell.set_isk(nks, isk.data());
#endif

    return;
}

// 1.reset the size of the K-point container according to nspin and nkstot
// 2.reserve space for nspin>2 (symmetry)
void K_Vectors::renew(const int& kpoint_number)
{
    kvec_c.resize(kpoint_number);
    kvec_d.resize(kpoint_number);
    wk.resize(kpoint_number);
    isk.resize(kpoint_number);
    ngk.resize(kpoint_number);

    /*ModuleBase::Memory::record("KV::kvec_c",sizeof(double) * kpoint_number*3);
    ModuleBase::Memory::record("KV::kvec_d",sizeof(double) * kpoint_number*3);
    ModuleBase::Memory::record("KV::wk",sizeof(double) * kpoint_number*3);
    ModuleBase::Memory::record("KV::isk",sizeof(int) * kpoint_number*3);
    ModuleBase::Memory::record("KV::ngk",sizeof(int) * kpoint_number*3);*/

    return;
}

// Read the KPT file, which contains K-point coordinates, weights, and grid size information
// Generate K-point grid according to different parameters of the KPT file
bool K_Vectors::read_kpoints(const UnitCell& ucell,
                             const std::string& fn)
{
    ModuleBase::TITLE("K_Vectors", "read_kpoints");
    if (GlobalV::MY_RANK != 0)
    {
        return true;
    }

    // 1. Overwrite the KPT file and default K-point information if needed
    // mohan add 2010-09-04
    if (PARAM.globalv.gamma_only_local)
    {
        GlobalV::ofs_warning << " Auto generating k-points file: " << fn << std::endl;
        std::ofstream ofs(fn.c_str());
        ofs << "K_POINTS" << std::endl;
        ofs << "0" << std::endl;
        ofs << "Gamma" << std::endl;
        ofs << "1 1 1 0 0 0" << std::endl;
        ofs.close();
    }
    else if (PARAM.inp.kspacing[0] > 0.0)
    {
        if (PARAM.inp.kspacing[1] <= 0 || PARAM.inp.kspacing[2] <= 0)
        {
            ModuleBase::WARNING_QUIT("K_Vectors", "kspacing should > 0");
        };
        // number of K points = max(1,int(|bi|/KSPACING+1))
        ModuleBase::Matrix3 btmp = ucell.G;
        double b1 = sqrt(btmp.e11 * btmp.e11 + btmp.e12 * btmp.e12 + btmp.e13 * btmp.e13);
        double b2 = sqrt(btmp.e21 * btmp.e21 + btmp.e22 * btmp.e22 + btmp.e23 * btmp.e23);
        double b3 = sqrt(btmp.e31 * btmp.e31 + btmp.e32 * btmp.e32 + btmp.e33 * btmp.e33);
        int nk1
            = std::max(1, static_cast<int>(b1 * ModuleBase::TWO_PI / PARAM.inp.kspacing[0] / ucell.lat0 + 1));
        int nk2
            = std::max(1, static_cast<int>(b2 * ModuleBase::TWO_PI / PARAM.inp.kspacing[1] / ucell.lat0 + 1));
        int nk3
            = std::max(1, static_cast<int>(b3 * ModuleBase::TWO_PI / PARAM.inp.kspacing[2] / ucell.lat0 + 1));

        GlobalV::ofs_warning << " Generate k-points file according to KSPACING: " << fn << std::endl;
        std::ofstream ofs(fn.c_str());
        ofs << "K_POINTS" << std::endl;
        ofs << "0" << std::endl;
        ofs << "Gamma" << std::endl;
        ofs << nk1 << " " << nk2 << " " << nk3 << " 0 0 0" << std::endl;
        ofs.close();
    }

    // 2. Generate the K-point grid automatically according to the KPT file
    // 2.1 read the KPT file
    std::ifstream ifk(fn.c_str());
    if (!ifk)
    {
        GlobalV::ofs_warning << " Can't find File name : " << fn << std::endl;
        return false;
    }

    ifk >> std::setiosflags(std::ios::uppercase);

    ifk.clear();
    ifk.seekg(0);

    std::string word;
    std::string kword;

    int ierr = 0;

    ifk.rdstate();

    while (ifk.good())
    {
        ifk >> word;
        ifk.ignore(150, '\n'); // LiuXh add 20180416, fix bug in k-point file when the first line with comments
        if (word == "K_POINTS" || word == "KPOINTS" || word == "K")
        {
            ierr = 1;
            break;
        }

        ifk.rdstate();
    }

    if (ierr == 0)
    {
        GlobalV::ofs_warning << " symbol K_POINTS not found." << std::endl;
        return false;
    }

    // input k-points are in 2pi/a units
    ModuleBase::GlobalFunc::READ_VALUE(ifk, nkstot);

    this->k_nkstot = nkstot; // LiuXh add 20180619

    // std::cout << " nkstot = " << nkstot << std::endl;
    ModuleBase::GlobalFunc::READ_VALUE(ifk, kword);

    this->k_kword = kword; // LiuXh add 20180619

    // mohan update 2021-02-22
    const int max_kpoints = 100000;
    if (nkstot > max_kpoints)
    {
        GlobalV::ofs_warning << " nkstot > MAX_KPOINTS" << std::endl;
        return false;
    }

    // 2.2 Select different methods and generate K-point grid
    int k_type = 0;
    if (nkstot == 0) // nkstot==0, use monkhorst_pack. add by dwan
    {
        if (kword == "Gamma") // MP(Gamma)
        {
            is_mp = true;
            k_type = 0;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Input type of k points", "Monkhorst-Pack(Gamma)");
        }
        else if (kword == "Monkhorst-Pack" || kword == "MP" || kword == "mp")
        {
            is_mp = true;
            k_type = 1;
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Input type of k points", "Monkhorst-Pack");
        }
        else
        {
            GlobalV::ofs_warning << " Error: neither Gamma nor Monkhorst-Pack." << std::endl;
            return false;
        }

        ifk >> nmp[0] >> nmp[1] >> nmp[2];

        koffset[0] = 0;
        koffset[1] = 0;
        koffset[2] = 0;
        if (!(ifk >> koffset[0] >> koffset[1] >> koffset[2]))
        {
            ModuleBase::WARNING("K_Vectors::read_kpoints", "Missing k-point offsets in the k-points file.");
        }

        this->Monkhorst_Pack(nmp, koffset, k_type);
    }
    else if (nkstot > 0) // nkstot>0, the K-point information is clearly set
    {
        if (kword == "Cartesian" || kword == "C") // Cartesian coordinates
        {
            this->renew(nkstot * nspin); // mohan fix bug 2009-09-01
            for (int i = 0; i < nkstot; i++)
            {
                ifk >> kvec_c[i].x >> kvec_c[i].y >> kvec_c[i].z;
                ModuleBase::GlobalFunc::READ_VALUE(ifk, wk[i]);
            }

            this->kc_done = true;
        }
        else if (kword == "Direct" || kword == "D") // Direct coordinates
        {
            this->renew(nkstot * nspin); // mohan fix bug 2009-09-01
            for (int i = 0; i < nkstot; i++)
            {
                ifk >> kvec_d[i].x >> kvec_d[i].y >> kvec_d[i].z;
                ModuleBase::GlobalFunc::READ_VALUE(ifk, wk[i]);
            }
            this->kd_done = true;
        }
        else if (kword == "Line_Cartesian")
        {
            if (ModuleSymmetry::Symmetry::symm_flag == 1)
            {
                ModuleBase::WARNING("K_Vectors::read_kpoints",
                                    "Line mode of k-points is open, please set symmetry to 0 or -1.");
                return false;
            }

            interpolate_k_between(ifk, kvec_c);

            std::for_each(wk.begin(), wk.end(), [](double& d) { d = 1.0; });

            this->kc_done = true;
        }

        else if (kword == "Line_Direct" || kword == "L" || kword == "Line")
        {
            if (ModuleSymmetry::Symmetry::symm_flag == 1)
            {
                ModuleBase::WARNING("K_Vectors::read_kpoints",
                                    "Line mode of k-points is open, please set symmetry to 0 or -1.");
                return false;
            }

            interpolate_k_between(ifk, kvec_d);

            std::for_each(wk.begin(), wk.end(), [](double& d) { d = 1.0; });

            this->kd_done = true;
        }

        else
        {
            GlobalV::ofs_warning << " Error : neither Cartesian nor Direct kpoint." << std::endl;
            return false;
        }
    }

    this->nkstot_full = this->nks = this->nkstot;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nkstot", nkstot);
    return true;
} // END SUBROUTINE

void K_Vectors::interpolate_k_between(std::ifstream& ifk, std::vector<ModuleBase::Vector3<double>>& kvec)
{
    // how many special points.
    int nks_special = this->nkstot;

    // number of points to the next k points
    std::vector<int> nkl(nks_special, 0);

    // coordinates of special points.
    std::vector<ModuleBase::Vector3<double>> ks(nks_special);

    // recalculate nkstot.
    nkstot = 0;
    /* ISSUE#3482: to distinguish different kline segments */
    std::vector<int> kpt_segids;
    kl_segids.clear();
    kl_segids.shrink_to_fit();
    int kpt_segid = 0;
    for (int iks = 0; iks < nks_special; iks++)
    {
        ifk >> ks[iks].x;
        ifk >> ks[iks].y;
        ifk >> ks[iks].z;
        ModuleBase::GlobalFunc::READ_VALUE(ifk, nkl[iks]);

        assert(nkl[iks] >= 0);
        nkstot += nkl[iks];
        /* ISSUE#3482: to distinguish different kline segments */
        if ((nkl[iks] == 1) && (iks != (nks_special - 1))) {
            kpt_segid++;
        }
        kpt_segids.push_back(kpt_segid);
    }
    assert(nkl[nks_special - 1] == 1);

    // std::cout << " nkstot = " << nkstot << std::endl;
    this->renew(nkstot * nspin); // mohan fix bug 2009-09-01

    int count = 0;
    for (int iks = 1; iks < nks_special; iks++)
    {
        double dxs = (ks[iks].x - ks[iks - 1].x) / nkl[iks - 1];
        double dys = (ks[iks].y - ks[iks - 1].y) / nkl[iks - 1];
        double dzs = (ks[iks].z - ks[iks - 1].z) / nkl[iks - 1];
        for (int is = 0; is < nkl[iks - 1]; is++)
        {
            kvec[count].x = ks[iks - 1].x + is * dxs;
            kvec[count].y = ks[iks - 1].y + is * dys;
            kvec[count].z = ks[iks - 1].z + is * dzs;
            kl_segids.push_back(kpt_segids[iks - 1]); /* ISSUE#3482: to distinguish different kline segments */
            ++count;
        }
    }

    // deal with the last special k point.
    kvec[count].x = ks[nks_special - 1].x;
    kvec[count].y = ks[nks_special - 1].y;
    kvec[count].z = ks[nks_special - 1].z;
    kl_segids.push_back(kpt_segids[nks_special - 1]); /* ISSUE#3482: to distinguish different kline segments */
    ++count;

    assert(count == nkstot);
    assert(kl_segids.size() == nkstot); /* ISSUE#3482: to distinguish different kline segments */
}

double K_Vectors::Monkhorst_Pack_formula(const int& k_type, const double& offset, const int& n, const int& dim)
{
    double coordinate;
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

// add by dwan
void K_Vectors::Monkhorst_Pack(const int* nmp_in, const double* koffset_in, const int k_type)
{
    const int mpnx = nmp_in[0];
    const int mpny = nmp_in[1];
    const int mpnz = nmp_in[2];

    this->nkstot = mpnx * mpny * mpnz;
    // only can renew after nkstot is estimated.
    this->renew(nkstot * nspin); // mohan fix bug 2009-09-01

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

void K_Vectors::update_use_ibz(const int& nkstot_ibz,
                               const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                               const std::vector<double>& wk_ibz)
{
    if (GlobalV::MY_RANK != 0) {
        return;
    }
    ModuleBase::TITLE("K_Vectors", "update_use_ibz");
    assert(nkstot_ibz > 0);
    assert(nkstot_ibz <= kvec_d_ibz.size());
    // update nkstot
    this->nks = this->nkstot = nkstot_ibz;

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nkstot now", nkstot);

    this->kvec_d.resize(this->nkstot * nspin); // qianrui fix a bug 2021-7-13 for nspin=2 in set_kup_and_kdw()

    for (int i = 0; i < this->nkstot; ++i)
    {
        this->kvec_d[i] = kvec_d_ibz[i];

        // update weight.
        this->wk[i] = wk_ibz[i];
    }

    this->kd_done = true;
    this->kc_done = false;
    return;
}

void K_Vectors::ibz_kpoint(const ModuleSymmetry::Symmetry& symm,
                           bool use_symm,
                           std::string& skpt,
                           const UnitCell& ucell,
                           bool& match)
{
    if (GlobalV::MY_RANK != 0) {
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
    if (this->is_mp)
    {
        k_vec1 = ModuleBase::Vector3<double>(recip_vec1.x / nmp[0], recip_vec1.y / nmp[0], recip_vec1.z / nmp[0]);
        k_vec2 = ModuleBase::Vector3<double>(recip_vec2.x / nmp[1], recip_vec2.y / nmp[1], recip_vec2.z / nmp[1]);
        k_vec3 = ModuleBase::Vector3<double>(recip_vec3.x / nmp[2], recip_vec3.y / nmp[2], recip_vec3.z / nmp[2]);
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
        GlobalV::ofs_running << "(for reciprocal lattice: )" << std::endl;
          ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "BRAVAIS TYPE", recip_brav_type);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "BRAVAIS LATTICE NAME", recip_brav_name);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ibrav", recip_brav_type);

        // the map of bravis lattice from real to reciprocal space
        // for example, 3(fcc) in real space matches 2(bcc) in reciprocal space
        std::vector<int> ibrav_a2b{1, 3, 2, 4, 5, 6, 7, 8, 10, 9, 11, 12, 13, 14};
        // check if the reciprocal lattice is compatible with the real space lattice
        auto ibrav_match = [&](int ibrav_b) -> bool {
            const int& ibrav_a = symm.real_brav;
            if (ibrav_a < 1 || ibrav_a > 14) {
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
                              nullptr);
            GlobalV::ofs_running << "(for k-lattice: )" << std::endl;
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
            if (!match) {
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
    else if (is_mp) // only include for mp grid
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
    if (this->is_mp) {
        symm.gmatrix_convert(kgmatrix.data(), kkmatrix, nrotkm, ucell.G, k_vec);
    }
    // direct coordinates of k-points in k-lattice
    std::vector<ModuleBase::Vector3<double>> kvec_d_k(nkstot);
    if (this->is_mp) {
        for (int i = 0; i < nkstot; ++i) {
            kvec_d_k[i] = kvec_d[i] * ucell.G * k_vec.Inverse();
        }
    }

    // use operation : kgmatrix to find
    // the new set kvec_d : ir_kpt
    int nkstot_ibz = 0;

    assert(nkstot > 0);
    std::vector<ModuleBase::Vector3<double>> kvec_d_ibz(this->nkstot);
    std::vector<double> wk_ibz(this->nkstot); // ibz kpoint wk ,weight of k points
    std::vector<int> ibz2bz(this->nkstot);

    // nkstot is the total input k-points number.
    const double weight = 1.0 / static_cast<double>(nkstot);

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
        if (std::abs(kvec.x) < symm.epsilon) {
            kvec.x = 0.0;
        }
        if (std::abs(kvec.y) < symm.epsilon) {
            kvec.y = 0.0;
        }
        if (std::abs(kvec.z) < symm.epsilon) {
            kvec.z = 0.0;
        }
        return;
    };
    // for output in kpoints file
    int ibz_index[nkstot];
    // search in all k-poins.
    for (int i = 0; i < nkstot; ++i)
    {
        // restrict to [0, 1)
        restrict_kpt(kvec_d[i]);

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
                kvec_rot = kvec_d[i] * kgmatrix[j]; // wrong for total energy, but correct for nonlocal force.
                // kvec_rot = kgmatrix[j] * kvec_d[i]; //correct for total energy, but wrong for nonlocal force.
                restrict_kpt(kvec_rot);
                if (this->is_mp)
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
			//if it's a new ibz kpoint.
			//nkstot_ibz indicate the index of ibz kpoint.
            kvec_d_ibz[nkstot_ibz] = kvec_d[i];
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

            double kmol_new = kvec_d[i].norm2();
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
                kvec_d_ibz[exist_number] = kvec_d[i];
            }
        }
        //		BLOCK_HERE("check k point");
    }

    delete[] kkmatrix;

#ifdef __EXX
    // setup kstars according to the final (max-norm) kvec_d_ibz
    this->kstars.resize(nkstot_ibz);
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        for (int i = 0; i < nkstot; ++i)
        {
            int exist_number = -1;
            int isym = 0;
            for (int j = 0; j < nrotkm; ++j)
            {
                kvec_rot = kvec_d[i] * kgmatrix[j];
                restrict_kpt(kvec_rot);
                for (int k = 0; k < nkstot_ibz; ++k)
                {
                    if (symm.equal(kvec_rot.x, kvec_d_ibz[k].x) &&
                        symm.equal(kvec_rot.y, kvec_d_ibz[k].y) &&
                        symm.equal(kvec_rot.z, kvec_d_ibz[k].z))
                    {
                        isym = j;
                        exist_number = k;
                        break;
                    }
                }
                if (exist_number != -1) { break;
}
            }
            this->kstars[exist_number].insert(std::make_pair(isym, kvec_d[i]));
        }
    }
#endif

    // output in kpoints file
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
                                 this->kvec_d[i].x,
                                 this->kvec_d[i].y,
                                 this->kvec_d[i].z,
                                 ibz_index[i] + 1,
                                 kvec_d_ibz[ibz_index[i]].x,
                                 kvec_d_ibz[ibz_index[i]].y,
                                 kvec_d_ibz[ibz_index[i]].z);
    }
    ss << table << std::endl;
    skpt = ss.str();
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nkstot_ibz", nkstot_ibz);

    table.clear();
    table += "K-POINTS REDUCTION ACCORDING TO SYMMETRY\n";
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
    if (use_symm || is_mp)
    {
        this->update_use_ibz(nkstot_ibz, kvec_d_ibz, wk_ibz);
    }

    return;
}

// complement coordinates of k-points according to existing coordinates
// if cartesian coordinates are given, then direct coordinates are calculated
// if direct coordinates are given, then cartesian coordinates are calculated
void K_Vectors::set_both_kvec(const ModuleBase::Matrix3& G, const ModuleBase::Matrix3& R, std::string& skpt)
{

    if (PARAM.inp.final_scf) // LiuXh add 20180606
    {
        if (k_nkstot == 0)
        {
            kd_done = true;
            kc_done = false;
        }
        else
        {
            if (k_kword == "Cartesian" || k_kword == "C")
            {
                kc_done = true;
                kd_done = false;
            }
            else if (k_kword == "Direct" || k_kword == "D")
            {
                kd_done = true;
                kc_done = false;
            }
            else
            {
                GlobalV::ofs_warning << " Error : neither Cartesian nor Direct kpoint." << std::endl;
            }
        }
    }

    // set cartesian k vectors.
    if (!kc_done && kd_done)
    {
        for (int i = 0; i < nkstot; i++)
        {
            // wrong!!   kvec_c[i] = G * kvec_d[i];
            //  mohan fixed bug 2010-1-10
            if (std::abs(kvec_d[i].x) < 1.0e-10) {
                kvec_d[i].x = 0.0;
            }
            if (std::abs(kvec_d[i].y) < 1.0e-10) {
                kvec_d[i].y = 0.0;
            }
            if (std::abs(kvec_d[i].z) < 1.0e-10) {
                kvec_d[i].z = 0.0;
            }

            kvec_c[i] = kvec_d[i] * G;

            // mohan add2012-06-10
            if (std::abs(kvec_c[i].x) < 1.0e-10) {
                kvec_c[i].x = 0.0;
            }
            if (std::abs(kvec_c[i].y) < 1.0e-10) {
                kvec_c[i].y = 0.0;
            }
            if (std::abs(kvec_c[i].z) < 1.0e-10) {
                kvec_c[i].z = 0.0;
            }
        }
        kc_done = true;
    }

    // set direct k vectors
    else if (kc_done && !kd_done)
    {
        ModuleBase::Matrix3 RT = R.Transpose();
        for (int i = 0; i < nkstot; i++)
        {
            //			std::cout << " ik=" << i
            //				<< " kvec.x=" << kvec_c[i].x
            //				<< " kvec.y=" << kvec_c[i].y
            //				<< " kvec.z=" << kvec_c[i].z << std::endl;
            // wrong!            kvec_d[i] = RT * kvec_c[i];
            // mohan fixed bug 2011-03-07
            kvec_d[i] = kvec_c[i] * RT;
        }
        kd_done = true;
    }
    std::string table;
    table += "K-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < nkstot; i++)
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
           << " = " << nkstot << std::endl;
        ss << table << std::endl;
        skpt = ss.str();
    }
    return;
}

void K_Vectors::normalize_wk(const int& degspin)
{
    if (GlobalV::MY_RANK != 0) {
        return;
    }
    double sum = 0.0;

    for (int ik = 0; ik < nkstot; ik++)
    {
        sum += this->wk[ik];
    }
    assert(sum > 0.0);

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

#ifdef __MPI
void K_Vectors::mpi_k()
{
    ModuleBase::TITLE("K_Vectors", "mpi_k");

    Parallel_Common::bcast_bool(kc_done);

    Parallel_Common::bcast_bool(kd_done);

    Parallel_Common::bcast_int(nspin);

    Parallel_Common::bcast_int(nkstot);

    Parallel_Common::bcast_int(nkstot_full);

    Parallel_Common::bcast_int(nmp, 3);

    kl_segids.resize(nkstot);
    Parallel_Common::bcast_int(kl_segids.data(), nkstot);

    Parallel_Common::bcast_double(koffset, 3);

    this->nks = this->para_k.nks_pool[GlobalV::MY_POOL];

    GlobalV::ofs_running << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "k-point number in this process", nks);
    int nks_minimum = this->nks;

    Parallel_Reduce::gather_min_int_all(GlobalV::NPROC, nks_minimum);

    if (nks_minimum == 0)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::mpi_k()", " nks == 0, some processor have no k point!");
    }
    else
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "minimum distributed K point number", nks_minimum);
    }

    std::vector<int> isk_aux(nkstot);
    std::vector<double> wk_aux(nkstot);
    std::vector<double> kvec_c_aux(nkstot * 3);
    std::vector<double> kvec_d_aux(nkstot * 3);

    // collect and process in rank 0
    if (GlobalV::MY_RANK == 0)
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
        }
    }

    // broadcast k point data to all processors
    Parallel_Common::bcast_int(isk_aux.data(), nkstot);

    Parallel_Common::bcast_double(wk_aux.data(), nkstot);
    Parallel_Common::bcast_double(kvec_c_aux.data(), nkstot * 3);
    Parallel_Common::bcast_double(kvec_d_aux.data(), nkstot * 3);

    // process k point data in each processor
    this->renew(this->nks * this->nspin);

    // distribute
    int k_index = 0;

    for (int i = 0; i < nks; i++)
    {
        // 3 is because each k point has three value:kx, ky, kz
        k_index = i + this->para_k.startk_pool[GlobalV::MY_POOL];
        kvec_c[i].x = kvec_c_aux[k_index * 3];
        kvec_c[i].y = kvec_c_aux[k_index * 3 + 1];
        kvec_c[i].z = kvec_c_aux[k_index * 3 + 2];
        kvec_d[i].x = kvec_d_aux[k_index * 3];
        kvec_d[i].y = kvec_d_aux[k_index * 3 + 1];
        kvec_d[i].z = kvec_d_aux[k_index * 3 + 2];
        wk[i] = wk_aux[k_index];
        isk[i] = isk_aux[k_index];
    }

#ifdef __EXX
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {//bcast kstars
        this->kstars.resize(nkstot);
        for (int ikibz = 0;ikibz < nkstot;++ikibz)
        {
            int starsize = this->kstars[ikibz].size();
            Parallel_Common::bcast_int(starsize);
            GlobalV::ofs_running << "starsize: " << starsize << std::endl;
            auto ks = this->kstars[ikibz].begin();
            for (int ik = 0;ik < starsize;++ik)
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
                GlobalV::ofs_running << "isym: " << isym << " ks_vec: " << ks_vec.x << " " << ks_vec.y << " " << ks_vec.z << std::endl;
                if (GlobalV::MY_RANK != 0)
                {
                    kstars[ikibz].insert(std::make_pair(isym, ks_vec));
                }
            }
        }
    }
#endif
} // END SUBROUTINE
#endif

//----------------------------------------------------------
// This routine sets the k vectors for the up and down spin
//----------------------------------------------------------
// from set_kup_and_kdw.f90
void K_Vectors::set_kup_and_kdw()
{
    ModuleBase::TITLE("K_Vectors", "setup_kup_and_kdw");

    //=========================================================================
    // on output: the number of points is doubled and xk and wk in the
    // first (nks/2) positions correspond to up spin
    // those in the second (nks/2) ones correspond to down spin
    //=========================================================================
    switch (nspin)
    {
    case 1:

        for (int ik = 0; ik < nks; ik++)
        {
            this->isk[ik] = 0;
        }

        break;

    case 2:

        for (int ik = 0; ik < nks; ik++)
        {
            this->kvec_c[ik + nks] = kvec_c[ik];
            this->kvec_d[ik + nks] = kvec_d[ik];
            this->wk[ik + nks] = wk[ik];
            this->isk[ik] = 0;
            this->isk[ik + nks] = 1;
        }

        this->nks *= 2;
        this->nkstot *= 2;

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nks(nspin=2)", nks);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nkstot(nspin=2)", nkstot);
        break;
    case 4:

        for (int ik = 0; ik < nks; ik++)
        {
            this->isk[ik] = 0;
        }

        break;
    }

    return;
} // end subroutine set_kup_and_kdw

void K_Vectors::print_klists(std::ofstream& ofs)
{
    ModuleBase::TITLE("K_Vectors", "print_klists");

    if (nkstot < nks)
    {
        std::cout << "\n nkstot=" << nkstot;
        std::cout << "\n nks=" << nks;
        ModuleBase::WARNING_QUIT("print_klists", "nkstot < nks");
    }
    std::string table;
    table += "K-POINTS CARTESIAN COORDINATES\n";
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
    table += "K-POINTS DIRECT COORDINATES\n";
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

// LiuXh add a new function here,
// 20180515
void K_Vectors::set_after_vc(const int& nspin_in,
                             const ModuleBase::Matrix3& reciprocal_vec,
                             const ModuleBase::Matrix3& latvec)
{
    ModuleBase::TITLE("K_Vectors", "set_after_vc");

    GlobalV::ofs_running << "\n SETUP K-POINTS" << std::endl;
    this->nspin = nspin_in;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nspin", nspin);

    // set cartesian k vectors.
    kd_done = true;
    kc_done = false;
    if (!kc_done && kd_done)
    {
        for (int i = 0; i < nks; i++)
        {
            // wrong!!   kvec_c[i] = G * kvec_d[i];
            //  mohan fixed bug 2010-1-10
            if (std::abs(kvec_d[i].x) < 1.0e-10) {
                kvec_d[i].x = 0.0;
            }
            if (std::abs(kvec_d[i].y) < 1.0e-10) {
                kvec_d[i].y = 0.0;
            }
            if (std::abs(kvec_d[i].z) < 1.0e-10) {
                kvec_d[i].z = 0.0;
            }

            kvec_c[i] = kvec_d[i] * reciprocal_vec;

            // mohan add2012-06-10
            if (std::abs(kvec_c[i].x) < 1.0e-10) {
                kvec_c[i].x = 0.0;
            }
            if (std::abs(kvec_c[i].y) < 1.0e-10) {
                kvec_c[i].y = 0.0;
            }
            if (std::abs(kvec_c[i].z) < 1.0e-10) {
                kvec_c[i].z = 0.0;
            }
        }
        kc_done = true;
    }

    // set direct k vectors
    else if (kc_done && !kd_done)
    {
        ModuleBase::Matrix3 RT = latvec.Transpose();
        for (int i = 0; i < nks; i++)
        {
            //			std::cout << " ik=" << i
            //				<< " kvec.x=" << kvec_c[i].x
            //				<< " kvec.y=" << kvec_c[i].y
            //				<< " kvec.z=" << kvec_c[i].z << std::endl;
            // wrong!            kvec_d[i] = RT * kvec_c[i];
            // mohan fixed bug 2011-03-07
            kvec_d[i] = kvec_c[i] * RT;
        }
        kd_done = true;
    }
    std::string table;
    table += "K-POINTS DIRECT COORDINATES\n";
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
    GlobalV::ofs_running << table << std::endl;
    // this->set_both_kvec(reciprocal_vec, latvec);

    this->print_klists(GlobalV::ofs_running);

    return;
}
