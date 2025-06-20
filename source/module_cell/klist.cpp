#include "klist.h"

#include "k_vector_utils.h"
#include "source_base/formatter.h"
#include "source_base/memory.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/berryphase.h"
#include "module_parameter/parameter.h"

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

    ofs << "\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " | Setup K-points                                                     |" << std::endl;
    ofs << " | We setup the k-points according to input parameters.               |" << std::endl;
    ofs << " | The reduced k-points are set according to symmetry operations.     |" << std::endl;
    ofs << " | We treat the spin as another set of k-points.                      |" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    ofs << "\n\n";

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
        KVectorUtils::kvec_ibz_kpoint(*this, symm, ModuleSymmetry::Symmetry::symm_flag, skpt1, ucell, match);
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
                KVectorUtils::kvec_ibz_kpoint(*this, symm, ModuleSymmetry::Symmetry::symm_flag, skpt1, ucell, match);
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
//    this->set_both_kvec(reciprocal_vec, latvec, skpt2);
    KVectorUtils::set_both_kvec(*this, reciprocal_vec, latvec, skpt2);

    if (GlobalV::MY_RANK == 0)
    {
        // output kpoints file
        std::stringstream skpt;
        skpt << PARAM.globalv.global_out_dir << "KPT.info"; //mohan modified 20250325
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
    KVectorUtils::kvec_mpi_k(*this);
#endif

    // set the k vectors for the up and down spin
    this->set_kup_and_kdw();

    // get ik2iktot
    this->cal_ik_global();

    KVectorUtils::print_klists(*this, ofs);

    // std::cout << " NUMBER OF K-POINTS   : " << nkstot << std::endl;

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