#include "klist.h"

#include "klist_io.h"
#include "source_base/formatter.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"
#include "source_cell/module_symmetry/symmetry.h"

void K_Vectors::set(const UnitCell& ucell,
                    const ModuleSymmetry::Symmetry& symm,
                    const std::string& k_file_name,
                    const int& nspin_in,
                    const ModuleBase::Matrix3& reciprocal_vec,
                    const ModuleBase::Matrix3& latvec,
                    std::ofstream& ofs,
                    std::ofstream& ofs_warning,
                    const bool use_ibz,
                    const std::string& global_out_dir,
                    const bool gamma_only_local,
                    const double kspacing[3],
                    const std::string& kmesh_type,
                    const double koffset[3])
{
    ModuleBase::TITLE("K_Vectors", "set");

    ofs << "\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " |                       #Setup K-points#                             |" << std::endl;
    ofs << " | We setup the k-points according to input parameters.               |" << std::endl;
    ofs << " | The reduced k-points are set according to symmetry operations.     |" << std::endl;
    ofs << " | We treat the spin as another set of k-points.                      |" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    ofs << "\n";

    ofs << "\n SETUP K-POINTS" << std::endl;

    const std::string global_out_dir_ = global_out_dir;
    const bool gamma_only_local_ = gamma_only_local;
    const std::string kmesh_type_ = kmesh_type;
    const int my_rank = GlobalV::MY_RANK;
    const int my_pool = GlobalV::MY_POOL;

    // (1) print nspin, set the k-point spin multiplicity, read kpoints.
    ModuleBase::GlobalFunc::OUT(ofs, "nspin", nspin_in);

    if (nspin_in != 1 && nspin_in != 2 && nspin_in != 4)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::set", "Only available for nspin 1, 2 or 4");
    }

    // non-collinear (nspin=4) does not double the k-point list, so its
    // k-point spin multiplicity is the same as for the unpolarized case.
    this->spin_mult = (nspin_in == 4) ? 1 : nspin_in;

    bool read_succesfully = this->read_kpoints(ucell,
                                               k_file_name,
                                               gamma_only_local_,
                                               kspacing,
                                               kmesh_type_,
                                               koffset,
                                               ofs,
                                               ofs_warning,
                                               my_rank);
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

    // complement the Cartesian coordinates of the full k-point list
    KListIO::fill_full_kvec(this->kc_done,
                            this->kd_done,
                            this->nkstot_nospin,
                            reciprocal_vec,
                            this->kvec_c,
                            this->kvec_d,
                            this->kvec_c_full);


    // (2)
    // reduce kpoints to IBZ according to symmetry operations
    if (use_ibz)
    {
        bool match = true;
        // calculate kpoints in IBZ and reduce kpoints according to symmetry
        this->reduce_by_symmetry(ucell, symm, ModuleSymmetry::Symmetry::symm_flag, skpt1, match, my_rank, ofs);
#ifdef __MPI
        Parallel_Common::bcast_bool(match);
#endif
        if (!match)
        {
            this->handle_symmetry_mismatch(ucell, symm, skpt1, match, my_rank, ofs);
        }
    }

    // (3)
    // Improve k point information

    // Complement the coordinates of k point
    this->set_both_kvec(reciprocal_vec, latvec, skpt2, ofs, ofs_warning);

    if (my_rank == 0)
    {
        // output kpoints file
        std::stringstream skpt;
        skpt << global_out_dir_ << "KPT.info"; //mohan modified 20250325
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
                       my_pool,
                       GlobalV::RANK_IN_POOL,
                       GlobalV::NPROC,
                       nspin_in); // assign k points to several process pools
#ifdef __MPI
    // distribute K point data to the corresponding process
    this->mpi_k(ofs, my_rank, my_pool);
#endif

    // set the k vectors for the up and down spin
    this->set_kup_and_kdw(ofs);

    // initialize ibz_index
    this->ibz_index.resize(this->nkstot_nospin);
    for (int ik = 0; ik < this->nkstot_nospin; ik++)
    {
        this->ibz_index[ik] = ik;
    }
    
    // get ik2iktot: map local k indices to global indices in the pool
    KListIO::build_ik2iktot(this->para_k.my_pool,
                            this->para_k.startk_pool,
                            this->spin_mult,
                            this->nks,
                            this->nkstot,
                            this->ik2iktot);

    this->print_klists(ofs);

    // std::cout << " NUMBER OF K-POINTS   : " << nkstot << std::endl;

    return;
}

void K_Vectors::handle_symmetry_mismatch(const UnitCell& ucell,
                                         const ModuleSymmetry::Symmetry& symm,
                                         std::string& skpt,
                                         bool& match,
                                         const int my_rank,
                                         std::ofstream& ofs)
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
        this->reduce_by_symmetry(ucell, symm, ModuleSymmetry::Symmetry::symm_flag, skpt, match, my_rank, ofs);
    }
    else
    {
        ModuleBase::WARNING_QUIT("K_Vectors::ibz_kpoint",
                                 "Possible solutions: \n \
1. Refine the lattice parameters in STRU;\n \
2. Use a different`symmetry_prec`.  \n \
3. Close symemtry: set `symmetry` to 0 in INPUT. \n \
4. Set `symmetry_autoclose` to 1 in INPUT to automatically close symmetry when this error occurs.");
    }
}

// Resize the k-point containers to kpoint_number. The base class resizes
// kvec_c/kvec_d/kvec_c_full/wk/ngk; here we additionally resize isk.
// Callers pass nkstot * spin_mult so spin-polarized (nspin=2) runs have
// room for the up/down doubling done in set_kup_and_kdw().
void K_Vectors::renew(const int& kpoint_number)
{
    ReciprocalGrid::renew(kpoint_number);
    isk.resize(kpoint_number);

    return;
}

// Read the KPT file, which contains K-point coordinates, weights, and grid size information
// Generate K-point grid according to different parameters of the KPT file
bool K_Vectors::read_kpoints(const UnitCell& ucell,
                             const std::string& fn,
                             const bool gamma_only_local,
                             const double kspacing[3],
                             const std::string& kmesh_type,
                             const double koffset[3],
                             std::ofstream& ofs_running,
                             std::ofstream& ofs_warning,
                             const int my_rank)
{
    ModuleBase::TITLE("K_Vectors", "read_kpoints");
    if (my_rank != 0)
    {
        return true;
    }

    // 1. Overwrite the KPT file and default K-point information if needed
    // mohan add 2010-09-04
    KListIO::write_auto_kfile(ucell, fn, gamma_only_local, kspacing, kmesh_type, koffset, ofs_warning);

    // 2. Read the KPT file and build the k-point list
    return this->parse_kfile(fn, ofs_running, ofs_warning);
}

// 2. Generate the K-point grid automatically according to the KPT file
bool K_Vectors::parse_kfile(const std::string& fn, std::ofstream& ofs_running, std::ofstream& ofs_warning)
{
    // 2.1 read the KPT file
    std::ifstream ifk(fn.c_str());
    if (!ifk)
    {
        ofs_warning << " Can't find File name : " << fn << std::endl;
        return false;
    }

    ifk >> std::setiosflags(std::ios::uppercase);

    ifk.clear();
    ifk.seekg(0);

    std::string kword;

    if (!KListIO::find_kpoints_header(ifk))
    {
        ofs_warning << " symbol K_POINTS not found." << std::endl;
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
        ofs_warning << " nkstot > MAX_KPOINTS" << std::endl;
        return false;
    }

    // 2.2 Select different methods and generate K-point grid
    bool kpts_ok = true;
    if (nkstot == 0) // nkstot==0, use monkhorst_pack. add by dwan
    {
        kpts_ok = this->read_mp_mesh(ifk, kword, ofs_running, ofs_warning);
    }
    else if (nkstot > 0) // nkstot>0, the K-point information is clearly set
    {
        kpts_ok = this->read_listed_kpoints(ifk, kword, ofs_warning);
    }

    if (!kpts_ok)
    {
        return false;
    }

    this->nkstot_nospin = this->nks = this->nkstot;

    ModuleBase::GlobalFunc::OUT(ofs_running, "nkstot", nkstot);
    return true;
} // END SUBROUTINE

bool K_Vectors::read_mp_mesh(std::ifstream& ifk,
                             const std::string& kword,
                             std::ofstream& ofs_running,
                             std::ofstream& ofs_warning)
{
    int k_type = 0;
    if (kword == "Gamma") // MP(Gamma)
    {
        is_mp = true;
        k_type = 0;
        ModuleBase::GlobalFunc::OUT(ofs_running, "Input type of k points", "Monkhorst-Pack(Gamma)");
    }
    else if (kword == "Monkhorst-Pack" || kword == "MP" || kword == "mp")
    {
        is_mp = true;
        k_type = 1;
        ModuleBase::GlobalFunc::OUT(ofs_running, "Input type of k points", "Monkhorst-Pack");
    }
    else
    {
        ofs_warning << " Error: neither Gamma nor Monkhorst-Pack." << std::endl;
        return false;
    }

    ifk >> nmp[0] >> nmp[1] >> nmp[2];

    this->koffset[0] = 0;
    this->koffset[1] = 0;
    this->koffset[2] = 0;
    if (!(ifk >> this->koffset[0] >> this->koffset[1] >> this->koffset[2]))
    {
        ofs_warning << " K_Vectors::read_kpoints  warning : "
                    << "Missing k-point offsets in the k-points file." << std::endl;
    }

    this->Monkhorst_Pack(nmp, this->koffset, k_type);
    return true;
}

bool K_Vectors::read_listed_kpoints(std::ifstream& ifk, const std::string& kword, std::ofstream& ofs_warning)
{
    if (kword == "Cartesian" || kword == "C") // Cartesian coordinates
    {
        this->renew(nkstot * this->spin_mult); // mohan fix bug 2009-09-01
        KListIO::read_kpt_list(ifk, nkstot, this->kvec_c, this->wk);
        this->kc_done = true;
        return true;
    }
    if (kword == "Direct" || kword == "D") // Direct coordinates
    {
        this->renew(nkstot * this->spin_mult); // mohan fix bug 2009-09-01
        KListIO::read_kpt_list(ifk, nkstot, this->kvec_d, this->wk);
        this->kd_done = true;
        return true;
    }
    if (kword == "Line_Cartesian")
    {
        return this->setup_line_kpoints(ifk, this->kvec_c, true, ofs_warning);
    }
    if (kword == "Line_Direct" || kword == "L" || kword == "Line")
    {
        return this->setup_line_kpoints(ifk, this->kvec_d, false, ofs_warning);
    }

    ofs_warning << " Error : neither Cartesian nor Direct kpoint." << std::endl;
    return false;
}

bool K_Vectors::setup_line_kpoints(std::ifstream& ifk,
                                   std::vector<ModuleBase::Vector3<double>>& kvec,
                                   const bool cartesian,
                                   std::ofstream& ofs_warning)
{
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        ofs_warning << " K_Vectors::read_kpoints  warning : "
                    << "Line mode of k-points is open, please set symmetry to 0 or -1."
                    << std::endl;
        return false;
    }

    this->interpolate_k_between(ifk, kvec);

    std::for_each(this->wk.begin(), this->wk.end(), [](double& d) { d = 1.0; });

    if (cartesian)
    {
        this->kc_done = true;
    }
    else
    {
        this->kd_done = true;
    }
    return true;
}

void K_Vectors::interpolate_k_between(std::ifstream& ifk, std::vector<ModuleBase::Vector3<double>>& kvec)
{
    // Thin wrapper: the interpolation itself is the this-free KListIO::interp_line;
    // here we only size the member containers and copy the results back.
    const KListIO::LineK line = KListIO::interp_line(ifk, this->nkstot);

    this->nkstot = line.nks_total;
    this->renew(this->nkstot * this->spin_mult); // mohan fix bug 2009-09-01

    for (int i = 0; i < this->nkstot; i++)
    {
        kvec[i] = line.kpts[i];
    }
    this->kl_segids = line.segids; /* ISSUE#3482: to distinguish different kline segments */
}

void K_Vectors::update_use_ibz(const int& nkstot_ibz,
                               const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                               const std::vector<double>& wk_ibz,
                               std::ofstream& ofs_running,
                               const int my_rank)
{
    if (my_rank != 0) {
        return;
    }
    ModuleBase::TITLE("K_Vectors", "update_use_ibz");
    assert(nkstot_ibz > 0);
    assert(nkstot_ibz <= kvec_d_ibz.size());
    // update nkstot
    this->nks = this->nkstot = nkstot_ibz;

    ModuleBase::GlobalFunc::OUT(ofs_running, "nkstot now", nkstot);

    // qianrui fix a bug 2021-7-13: size for the spin_mult=2 doubling in set_kup_and_kdw()
    this->kvec_d.resize(this->nkstot * this->spin_mult);

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

//----------------------------------------------------------
// This routine sets the k vectors for the up and down spin
//----------------------------------------------------------
// from set_kup_and_kdw.f90
void K_Vectors::set_kup_and_kdw(std::ofstream& ofs_running)
{
    ModuleBase::TITLE("K_Vectors", "setup_kup_and_kdw");

    KListIO::expand_spin_kpoints(this->spin_mult,
                                 this->kvec_c,
                                 this->kvec_d,
                                 this->wk,
                                 this->isk,
                                 this->nks,
                                 this->nkstot);

    if (this->spin_mult == 2)
    {
        ModuleBase::GlobalFunc::OUT(ofs_running, "nks(nspin=2)", this->nks);
        ModuleBase::GlobalFunc::OUT(ofs_running, "nkstot(nspin=2)", this->nkstot);
    }

    return;
} // end subroutine set_kup_and_kdw

void K_Vectors::reduce_by_symmetry(const UnitCell& ucell,
                                   const ModuleSymmetry::Symmetry& symm,
                                   bool use_symm,
                                   std::string& skpt,
                                   bool& match,
                                   const int my_rank,
                                   std::ofstream& ofs_running)
{
    if (my_rank != 0)
    {
        return;
    }
    ModuleBase::TITLE("K_Vectors", "reduce_by_symmetry");

    //===============================================
    // search in all space group operations
    // if the operations does not already included
    // inverse operation, double it.
    //===============================================
    std::vector<ModuleBase::Matrix3> kgmatrix(48 * 2);

    ModuleBase::Matrix3 k_vec;
    int nrotkm = 0;
    if (!this->build_star_ops(ucell, symm, use_symm, k_vec, kgmatrix, nrotkm))
    {
        match = false;
        return;
    }
    if (nrotkm == 0)
    {
        return;
    }

    // append time-reversal-related operations (Theta*g for magnetic
    // nspin=4; -g otherwise unless inversion is already present)
    nrotkm = KListIO::append_time_reversal_ops(symm, kgmatrix, nrotkm);

    // convert kgmatrix to k-lattice
    std::vector<ModuleBase::Matrix3> kkmatrix(nrotkm);
    if (this->get_is_mp())
    {
        symm.gmatrix_convert(kgmatrix.data(), kkmatrix.data(), nrotkm, ucell.G, k_vec);
    }

    // use operation : kgmatrix to find
    // the new set kvec_d : ir_kpt
    std::vector<ModuleBase::Vector3<double>> kvec_d_ibz;
    std::vector<double> wk_ibz;
    std::vector<int> ibz2bz;
    this->reduce_ibz(kgmatrix.data(),
                     nrotkm,
                     ucell.G,
                     k_vec,
                     kkmatrix.data(),
                     symm.epsilon,
                     kvec_d_ibz,
                     wk_ibz,
                     this->ibz_index,
                     ibz2bz);
    const int nkstot_ibz = kvec_d_ibz.size();

#ifdef __EXX
    // setup kstars according to the final (max-norm) kvec_d_ibz
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        KListIO::build_kstars(this->kvec_d,
                              kgmatrix,
                              nrotkm,
                              kvec_d_ibz,
                              symm.epsilon,
                              [&symm](double a, double b) { return symm.equal(a, b); },
                              this->kstars);
    }
#endif

    // output in kpoints file
    skpt = KListIO::ibz_kpt_table(this->nkstot, this->kvec_d, this->ibz_index, kvec_d_ibz);
    ModuleBase::GlobalFunc::OUT(ofs_running, "Number of irreducible k-points", nkstot_ibz);

    ofs_running << KListIO::ibz_wk_table(nkstot_ibz, kvec_d_ibz, wk_ibz, ibz2bz) << std::endl;

    // resize the kpoint container according to nkstot_ibz
    if (use_symm || this->get_is_mp())
    {
        this->update_use_ibz(nkstot_ibz, kvec_d_ibz, wk_ibz, ofs_running, my_rank);
    }

    return;
}

void K_Vectors::set_after_vc(const ModuleBase::Matrix3& G, std::ofstream& ofs_running)
{
    ofs_running << "\n SETUP K-POINTS" << std::endl;

    // set cartesian k vectors.
    this->kvec_d2c(G);

    std::string table;
    table += "K-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < this->nks; i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 this->kvec_d[i].x,
                                 this->kvec_d[i].y,
                                 this->kvec_d[i].z,
                                 this->wk[i]);
    }
    ofs_running << table << std::endl;

    this->kd_done = true;
    this->kc_done = true;

    this->print_klists(ofs_running);
}

#ifdef __MPI
void K_Vectors::mpi_k(std::ofstream& ofs_running, const int my_rank, const int my_pool)
{
    ModuleBase::TITLE("K_Vectors", "mpi_k");

    Parallel_Common::bcast_bool(this->kc_done);

    Parallel_Common::bcast_bool(this->kd_done);

    Parallel_Common::bcast_int(this->spin_mult);

    Parallel_Common::bcast_int(this->nkstot);

    Parallel_Common::bcast_int(this->nkstot_nospin);

    Parallel_Common::bcast_int(this->nmp, 3);

    this->kl_segids.resize(this->nkstot);
    Parallel_Common::bcast_int(this->kl_segids.data(), this->nkstot);

    Parallel_Common::bcast_double(this->koffset, 3);

    this->nks = this->para_k.nks_pool[my_pool];

    ofs_running << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_running, "Number of k-points in this process", this->nks);
    int nks_minimum = this->nks;

    Parallel_Reduce::reduce_min(nks_minimum);

    if (nks_minimum == 0)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::mpi_k()", " nks == 0, some processor have no k points!");
    }
    else
    {
        ModuleBase::GlobalFunc::OUT(ofs_running, "Minimum distributed k-point number", nks_minimum);
    }

    std::vector<int> isk_aux(this->nkstot);
    std::vector<double> wk_aux(this->nkstot);
    std::vector<double> kvec_c_aux(this->nkstot * 3);
    std::vector<double> kvec_d_aux(this->nkstot * 3);
    std::vector<double> kvec_c_full_aux(this->nkstot_nospin * 3);

    // collect and process in rank 0
    if (my_rank == 0)
    {
        KListIO::pack_kpts(this->isk,
                           this->wk,
                           this->kvec_c,
                           this->kvec_d,
                           this->kvec_c_full,
                           this->nkstot,
                           isk_aux,
                           wk_aux,
                           kvec_c_aux,
                           kvec_d_aux,
                           kvec_c_full_aux);
    }

    // broadcast k point data to all processors
    Parallel_Common::bcast_int(isk_aux.data(), this->nkstot);

    Parallel_Common::bcast_double(wk_aux.data(), this->nkstot);
    Parallel_Common::bcast_double(kvec_c_aux.data(), this->nkstot * 3);
    Parallel_Common::bcast_double(kvec_d_aux.data(), this->nkstot * 3);
    Parallel_Common::bcast_double(kvec_c_full_aux.data(), this->nkstot_nospin * 3);

    // process k point data in each processor
    this->renew(this->nks * this->spin_mult);

    // distribute
    KListIO::unpack_kpts(isk_aux,
                         wk_aux,
                         kvec_c_aux,
                         kvec_d_aux,
                         kvec_c_full_aux,
                         this->nks,
                         this->para_k.startk_pool[my_pool],
                         this->isk,
                         this->wk,
                         this->kvec_c,
                         this->kvec_d,
                         this->kvec_c_full);

#ifdef __EXX
    // bcast kstars (rank 0 holds the filled maps; other ranks rebuild them)
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        KListIO::bcast_kstars(this->kstars, this->nkstot, my_rank);
    }
#endif
} // END SUBROUTINE mpi_k
#endif
