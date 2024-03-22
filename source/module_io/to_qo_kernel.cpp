#include "module_io/to_qo.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_base/ylm.h"
#include "module_base/libm/libm.h"
#ifdef __MPI
#include "module_base/parallel_common.h"
#endif
// how define QO task, how create QO instance
toQO::toQO(const std::string& qo_basis, 
           const std::vector<std::string>& strategies,
           const double& qo_thr,
           const std::vector<double>& screening_coeffs)
{
    // totally the same as what defined in INPUT
    // qo_switch_ = 1 // certainly, this constructor will only be called when qo_switch_ == 1
    qo_basis_ = qo_basis;
    strategies_ = strategies;
    qo_thr_ = qo_thr;
    screening_coeffs_ = screening_coeffs;
}

toQO::~toQO() {} // we dont use any new or malloc, so no need to delete or free

// initialize means to initialize in actual program environment, therefore most of program
// -dependent vairables import here.
void toQO::initialize(const std::string& out_dir,
                      const std::string& pseudo_dir,
                      const std::string& orbital_dir,
                      const UnitCell* p_ucell,
                      const std::vector<ModuleBase::Vector3<double>>& kvecs_d,
                      std::ofstream& ofs_running,
                      const int& rank,
                      const int& nranks)
{
    // print parameter settings for QO
    if(rank == 0)
    {
        std::string init_info = "\nQuasiatomic orbital analysis activated.\n";
        init_info += "Parameters settings check:\n";
        init_info += "qo_basis: " + qo_basis_ + "\n";
        init_info += "qo_thr: " + std::to_string(qo_thr_) + "\n";
        init_info += "qo_strategies: ";
        for(auto s: strategies_) init_info += s + " ";
        init_info += "\n";
        init_info += "Output directory: " + out_dir + "\n";
        init_info += "Pseudopotential directory: " + pseudo_dir + "\n";
        init_info += "Numerical atomic orbital directory: " + orbital_dir + "\n";
        init_info += "Number of kpoints: " + std::to_string(kvecs_d.size()) + "\n";
        init_info += "Parallelized on " + std::to_string(nranks) + " MPI processes\n";
        printf("%s", init_info.c_str());
    }
#ifdef __MPI
    // this MPI_Barrier is to ensure the above information is printed before the following
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // initialize the variables defining I/O
    out_dir_ = out_dir;
    pseudo_dir_ = pseudo_dir;
    orbital_dir_ = orbital_dir;

    // initialize the variables defining parallelism
    iproc_ = rank;
    nprocs_ = nranks;

    // initialize the variables defining the structure
    // review what are needed by toQO module from UnitCell
    // 
    read_structures(p_ucell, kvecs_d, iproc_, nprocs_);

    // build orbitals, for all processes
    // orbitals building operation is to (for example) read all orbital files and save to
    // radial collection data structure. These kinds of information is irrelevant to the
    // structure simulated. Once two atoms' one certain overlap integral is needed, get
    // orbitals from radial collection by type-l-izeta, so that get the radial function.
    // Then specify the angular part Ylm l and m, also l' and m', and the correct distance.
    // Then the overlap integral is calculated.
    // build two-center overlap calculator
    overlap_calculator_ = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    // build the numerical atomic orbital basis
    build_nao(ntype_, orbital_dir_, p_ucell_->orbital_fn, iproc_);
    // build another atomic orbital
    build_ao(ntype_, pseudo_dir_, p_ucell_->pseudo_fn, screening_coeffs_, qo_thr_, ofs_running, iproc_);

    // neighbor list search, based on built RadialCollection(s)
    scan_supercell(iproc_, nprocs_);

    // build grids, for all processes
    double rcut_max = std::max(nao_->rcut_max(), ao_->rcut_max());
    int ngrid = int(rcut_max / 0.01) + 1;
    double cutoff = 2.0*rcut_max;
    nao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    ao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    overlap_calculator_->tabulate(*ao_, *nao_, 'S', ngrid, cutoff);

    // prepare for Ylm, if this is not called, the Ylm will not be available and always
    // return 0.0
    ModuleBase::Ylm::set_coefficients();

    // based on built RadialCollection(s), allocate memory for overlap integrals
    allocate_ovlp(true); allocate_ovlp(false);
}

void toQO::build_nao(const int ntype, 
                     const std::string orbital_dir,
                     const std::string* const orbital_fn,
                     const int rank)
{
    // build the numerical atomic orbital basis
    ModuleBase::SphericalBesselTransformer sbt;
    nao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    int ntype_ = ntype;
#ifdef __MPI
    Parallel_Common::bcast_int(ntype_);
#endif
    std::string* orbital_fn_ = new std::string[ntype_];
    if(rank == 0)
    {
        for(int it = 0; it < ntype_; it++)
        {
            orbital_fn_[it] = orbital_dir + orbital_fn[it];
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_string(orbital_fn_, ntype_);
#endif
    // this method is done by RANK-0, then broadcast to other processes.
    nao_->build(ntype_, orbital_fn_, 'o');
    nao_->set_transformer(sbt);
    // indexing, all processes can do together. It is also cumbersome to broadcast a std::map of std::tuple
    // by ABACUS self-built MPI broadcast function, so the following indexing is done by all processes,
    // directly.
    radialcollection_indexing(*nao_, na_, false, index_nao_, rindex_nao_);
    nphi_ = index_nao_.size();

    delete[] orbital_fn_;
    if(rank == 0)
    {
        std::string nao_build_info = "toQO::build_nao: built numerical atomic orbitals for calculating QO overlap integrals\n";
        nao_build_info += "Number of columns in QO_ovlp_*.dat: " + std::to_string(nphi_) + "\n";
        nao_build_info += "Orbitals arrange in sequence of (it, ia, l, zeta, m), m in order of 0, 1, -1, 2, -2, ...\n";
        printf("%s", nao_build_info.c_str());
    }
}

bool toQO::orbital_filter_out(const int& itype,
                              const int& l,
                              const int& izeta)
{
    // true = filter out, false = not filter out = keep
    // this function works for RadialCollection, to select the orbitals of interest
    std::vector<std::string> l2symbol = {"s", "p", "d", "f", "g"}; // seems enough
    if(qo_basis_ == "pswfc")
    {
        // for pswfc, what supports is specifying the name of subshell layer, like
        // qo_strategy sp spdf
        // means for the first atom type, use s and p orbitals, for the second, use
        // s, p, d, and f orbitals
        // default is `all` for all types, and for each type, all orbitals are used
        if(strategies_[itype] == "all") return false;
        else if(l >= l2symbol.size()) return true;
        else if(strategies_[itype].find_first_of(l2symbol[l]) != std::string::npos) return false;
        else return true;
    }
    else if(qo_basis_ == "szv")
    {
        // use two individual logic branch allows them have different orbital filtering logic,
        // although presently they are almost the same
        if(izeta != 0) return true; // filter out
        else if(strategies_[itype] == "all") return false; // keep
        else if(l >= l2symbol.size()) return true; // filter out
        else if(strategies_[itype].find_first_of(l2symbol[l]) != std::string::npos) return false; // keep
        else return true; // filter out
    }
    else return false;
}

void toQO::build_hydrogen(const int ntype, 
                          const double* const charges, 
                          const bool slater_screening,
                          const int* const nmax,
                          const double qo_thr,
                          const int rank)
{
    ModuleBase::SphericalBesselTransformer sbt;
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    // for this method, all processes CAN do together, there is no apparent conflict between processes
    ao_->build(ntype, charges, slater_screening, nmax, symbols_.data(), qo_thr, strategies_.data(), rank);
    ao_->set_transformer(sbt);
    // indexing, all processes can do together. It is also cumbersome to broadcast a std::map of std::tuple
    // by ABACUS self-built MPI broadcast function, so the following indexing is done by all processes,
    // directly.
    radialcollection_indexing(*ao_, na_, true, index_ao_, rindex_ao_);
    nchi_ = index_ao_.size();
}

void toQO::build_pswfc(const int ntype, 
                       const std::string pseudo_dir,
                       const std::string* const pspot_fn, 
                       const double* const screening_coeffs,
                       const double qo_thr,
                       const int rank)
{
    ModuleBase::SphericalBesselTransformer sbt;
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    int ntype_ = ntype;
#ifdef __MPI
    Parallel_Common::bcast_int(ntype_);
#endif
    std::string* pspot_fn_ = new std::string[ntype_];
    for(int it = 0; it < ntype; it++)
    {
        pspot_fn_[it] = pseudo_dir + pspot_fn[it];
    }
#ifdef __MPI
    Parallel_Common::bcast_string(pspot_fn_, ntype_);
#endif
    // for this method, all processes MIGHT NOT do together, because of possible conflict of reading files
    // in the following build function, the file reading is done by RANK-0, then broadcast to other processes
    // this MPI strategy is done by refactoring PswfcRadials instance. For details, see the impelementation
    // of build function of PswfcRadials: source/module_basis/module_nao/pswfc_radials.cpp
    ao_->build(ntype, pspot_fn_, screening_coeffs, qo_thr, rank);
    ao_->set_transformer(sbt);
    // indexing, all processes can do together. It is also cumbersome to broadcast a std::map of std::tuple
    // by ABACUS self-built MPI broadcast function, so the following indexing is done by all processes,
    // directly.
    radialcollection_indexing(*ao_, na_, true, index_ao_, rindex_ao_);
    nchi_ = index_ao_.size();
    delete[] pspot_fn_;
}

void toQO::build_szv()
{
    // build the numerical atomic orbital basis
    ModuleBase::SphericalBesselTransformer sbt;
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection(*nao_));
    ao_->set_transformer(sbt);
    // indexing, all processes can do together. It is also cumbersome to broadcast a std::map of std::tuple
    // by ABACUS self-built MPI broadcast function, so the following indexing is done by all processes,
    // directly.
    radialcollection_indexing(*ao_, na_, true, index_ao_, rindex_ao_);
    nchi_ = index_ao_.size();
}

void toQO::build_ao(const int ntype, 
                    const std::string pseudo_dir,
                    const std::string* const pspot_fn,
                    const std::vector<double> screening_coeffs,
                    const double qo_thr,
                    const std::ofstream& ofs_running,
                    const int rank)
{
    if(qo_basis_ == "hydrogen")
    {
        bool with_slater_screening = std::find_if(screening_coeffs.begin(), screening_coeffs.end(), 
            [](double sc) { return sc > 1e-10; }) != screening_coeffs.end();
        build_hydrogen(ntype_,                  /// ntype
                       charges_.data(),         /// charges
                       with_slater_screening,   /// slater_screening
                       nmax_.data(),            /// nmax
                       qo_thr,                  /// qo_thr
                       rank);                   /// rank
    }
    else if(qo_basis_ == "pswfc")
    {
        build_pswfc(ntype_,                     /// ntype
                    pseudo_dir,                 /// pseudo_dir
                    pspot_fn,                   /// pspot_fn
                    screening_coeffs.data(),    /// screening_coeffs
                    qo_thr,                     /// qo_thr
                    rank);                      /// rank
    }
    else if(qo_basis_ == "szv") build_szv();
    if(rank == 0)
    {
        std::string ao_build_info = "toQO::build_ao: built atomic orbitals for calculating QO overlap integrals\n";
        ao_build_info += "Atom-centered orbital to project is: " + qo_basis_ + "\n";
        ao_build_info += "Number of rows in QO_ovlp_*.dat: " + std::to_string(nchi_) + "\n";
        ao_build_info += "Orbitals arrange in sequence of (it, ia, l, zeta, m), m in order of 0, 1, -1, 2, -2, ...\n";
        printf("%s", ao_build_info.c_str());
    }
}

void toQO::calculate_ovlpR(const int iR)
{
    assert (rindex_ao_.size() == nchi_);
    assert (rindex_nao_.size() == nphi_);
    for(int irow = 0; irow < nchi_; irow++)
    {
        //         it,  ia,  li,  izeta, mi
        std::tuple<int, int, int, int, int> orb1 = rindex_ao_[irow];
        int it = std::get<0>(orb1);
        int ia = std::get<1>(orb1);
        int li = std::get<2>(orb1);
        int izeta = std::get<3>(orb1);
        int mi = std::get<4>(orb1);
        for(int icol = 0; icol < nphi_; icol++)
        {
            //         jt,  ja,  lj,  jzeta, mj
            std::tuple<int, int, int, int, int> orb2 = rindex_nao_[icol];
            int jt = std::get<0>(orb2);
            int ja = std::get<1>(orb2);
            int lj = std::get<2>(orb2);
            int jzeta = std::get<3>(orb2);
            int mj = std::get<4>(orb2);
            ModuleBase::Vector3<double> rij = p_ucell_->atoms[jt].tau[ja] - p_ucell_->atoms[it].tau[ia];
            // there is waste here, but for easy to understand, I don't optimize it.
            ModuleBase::Vector3<int> R = supercells_[iR];
            ModuleBase::Vector3<double> Rij;
            Rij.x = rij.x + double(R.x) * p_ucell_->a1.x 
                          + double(R.y) * p_ucell_->a2.x 
                          + double(R.z) * p_ucell_->a3.x;
            Rij.y = rij.y + double(R.x) * p_ucell_->a1.y 
                          + double(R.y) * p_ucell_->a2.y 
                          + double(R.z) * p_ucell_->a3.y;
            Rij.z = rij.z + double(R.x) * p_ucell_->a1.z 
                          + double(R.y) * p_ucell_->a2.z 
                          + double(R.z) * p_ucell_->a3.z;
            Rij *= p_ucell_->lat0; // convert to Bohr
            overlap_calculator_->calculate(
                it, li, izeta, mi,
                jt, lj, jzeta, mj,
                Rij, &ovlpR_[irow*nphi_+icol]
            );
        }
    }
}

void toQO::calculate_ovlpk(int ik)
{
    // On the parallelization of toQO module
    // 2024-03-12, kirk0830
    // Let's plant trees!
    // present parallelization strategy will brings a little bit of burden to the disk IO
    // another disk-friendly parallelization strategy is to save all two-center integrals
    // on one single process and broadcast to other processes. This is more memory-consuming,
    // but less disk IO.

    // for the first run, need to calculate all two-center integrals in realspace.
    // the calculation of two-center integrals in realspace is also parallelized.
    // but to avoid high frequency of sending and receiving messages, use file to
    // store the results of two-center integrals in realspace. Then read the files
    // and calculate the two-center integrals in k-space.
    if(ik == iks_[0])
    {
        for(auto iR: iRs_)
        {
            calculate_ovlpR(iR);
            write_ovlp<double>(out_dir_, ovlpR_, nchi_, nphi_, true, iR);
        }
    }
    // then read and calculate all ks corresponding two-center integrals on own process
    // there is a way to avoid the conflict of reading files between different processes,
    // that is, when proc1 is reading R1, let proc2 read R2, and so on. After the first
    // round of reading, then proc1 reads R2, proc2 reads R3, and so on. This way, everytime
    // before a new-reading, should all processes be synchronized.
    for(int iR = 0; iR < nR_tot_; iR++)
    {
        int barrier_iR = (iR + iproc_) % nR_tot_;
#ifdef __MPI
        // the following MPI_Barrier ensures two unexpected things, the first has been mentioned
        // above, the second is that, avoid the case for one process has finished calculation
        // of integral in realspace and ready for calculating in kspace but other processes are
        // still calculating in realspace, if that happens, the faster process can not find all
        // the files needed, which are, all QO_ovlpR_*.dat files.
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        // ik == -1 corresponds to the case of those processes with less kpoints than others
        if(ik != -1) read_ovlp(out_dir_, nchi_, nphi_, true, barrier_iR);
        if(ik != -1) append_ovlpR_eiRk(ik, barrier_iR);
    }
}

void toQO::calculate()
{
    for(auto ik: iks_)
    {
#ifdef __MPI
        // we must enforce the process to be highly synchronized and as much as possible near the REAL
        // kpoint parallelism because need to avoid the possible conflict of reading files between 
        // different processes
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        // while for zero_out overlap, it can be not so strictly synchronized
        zero_out_ovlps(false);
        calculate_ovlpk(ik);
        if(ik != -1) write_ovlp<std::complex<double>>(out_dir_, ovlpk_, nchi_, nphi_, false, ik);
    }
#ifdef __MPI
    // once the calculation of S(k) is finished, prone to delete all QO_ovlpR_*.dat files. But the most
    // important thing is to ensure all processes have finished the calculation of S(k), so MPI_Barrier.
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // delete all QO_ovlpR_* files
    if(iproc_ == 0)
    {
        for(int iR = 0; iR < nR_tot_; iR++)
        {
            std::string filename = out_dir_ + "/QO_ovlpR_" + std::to_string(iR) + ".dat";
            std::remove(filename.c_str());
        }
        printf("toQO::calculate: calculation of S(k) done, run /tools/qo/postprocess.py to do representation transform.\n");
    }
}

void toQO::append_ovlpR_eiRk(int ik, int iR)
{
    // calculate sum S(R)*eiRk = S(k)
    ModuleBase::Vector3<double> R(double(supercells_[iR].x), double(supercells_[iR].y), double(supercells_[iR].z));
    double arg = (kvecs_d_[ik] * R) * ModuleBase::TWO_PI;
    double sinp, cosp;
    ModuleBase::libm::sincos(arg, &sinp, &cosp);
    std::complex<double> phase = std::complex<double>(cosp, sinp);
    // add all values of ovlpR_ to ovlpk_ with multiplication of phase
    for(int i = 0; i < nchi_ * nphi_; i++) ovlpk_[i] += ovlpR_[i] * phase;
}

void toQO::allocate_ovlp(const bool& is_R)
{
    if(is_R) ovlpR_.resize(nchi_ * nphi_, 0.0);
    else ovlpk_.resize(nchi_ * nphi_, std::complex<double>(0.0, 0.0));
}

void toQO::deallocate_ovlp(const bool& is_R)
{
    if(is_R) {ovlpR_.clear(); ovlpR_.shrink_to_fit();}
    else {ovlpk_.clear(); ovlpk_.shrink_to_fit();}
}

void toQO::zero_out_ovlps(const bool& is_R)
{
    if(is_R) std::fill(ovlpR_.begin(), ovlpR_.end(), 0.0);
    else std::fill(ovlpk_.begin(), ovlpk_.end(), std::complex<double>(0.0, 0.0));
}

void toQO::radialcollection_indexing(const RadialCollection& radcol,
                                     const std::vector<int>& natoms,
                                     const bool& with_filter,
                                     std::map<std::tuple<int,int,int,int,int>,int>& index_map,
                                     std::map<int,std::tuple<int,int,int,int,int>>& index_map_reverse)
{
    // in RadialCollection, radials are stored type by type and actually not relevant with exact atom index,
    // so the number of atom of each type is external information.
    // the map should be like: (itype, iatom, l, izeta, m) -> index and the reverse one is index -> (itype, iatom, l, izeta, m)
    int index = 0;
    for(int itype = 0; itype < radcol.ntype(); itype++)
    {
        for(int iatom = 0; iatom < natoms[itype]; iatom++)
        {
            for(int l = 0; l <= radcol.lmax(itype); l++)
            {
                std::vector<int> ms;
                for(int m_abs = 0; m_abs <= l; m_abs++)
                {
                    ms.push_back(m_abs);
                    if(m_abs != 0) ms.push_back(-m_abs);
                }
                for(int izeta = 0; izeta < radcol.nzeta(itype, l); izeta++)
                {
                    // usually, the orbital is distinguished by it, l and zeta, the ia and m are not
                    // commonly used.
                    if(orbital_filter_out(itype, l, izeta)&&with_filter) continue;
                    for(int m: ms)
                    {
                        index_map[std::make_tuple(itype, iatom, l, izeta, m)] = index;
                        index_map_reverse[index] = std::make_tuple(itype, iatom, l, izeta, m);
                        index++;
                    }
                }
            }
        }
    }
}

// template function definition
// 20240310 make compatible with R space matrices
template <typename T>
void toQO::write_ovlp(const std::string& dir,
                      const std::vector<T> &ovlp, 
                      const int& nrows,
                      const int& ncols,
                      const bool& is_R,
                      const int& i)
{
    std::string filename = is_R? "QO_ovlpR_" + std::to_string(i) + ".dat": "QO_ovlp_" + std::to_string(i) + ".dat";
    std::ofstream ofs(dir + filename);
    if(!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("toQO::write_ovlp", "can not open file: " + filename);
    }
    if(is_R)
    {
        ofs << "SUPERCELL_COORDINATE: " << std::setw(5) << std::right << supercells_[i].x << " "
                                        << std::setw(5) << std::right << supercells_[i].y << " "
                                        << std::setw(5) << std::right << supercells_[i].z << std::endl;
    }
    else
    {
        ofs << "KPOINT_COORDINATE: " << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[i].x << " "
                                     << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[i].y << " "
                                     << std::setw(22) << std::setprecision(14) << std::right << std::scientific << kvecs_d_[i].z << std::endl;
    }
    for(int irow = 0; irow < nrows; irow++)
    {
        for(int icol = 0; icol < ncols; icol++)
        {
            ofs << std::setw(22) << std::setprecision(14) << std::right << std::scientific << ovlp[irow*ncols+icol] << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}
// explicit instantiation
template void toQO::write_ovlp<double>(const std::string& dir, const std::vector<double>& ovlp, const int& nrows, const int& ncols,
                                       const bool& is_R, const int& ik);
template void toQO::write_ovlp<std::complex<double>>(const std::string& dir, const std::vector<std::complex<double>>& ovlp, 
                                                     const int& nrows, const int& ncols, const bool& is_R, const int& ik);
// a free function to convert string storing C++ std::complex to std::complex
// format: (real,imag), both part in scientific format
std::complex<double> str2complex(const std::string& str)
{
    std::string real_str, imag_str;
    int i = 1; // skip '('
    while(str[i] != ',') real_str += str[i]; i++;
    i++; // skip ','
    while(str[i] != ')') imag_str += str[i]; i++;
    return std::complex<double>(std::stod(real_str), std::stod(imag_str));
}
// complete I/O of QO module
void toQO::read_ovlp(const std::string& dir,
                     const int& nrows,
                     const int& ncols,
                     const bool& is_R,
                     const int& ik)
{
    zero_out_ovlps(is_R); // clear the ovlp vector before reading
    assert (nrows * ncols == nchi_ * nphi_);
    std::string filename = is_R? "QO_ovlpR_" + std::to_string(ik) + ".dat": "QO_ovlp_" + std::to_string(ik) + ".dat";
    std::ifstream ifs(dir + "/" + filename);
    if(!ifs.is_open())
    {
        ModuleBase::WARNING_QUIT("toQO::read_ovlp", "can not open file: " + filename);
    }
    // read header
    std::string line;
    std::getline(ifs, line);
    // read ovlp values
    int inum = 0;
    while(ifs.good())
    {
        if(is_R)
        {
            double val;
            ifs >> val; inum++;
            if(inum <= nchi_ * nphi_) ovlpR_[inum-1] = val;
            else break;
        }
        else
        {
            std::string val_str;
            ifs >> val_str; inum++;
            if(inum <= nchi_ * nphi_) ovlpk_[inum-1] = str2complex(val_str);
            else break;
        }
    }
}
