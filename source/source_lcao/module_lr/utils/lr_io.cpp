#include "lr_io.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"
#include <algorithm>
#include <cassert>
#include <dirent.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace LR_IO {

void parse_band_out_file(const std::string& in_dir, int& nbands_file, int& nk_file, int& nspin_file, int& nocc_file)
{
    std::string file = in_dir + "band_out";
    std::ifstream ifs(file);
    if (!ifs) throw std::runtime_error(file + " not found");
    std::string tmp, line;
    double occ;
    int nocc_count = 0;

    ifs >> nk_file >> nspin_file >> nbands_file;
    std::cout << "band_out: nk = " << nk_file << std::endl;
    std::cout << "band_out: nspin = " << nspin_file << std::endl;
    std::cout << "band_out: nbands = " << nbands_file << std::endl;
    for (int i = 0; i < 4; ++i) {std::getline(ifs, tmp); } //skip 4 lines

    while (ifs.peek() != EOF)
    {
        std::getline(ifs, line);
        std::istringstream iss(line);

        iss >> tmp >> occ;
        if (occ > 0.1) nocc_count++;
        else if (occ < 0.1) break;
    }
    nocc_file = nocc_count;
    std::cout << "nocc in band_out: " << nocc_file << std::endl;

    if (PARAM.inp.bse_use_fine_kgrid)
    {
        file = in_dir + "band_kpath_info";
        std::ifstream ifs(file);
        if (!ifs) throw std::runtime_error(file + " not found");
        std::string tmp;
        ifs >> tmp >> nbands_file >> nspin_file >> nk_file;
        std::cout << "band_kpath_info: nk = " << nk_file << std::endl;
        std::cout << "band_kpath_info: nspin = " << nspin_file << std::endl;
        std::cout << "band_kpath_info: nbands = " << nbands_file << std::endl;
        std::cout << "BSE will use fine kgrid." << std::endl;
        ifs.close();
    }
}

std::vector<double> read_energy_qp(const int nocc,
                                   const int nvirt,
                                   const std::string& in_dir,
                                   int& ncore,
                                   const int nk,
                                   const int nspin_tmp,
                                   const int nspin_file)
{
    const std::string file = in_dir + "energy_qp";
    std::cout << "in read_energy_qp, nbands(nocc+nvir): " << (nocc+nvirt) << std::endl;
    std::vector<double> eig_info( nspin_tmp * nk * (nocc + nvirt) * 3 ); // occ, eig_ks, eig_gw
    std::ifstream ifs_gw (file);
    if (!ifs_gw) throw std::runtime_error(file + " not found");
    std::string temp;
    int read_ik;
    double occ, eig_ks, eig_gw;

    for (int is =0; is < nspin_file; ++is){
        for (int ik = 0; ik < nk; ++ik){
            for (int i = 0;i < 2;++i) { std::getline(ifs_gw, temp); } // skip the first 2 lines
            ifs_gw >> temp >> read_ik ;
            // std::cout << "ik: " << ik <<" is:" << is << std::endl;
            assert(ik == (read_ik-1));
            int ivirt = 0;
            std::getline(ifs_gw, temp); // skip the interval line
            std::getline(ifs_gw, temp); // skip the interval line
            std::vector<double> ks_temps;
            std::vector<double> gw_temps;
            std::vector<double> occ_temps;
            while (ifs_gw.peek() != '-')
            {
                std::getline(ifs_gw, temp);
                std::istringstream iss(temp);
                iss >> temp >> occ >> eig_ks >> eig_gw;
                ks_temps.push_back(eig_ks * 2); // Ha to Ry
                gw_temps.push_back(eig_gw * 2); // Ha to Ry
                occ_temps.push_back(occ);
                if (occ < 0.1) { ivirt++;}
                if (ivirt == nvirt) { break; }
            }
            ncore = gw_temps.size() - nocc - nvirt;
            for (int ib = 0;ib < nocc + nvirt;++ib)
            {   
                int ikstep = (ik + is * nk) * (nocc + nvirt);
                eig_info[(ikstep + ib)*3] = occ_temps[ncore + ib];
                eig_info[(ikstep + ib)*3 + 1] = ks_temps[ncore + ib];
                eig_info[(ikstep + ib)*3 + 2] = gw_temps[ncore + ib];
                // std::cout <<"GW_info: ik=" << std::setw(5) << ik << " ib=" << std::setw(5) << ib
                //         << std::setw(9) << eig_info[(ikstep + ib)*3] << std::setw(11)
                //         << eig_info[(ikstep + ib)*3 + 1] << std::setw(11)
                //         << eig_info[(ikstep + ib)*3 + 2] << std::endl; //check
            }
            while (ifs_gw.peek() != '-' && ifs_gw.peek() != EOF)
            {
                std::getline(ifs_gw, temp); // skip the virtual bands to next k-point
            }            
        }
    }
    if (nspin_file == 1 && nspin_tmp == 2) {
        std::cout << "duplicate the spin channel since the gw file only has one spin channel" << std::endl;
        int spin_block = nk * (nocc + nvirt) * 3;
        assert(eig_info.size() == 2 * spin_block);
        std::copy_n(eig_info.data(), spin_block, eig_info.data() + spin_block);
    }        
    ifs_gw.close();
    std::cout << "Finish read gw, ncore=" << ncore << std::endl;
    return eig_info;
}

std::vector<double> read_energy_qp_from_band_files(const K_Vectors& kv,
                                                   const int nocc,
                                                   const int nvirt,
                                                   int& ncore,
                                                   const std::string& in_dir,
                                                   const int nk,
                                                   const int nspin_tmp,
                                                   const int nspin_file)
{
    const std::string ks_prefix = in_dir + "KS_band_spin_";
    const std::string gw_prefix = in_dir + "GW_band_spin_";
    std::cout << "in read_energy_qp_from_band_files, nbands(nocc+nvir): " << (nocc+nvirt) << std::endl;
    std::vector<double> eig_info( nspin_tmp * nk * (nocc + nvirt) * 3 ); // occ, eig_ks, eig_gw
    for (int is =0; is < nspin_file; ++is)
    {
        std::string ks_file = ks_prefix + std::to_string(is + 1) + ".dat";
        std::string gw_file = gw_prefix + std::to_string(is + 1) + ".dat";
        std::ifstream ifs_ks (ks_file);
        if (!ifs_ks) throw std::runtime_error(ks_file + " not found");
        std::ifstream ifs_gw (gw_file);
        if (!ifs_gw) throw std::runtime_error(gw_file + " not found");
        std::string temp;
        int read_ik;
        double occ, eig_ks, eig_gw;
        double kx, ky, kz;
        for (int ik = 0; ik < nk; ++ik)
        {
            ifs_ks >> read_ik >> kx >> ky >> kz;
            // std::cout << "ik: " << ik <<" is:" << is << std::endl;
            double thread = 1.e-6;
            assert(ik == (read_ik-1));
            assert(std::abs(kx - kv.kvec_d[ik + is * nk].x) < thread);
            assert(std::abs(ky - kv.kvec_d[ik + is * nk].y) < thread);
            assert(std::abs(kz - kv.kvec_d[ik + is * nk].z) < thread);
            ifs_gw >> read_ik >> kx >> ky >> kz;
            assert(ik == (read_ik-1));
            assert(std::abs(kx - kv.kvec_d[ik + is * nk].x) < thread);
            assert(std::abs(ky - kv.kvec_d[ik + is * nk].y) < thread);
            assert(std::abs(kz - kv.kvec_d[ik + is * nk].z) < thread);
            int ivirt = 0;
            std::vector<double> ks_temps;
            std::vector<double> gw_temps;
            std::vector<double> occ_temps;

            std::getline(ifs_ks, temp);
            std::istringstream ks_line(temp);
            std::getline(ifs_gw, temp);
            std::istringstream gw_line(temp);
            while (true)
            {
                ks_line >> occ >> eig_ks;
                gw_line >> occ >> eig_gw;
                ks_temps.push_back(eig_ks / ModuleBase::Ry_to_eV);
                gw_temps.push_back(eig_gw / ModuleBase::Ry_to_eV);
                occ_temps.push_back(occ);
                if (occ*nk < 0.1) { ivirt++;}
                if (ivirt == nvirt) { break; }                
            }

            ncore = gw_temps.size() - nocc - nvirt;
            for (int ib = 0;ib < nocc + nvirt;++ib)
            {   
                int ikstep = (ik + is * nk) * (nocc + nvirt);
                eig_info[(ikstep + ib)*3] = occ_temps[ncore + ib];
                eig_info[(ikstep + ib)*3 + 1] = ks_temps[ncore + ib];
                eig_info[(ikstep + ib)*3 + 2] = gw_temps[ncore + ib];
                // std::cout <<"GW_info: ib=" << std::setw(5) << ib
                //         << std::setw(9) << eig_info[(ikstep + ib)*3] << std::setw(11)
                //         << eig_info[(ikstep + ib)*3 + 1] << std::setw(11)
                //         << eig_info[(ikstep + ib)*3 + 2] << std::endl; //check
            }          
        }
        ifs_ks.close();
        ifs_gw.close();
    }
    if (nspin_file == 1 && nspin_tmp == 2) {
        std::cout << "duplicate the spin channel since the gw file only has one spin channel" << std::endl;
        int spin_block = nk * (nocc + nvirt) * 3;
        assert(eig_info.size() == 2 * spin_block);
        std::copy_n(eig_info.data(), spin_block, eig_info.data() + spin_block);
    }
    std::cout << "Finish read gw, ncore=" << ncore << std::endl;
    return eig_info;
}

template <typename TK>
void read_librpa_eigenvectors(psi::Psi<TK>& wfc_ks,
                              psi::Psi<TK>& wfc_ks_global,
                              const std::string& in_dir,
                              const int ncore,
                              const int nbands_file,
                              const int nspin_tmp,
                              const int nspin_file,
                              const int my_rank,
                              Parallel_Orbitals& pmat)
{
    int nbands = pmat.get_wfc_global_nbands();// nbands = nocc + nvirt
    int nbasis = pmat.get_wfc_global_nbasis();
    assert(nbands == wfc_ks_global.get_nbands());
    assert(nbasis == wfc_ks_global.get_nbasis());
    assert((ncore + nbands) <= nbands_file);
    const size_t nk = PARAM.inp.nspin == 2 ? wfc_ks.get_nk() / 2 : wfc_ks.get_nk();

    if (my_rank == 0)
    {
        struct dirent *ptr;
        DIR *dir;
        dir = opendir(in_dir.c_str());
        std::vector<bool> readen_k(nk, false);

        while ((ptr = readdir(dir)) != NULL){// read all the files in the directory
            std::string fm(ptr->d_name);
            if (fm.find("KS_eigenvector") == 0)// find file KS_eigenvectorXXX
            {
                //std::cout << "found librpa_eigenvector file:" << fm << std::endl;
                std::ifstream file_librpa_ks(in_dir + fm);
                std::string tmp;
                while (file_librpa_ks.peek() != EOF)
                {
                    int ik;
                    file_librpa_ks >> ik;
                    ik = ik - 1;
                    assert(readen_k[ik] == false);
                    file_librpa_ks >> std::ws; // skip the blank and '\n' to get the next content
                    for (int iw = 0; iw < nbasis; ++iw) {
                        for (int ib = 0; ib < nbands_file; ++ib) {
                            for (int is = 0; is < nspin_file; ++is) {
                                if (ib >= ncore && ib< (ncore+nbands)) {
                                    LR_IO::read_one_data(file_librpa_ks, wfc_ks_global(ik+is*nk, ib-ncore, iw));
                                    file_librpa_ks >> std::ws;
                                }
                                else {
                                    std::getline(file_librpa_ks, tmp); //skip the useless bands
                                }
                            }
                        }
                    }
                    if (nspin_tmp == 2 && nspin_file == 1) {
                        std::copy_n(&wfc_ks_global(ik,0,0), nbands*nbasis, &wfc_ks_global(ik + nk,0,0));
                    }
                    readen_k[ik] = true;
                }
                file_librpa_ks.close();
            }
        }
        closedir(dir);
        for(int ik = 0; ik < nk; ++ik) {
            if (!readen_k[ik])
                throw std::runtime_error("librpa_eigenvector file not found for k-point " + std::to_string(ik+1));
        }
        
    }// end of if (my_rank == 0); next MPI_comm to other ranks

    // change wfc_ks_global phase to make arg(<psi(k)|psi(k'=0)>) = 0
    std::cout << "Do phase correction" << std::endl;
    for (int iks = 0; iks < wfc_ks.get_nk(); ++iks){
        if (my_rank == 0) {
            // test: output wfc            
            // std::cout << "wfc_gs_read_from_librpa for iks:" << iks << std::endl;
            // for (int ib = 0; ib < nbands; ++ib)
            // {
            //     std::cout << "band " << ib << ": ";
            //     for (int iw = 0; iw < nbasis; ++iw)
            //     {
            //         std::cout << wfc_ks_global(iks, ib, iw) << "  ";
            //     }
            //     std::cout << std::endl;
            // }
            if (iks != 0)
            {
                for(int ib = 0; ib < nbands; ++ib)
                {
                    TK phase = LR_Util::inner_product(&wfc_ks_global(iks,ib,0), &wfc_ks_global(0,ib,0), nbasis);
                    phase = phase / std::abs(phase);
                    for (int iw = 0; iw < nbasis; ++iw)
                    {
                        wfc_ks_global(iks, ib, iw) *= phase;
                    }
                    // TK test_phase = LR_Util::inner_product(&wfc_ks_global(iks,ib,0), &wfc_ks_global(0,ib,0), nbasis);
                    // std::cout << "After phase correction, iks, ib, phase: " << iks << " " << ib << " " << test_phase << std::endl;
                }
            }
        }

        wfc_ks_global.fix_k(iks);
        wfc_ks.fix_k(iks);
#ifdef __MPI
        MPI_Bcast(wfc_ks_global.get_pointer(), nbands * nbasis, LR_Util::MPIType<TK>::value(), 0, MPI_COMM_WORLD);
        Parallel_2D pv_glb;
        pv_glb.set(nbasis, nbands, std::max(nbasis, nbands), pmat.blacs_ctxt);
        Cpxgemr2d(nbasis, nbands, wfc_ks_global.get_pointer(), 1, 1, pv_glb.desc,
                    wfc_ks.get_pointer(), 1, 1, const_cast<int*>(pmat.desc_wfc)/*nbasis×nbands*/,
                    pv_glb.blacs_ctxt);
#else
        BlasConnector::copy(nbands*nlocal, wfc_ks_global.get_pointer(), 1, wfc_ks.get_pointer(), 1);
#endif
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read librpa eigenvectors.");
}

template <typename TK>
void read_librpa_eigenvectors_from_band_files(psi::Psi<TK>& wfc_ks,
                                              psi::Psi<TK>& wfc_ks_global,
                                              const std::string& in_dir,
                                              const int ncore,
                                              const int nbands_file,
                                              const int nspin_tmp,
                                              const int nspin_file,
                                              const int my_rank,
                                              Parallel_Orbitals& pmat)
{
    int nbands = pmat.get_wfc_global_nbands();// nbands = nocc + nvirt
    int nbasis = pmat.get_wfc_global_nbasis();
    assert(nbands == wfc_ks_global.get_nbands());
    assert(nbasis == wfc_ks_global.get_nbasis());
    assert((ncore + nbands) <= nbands_file);
    const size_t nk = PARAM.inp.nspin == 2 ? wfc_ks.get_nk() / 2 : wfc_ks.get_nk();
    
    if (my_rank == 0)
    {
        for (int ik = 0; ik < nk; ++ik)
        {
            std::stringstream ss;
            ss << in_dir << "band_KS_eigenvector_k_" << std::setfill('0') << std::setw(5) << ik + 1 << ".txt";
            std::ifstream infile(ss.str(), std::ios::binary);
            if (!infile)
            {
                throw std::runtime_error("Error: Cannot open file " + ss.str());
            }
            size_t total_complex = static_cast<size_t>(nspin_file * nbands_file * nbasis);
            size_t total_doubles = total_complex * 2;

            std::vector<double> double_buffer(total_doubles);
            infile.read(reinterpret_cast<char *>(double_buffer.data()), total_doubles * sizeof(double));
            if (infile.gcount() != static_cast<ptrdiff_t>(total_doubles * sizeof(double)))
            {
                throw std::runtime_error("Error: failed to read " + ss.str());
            }
            for (int is = 0; is < nspin_file; ++is) {
                for (int ib = ncore; ib < (ncore+nbands); ++ib) {
                    for (int iw = 0; iw < nbasis; ++iw) {
                        int index = 2 * (is * nbands_file * nbasis + ib * nbasis + iw);
                        LR_IO::read_one_data(double_buffer[index], double_buffer[index+1], wfc_ks_global(ik+is*nk, ib-ncore, iw));
                    }
                }
            }
            if (nspin_tmp == 2 && nspin_file == 1) {
                std::copy_n(&wfc_ks_global(ik,0,0), nbands*nbasis, &wfc_ks_global(ik + nk,0,0));
            }
            infile.close();
        }
    }// end of if (my_rank == 0); next MPI_comm to other ranks

    // change wfc_ks_global phase to make arg(<psi(k)|psi(k'=0)>) = 0
    std::cout << "Do phase correction" << std::endl;
    for (int iks = 0; iks < wfc_ks.get_nk(); ++iks){
        if (my_rank == 0) {
            // test: output wfc            
            // std::cout << "wfc_gs_read_from_librpa for iks:" << iks << std::endl;
            // for (int ib = 0; ib < nbands; ++ib)
            // {
            //     std::cout << "band " << ib << ": ";
            //     for (int iw = 0; iw < nbasis; ++iw)
            //     {
            //         std::cout << wfc_ks_global(iks, ib, iw) << "  ";
            //     }
            //     std::cout << std::endl;
            // }
            if (iks != 0)
            {
                for(int ib = 0; ib < nbands; ++ib)
                {
                    TK phase = LR_Util::inner_product(&wfc_ks_global(iks,ib,0), &wfc_ks_global(0,ib,0), nbasis);
                    phase = phase / std::abs(phase);
                    for (int iw = 0; iw < nbasis; ++iw)
                    {
                        wfc_ks_global(iks, ib, iw) *= phase;
                    }
                    // TK test_phase = LR_Util::inner_product(&wfc_ks_global(iks,ib,0), &wfc_ks_global(0,ib,0), nbasis);
                    // std::cout << "After phase correction, iks, ib, phase: " << iks << " " << ib << " " << test_phase << std::endl;
                }
            }
        }

        wfc_ks_global.fix_k(iks);
        wfc_ks.fix_k(iks);
#ifdef __MPI
        MPI_Bcast(wfc_ks_global.get_pointer(), nbands * nbasis, LR_Util::MPIType<TK>::value(), 0, MPI_COMM_WORLD);
        Parallel_2D pv_glb;
        pv_glb.set(nbasis, nbands, std::max(nbasis, nbands), pmat.blacs_ctxt);
        Cpxgemr2d(nbasis, nbands, wfc_ks_global.get_pointer(), 1, 1, pv_glb.desc,
                    wfc_ks.get_pointer(), 1, 1, const_cast<int*>(pmat.desc_wfc)/*nbasis×nbands*/,
                    pv_glb.blacs_ctxt);
#else
        BlasConnector::copy(nbands*nlocal, wfc_ks_global.get_pointer(), 1, wfc_ks.get_pointer(), 1);
#endif
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read librpa eigenvectors.");
}

template void read_librpa_eigenvectors<double>(
    psi::Psi<double>& wfc_ks, psi::Psi<double>& wfc_ks_global,
    const std::string& in_dir, const int ncore, const int nbands_file,
    const int nspin_tmp, const int nspin_file, const int my_rank, Parallel_Orbitals& pmat);
template void read_librpa_eigenvectors<std::complex<double>>(
    psi::Psi<std::complex<double>>& wfc_ks, psi::Psi<std::complex<double>>& wfc_ks_global,
    const std::string& in_dir, const int ncore, const int nbands_file,
    const int nspin_tmp, const int nspin_file, const int my_rank, Parallel_Orbitals& pmat);

template void read_librpa_eigenvectors_from_band_files<double>(
    psi::Psi<double>& wfc_ks, psi::Psi<double>& wfc_ks_global,
    const std::string& in_dir, const int ncore, const int nbands_file,
    const int nspin_tmp, const int nspin_file, const int my_rank, Parallel_Orbitals& pmat);
template void read_librpa_eigenvectors_from_band_files<std::complex<double>>(
    psi::Psi<std::complex<double>>& wfc_ks, psi::Psi<std::complex<double>>& wfc_ks_global,
    const std::string& in_dir, const int ncore, const int nbands_file,
    const int nspin_tmp, const int nspin_file, const int my_rank, Parallel_Orbitals& pmat);

// ================= LRI related functions =========================

#ifdef __EXX
template <typename TCs, typename TVs> // only for blocking by atom pairs
TLRI<TVs> read_coulomb_mat_k(const std::string& in_dir, const TLRI<TCs>& Cs, LR_IO::RI_kRlist& kRlist)
{
    struct dirent *ptr;
    DIR *dir;
    bool has_unshrinked = false;
    dir = opendir(in_dir.c_str());
    if (!dir)
    {
        throw std::runtime_error("Cannot open directory: " + in_dir);
    }
    while ((ptr = readdir(dir)) != nullptr)
    {
        std::string fm(ptr->d_name);
        if (fm.find("coulomb_unshrinked_cut_") == 0)
        {
            has_unshrinked = true;
            break;
        }
    }
    rewinddir(dir);
    const std::string prefix = has_unshrinked ? "coulomb_unshrinked_cut_" : "coulomb_cut_";
    std::cout << "read_coulomb_mat_k: using prefix \"" << prefix << "\" in directory \"" << in_dir << "\"" << std::endl;

    size_t nk = 0, nabf = 0, istart = 0, jstart = 0, iend = 0, jend = 0;

    K_Vectors* const klist = &(kRlist.klist_coarse);
    int klist_nk = klist->nmp[0] * klist->nmp[1] * klist->nmp[2];
    int ik_readin = -1;
    TLRI<TVs> Vs;
    std::map<int, std::map<std::pair<int,int>, RI::Tensor<std::complex<double>>>> Vq; // <iat1, <<iat2,ik>, T>>
    const int nat = Cs.size();
    std::vector<size_t> abf_start_index(nat+1, 1);
    for (int iat = 0; iat < nat; ++iat)
    {
        abf_start_index[iat+1] = abf_start_index[iat] + Cs.at(iat).at({ 0, {0,0,0} }).shape[0];
    }

    auto to_atom = [&](const int start, const int end) -> int
    {
        for (int iat = 0;iat < nat;++iat)
        {
            size_t abf_start = abf_start_index[iat];
            size_t abf_end = abf_start_index[iat+1] - 1;
            if (start == abf_start && end == abf_end)
            {
                return iat;
            }
        }
        throw std::runtime_error("Error in read_coulomb_mat_k: cannot find the atom for given auxiliary basis set range");
    };

    while ((ptr = readdir(dir)) != NULL){// read all the files in the directory
        std::string fm(ptr->d_name);
        if (fm.find(prefix) == 0)// find file coulomb_cut_xxx
        {
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "found Coulomb file: " + fm + ", start reading...");
            std::ifstream ifs(in_dir  + fm);
            ifs >> nk;//   actual nk
            assert(nk == klist_nk);

            while (ifs.peek() != EOF)
            {
                ifs >> nabf >> istart >> iend >> jstart >> jend >> ik_readin >> klist->wk[ik_readin-1];
                if (ifs.peek() == EOF) { break; }
                int ik = ik_readin - 1;
                int iat1 = to_atom(istart, iend);
                int iat2 = to_atom(jstart, jend);
                const size_t nabf1 = iend - istart + 1;
                const size_t nabf2 = jend - jstart + 1;             

                RI::Tensor<std::complex<double>> t({ nabf1, nabf2 });
                for (int i = 0;i < nabf1;++i)
                {
                    for (int j = 0;j < nabf2;++j)
                    {
                        LR_IO::read_one_data(ifs, t(i, j));
                    }
                }
                Vq[iat1][{iat2, ik}] = t;

                if (iat1 != iat2)
                {   // coulomb_mat has only the upper triangle part
                    Vq[iat2][{iat1, ik}] = t.dagger();
                    continue;
                }   
            }
            ifs.close();
        }
    }
    closedir(dir);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read Vq files. Now convert Vq to VR.");
    for (const TC& R : kRlist.Rlist )
    {
        for (int iat1 = 0;iat1 < nat;++iat1)
        {
            for (int iat2 = 0;iat2 < nat;++iat2)
            {
                Vs[iat1][{iat2, R}] = RI::Tensor<TVs>(Vq[iat1][{iat2, 0}].shape);
            }
        }
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "VR keys has been prepared.");
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) collapse(3)
#endif
    for (const TC& R : kRlist.Rlist )
    {
// #ifdef _OPENMP
// #pragma omp critical
//         std::cout<<"thread"<< omp_get_thread_num() << " convert V: R="<<R[0]<<" "<<R[1]<<" "<<R[2]<<std::endl;
// #endif
        for (int iat1 = 0;iat1 < nat;++iat1)
        {
            for (int iat2 = 0;iat2 < nat;++iat2)
            {
                Vs[iat1][{iat2, R}] = RI::Tensor<TVs>(Vq[iat1][{iat2, 0}].shape);
                for (int ik = 0;ik < nk;++ik)
                {
                    const ModuleBase::Vector3<double>& kvec = klist->kvec_d.at(ik);
                    const double arg = -1.0 * ModuleBase::TWO_PI * (kvec.x * R[0] + kvec.y * R[1] + kvec.z * R[2]);
                    const std::complex<double> kphase (cos(arg), sin(arg));
                    Vs[iat1][{iat2, R}] += RI::Global_Func::convert<TVs> (Vq[iat1][{iat2, ik}] * kphase) * RI::Global_Func::convert<TVs>(klist->wk[ik]);
                }
            }
        }
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "convert Vq to VR.");
    ModuleBase::TITLE("LR_IO", "read_Vs done.");
    return Vs;
}

template <typename TCs, typename TVs> // any blocking
TLRI<TVs> read_coulomb_mat_general_k(const std::string& in_dir, const TLRI<TCs>& Cs, LR_IO::RI_kRlist& kRlist)
{
    struct dirent *ptr;
    DIR *dir;
    bool has_unshrinked = false;
    dir = opendir(in_dir.c_str());
    if (!dir)
    {
        throw std::runtime_error("Cannot open directory: " + in_dir);
    }
    while ((ptr = readdir(dir)) != nullptr)
    {
        std::string fm(ptr->d_name);
        if (fm.find("coulomb_unshrinked_cut_") == 0)
        {
            has_unshrinked = true;
            break;
        }
    }
    rewinddir(dir);
    const std::string prefix = has_unshrinked ? "coulomb_unshrinked_cut_" : "coulomb_cut_";
    std::cout << "read_coulomb_mat_k: using prefix \"" << prefix << "\" in directory \"" << in_dir << "\"" << std::endl;
    
    TLRI<TVs> Vs;
    std::map<int, std::map<std::pair<int,int>, RI::Tensor<std::complex<double>>>> Vq; // <iat1, <<iat2,ik>, T>>
    std::map<int,std::vector<std::complex<double>>> Vq_tmp; //<ik, vector> 

    size_t nabf = 0, istart = 0, jstart = 0, iend = 0, jend = 0;
    int nk = 0;
    K_Vectors* const klist = &(kRlist.klist_coarse);
    int klist_nk = klist->nmp[0] * klist->nmp[1] * klist->nmp[2];
    int ik_readin = -1;

    while ((ptr = readdir(dir)) != NULL){// read all the files in the directory
        std::string fm(ptr->d_name);
        if (fm.find(prefix) == 0)// find file coulomb_cut_xxx
        {
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "found Coulomb file: " + fm + ", start reading...");
            std::ifstream ifs(in_dir  + fm);
            ifs >> nk;  //   actual nk            
            assert(nk == klist_nk);
            
            while (ifs.peek() != EOF)
            {
                ifs >> nabf >> istart >> iend >> jstart >> jend >> ik_readin >> klist->wk[ik_readin-1];// wk is not used
                if (ifs.peek() == EOF) { break; }
                int ik = ik_readin - 1;
                if (Vq_tmp[ik].empty()) { Vq_tmp[ik].resize(nabf * nabf, 0.0); }
                auto& Vq_tmp_k = Vq_tmp.at(ik);
                for (int i = istart - 1;i < iend;++i)
                {
                    for (int j = jstart - 1;j < jend;++j)
                    {
                        LR_IO::read_one_data(ifs, Vq_tmp_k[i * nabf + j]);
                    }
                }
            }
            ifs.close();
        }
    }
    closedir(dir);

    const int nat = Cs.size();
    istart = 0;
    for (int iat1 = 0;iat1 < nat;++iat1)
    {
        const size_t nabf1 = Cs.at(iat1).at({ 0, {0,0,0} }).shape[0];
        jstart = 0;
        for (int iat2 = 0;iat2 < nat;++iat2)
        {
            const size_t nabf2 = Cs.at(iat2).at({ 0, {0,0,0} }).shape[0];
            for (int ik = 0; ik < nk; ++ik){                    
                if (iat1 > iat2)
                {   // coulomb_mat has only the upper triangle part
                    Vq[iat1][{iat2, ik}] = Vq[iat2][{iat1, ik}].dagger();
                }
                else
                {
                    RI::Tensor<std::complex<double>> t({ nabf1, nabf2 });
                    for (int i = 0;i < nabf1;++i)
                    {
                        for (int j = 0;j < nabf2;++j)
                        {
                            t(i, j) = Vq_tmp[ik][(istart + i) * nabf + jstart + j];
                        }
                    }
                    Vq[iat1][{iat2, ik}] = t;
                }
            }
            jstart += nabf2;
        }
        assert(jstart == nabf);
        istart += nabf1;
    }
    assert(istart == nabf);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read Vq files. Now convert Vq to VR.");
    for (const TC& R : kRlist.Rlist )
    {
        for (int iat1 = 0;iat1 < nat;++iat1)
        {
            for (int iat2 = 0;iat2 < nat;++iat2)
            {
                Vs[iat1][{iat2, R}] = RI::Tensor<TVs>({ Vq[iat1][{iat2, 0}].shape[0], Vq[iat1][{iat2, 0}].shape[1] });
            }
        }
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "VR keys has been prepared.");

    double reciprocal_nk = 1.0 / double(nk);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) collapse(3)
#endif
    for (const TC& R : kRlist.Rlist )
    {
// #ifdef _OPENMP
// #pragma omp critical
//         std::cout<<"thread"<< omp_get_thread_num() << "convert V: R="<<R[0]<<" "<<R[1]<<" "<<R[2]<<std::endl;
// #endif
        for (int iat1 = 0;iat1 < nat;++iat1)
        {
            for (int iat2 = 0;iat2 < nat;++iat2)
            {
                for (int ik = 0; ik < nk; ++ik)
                {
                    const ModuleBase::Vector3<double>& kvec = klist->kvec_d.at(ik);
                    const double arg = -1.0 * ModuleBase::TWO_PI * (kvec.x * R[0] + kvec.y * R[1] + kvec.z * R[2]);
                    const std::complex<double> kphase (cos(arg), sin(arg));
                    Vs[iat1][{iat2, R}] += RI::Global_Func::convert<TVs> (Vq[iat1][{iat2, ik}] * kphase) * RI::Global_Func::convert<TVs>(reciprocal_nk);
                }
            }
        }
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "convert Vq to VR.");
    ModuleBase::TITLE("LR_IO", "read_Vs done.");
    return Vs;
}

template<typename Tdata, typename TVs>
TLRI<Tdata> read_Ws(const TLRI<TVs>& Vs, const std::vector<TC>& Rlist)
{
    ModuleBase::TITLE("LR_IO", "read_Ws");
    std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> Ws;
    
    const int nat = Vs.size();
    std::string temp;
    int nk, istart, iend, jstart, jend, ik;
    size_t nabfmu, nabfnu, non_zero, mu, nu; //I.nab, J.nab
    size_t nR = Rlist.size();
    for(int iat = 0; iat != nat; ++iat)//loop atom I
    {
        for(int jat = 0; jat != nat; ++jat)//loop atom J
        {
            for(int iR = 0; iR < nR; ++iR)
            {
                std::ifstream infileW;
                std::string filename = "librpa.d/Wc_Mu_"+std::to_string(iat)+"_Nu_"+std::to_string(jat)+"_iR_"+std::to_string(iR)+"_ifreq_0.mtx";
                infileW.open(filename);
                if(!infileW) throw std::runtime_error( filename + " not found!");
                // else std::cout << "reading Wc file: " << filename ;
                int nabf1 = Vs.at(iat).at({jat,{0,0,0}}).shape[0];
                int nabf2 = Vs.at(iat).at({jat,{0,0,0}}).shape[1];

                TC R; // iR of Wc file is not equal to iR in Rlist !!!
                infileW.ignore(2048, '\n'); // skip line 1: %%MatrixMarket...
                std::getline(infileW, temp); // read line 2: "%"
                std::getline(infileW, temp); // read line 3: "% Wc at iR N ( Rx Ry Rz ) ..."
                size_t lparen = temp.find('(');
                size_t rparen = temp.find(')', lparen);
                if (lparen == std::string::npos || rparen == std::string::npos)
                    throw std::runtime_error("Failed to parse R coordinates in " + filename);
                std::istringstream riss(temp.substr(lparen + 1, rparen - lparen - 1));
                riss >> R[0] >> R[1] >> R[2];
                while(infileW.peek() == '%') infileW.ignore(2048, '\n');	//skip comments

                infileW >> nabfmu >> nabfnu >> non_zero;
                assert(nabfmu == nabf1);
                assert(nabfnu == nabf2);
                RI::Tensor<Tdata> tensor_W({ nabfmu, nabfnu });
                for (int index = 0; index < non_zero; ++index)
                {
                    infileW >> mu >> nu ;
                    LR_IO::read_one_data(infileW, tensor_W(mu-1, nu-1));
                }
                infileW.close();
                tensor_W += Vs.at(iat).at({jat, R});
                // for(int i = 0; i != nabf1; ++i)
                //     for(int j = 0; j != nabf2; ++j)
                //     {
                //         tensor_W(i, j) += Vs.at(iat).at({jat, R})(i,j);
                //         std::cout << "Wxc: " << i << " " << j << " " << tensor_W(i,j) << std::endl; //check
                //     }
                Ws[iat][{jat, R}] = std::move(tensor_W);
                // std::cout << " Finished. R: " << "( " << R[0] << " " << R[1] << " " << R[2] << " )" << std::endl;
            }
        }
    }
    ModuleBase::TITLE("LR_IO", "read_Ws done.");
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "read WR files.");
    return Ws;
}
template<typename T>
void write_lri_R_max_norm(const TLRI<T>& tensors,
                            const UnitCell& ucell,
                            const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs) { throw std::runtime_error("Cannot open " + filename); }
    ofs << "# iat jat Rx Ry Rz Rnorm_bohr tensor_max_abs\n";
    ofs << std::setprecision(16);
    for (const auto& iat_blocks : tensors)
    {
        for (const auto& pair_tensor : iat_blocks.second)
        {
            const int jat = pair_tensor.first.first;
            const auto& R = pair_tensor.first.second;
            const ModuleBase::Vector3<double> R_cart =
                (static_cast<double>(R[0]) * ucell.a1
                + static_cast<double>(R[1]) * ucell.a2
                + static_cast<double>(R[2]) * ucell.a3) * ucell.lat0;
            ofs << iat_blocks.first << ' ' << jat << ' '
                << R[0] << ' ' << R[1] << ' ' << R[2] << ' ' << R_cart.norm() << ' '
                << pair_tensor.second.norm(std::numeric_limits<double>::max()) << '\n';
        }
    }
}

template TLRI<double> read_coulomb_mat_k<double, double>
    (const std::string& in_dir, const TLRI<double>& Cs, LR_IO::RI_kRlist& kRlist);
template TLRI<std::complex<double>> read_coulomb_mat_k<std::complex<double>, std::complex<double>>
    (const std::string& in_dir, const TLRI<std::complex<double>>& Cs, LR_IO::RI_kRlist& kRlist);

template TLRI<double> read_coulomb_mat_general_k<double, double>
    (const std::string& in_dir, const TLRI<double>& Cs, LR_IO::RI_kRlist& kRlist);
template TLRI<std::complex<double>> read_coulomb_mat_general_k<std::complex<double>, std::complex<double>>
    (const std::string& in_dir, const TLRI<std::complex<double>>& Cs, LR_IO::RI_kRlist& kRlist);

template TLRI<double> read_Ws<double, double>
    (const TLRI<double>& Vs, const std::vector<TC>& Rlist);
template TLRI<std::complex<double>> read_Ws<std::complex<double>, std::complex<double>>
    (const TLRI<std::complex<double>>& Vs, const std::vector<TC>& Rlist);

template void write_lri_R_max_norm<double>
    (const TLRI<double>& tensors, const UnitCell& ucell, const std::string& filename);
template void write_lri_R_max_norm<std::complex<double>>
    (const TLRI<std::complex<double>>& tensors, const UnitCell& ucell, const std::string& filename);

#endif // __EXX

} // namespace LR_IO
