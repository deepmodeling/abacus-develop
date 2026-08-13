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
    const std::string FILE_COARSE = "stru_out";
    const std::string FILE_FINE_UNIFORM = "band_kpath_info";
    const std::string FILE_FINE_NONUNIFORM = "KPT_bse";
    const std::string FILE_BAND_OUT = "band_out";
    const std::string FILE_BAND_KPATH = "band_kpath_info";

void parse_band_out_file(const std::string& in_dir, int& nbands_file, int& nk_file, int& nspin_file, int& nocc_file)
{
    std::string file = in_dir + FILE_BAND_OUT;
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
        file = in_dir + FILE_BAND_KPATH;
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

#ifdef __EXX

RI_kRlist::RI_kRlist(const UnitCell& ucell, K_Vectors* const pkv,
    const std::string& in_dir, const int use_fine_kgrid, const std::string& out_dir)
    : klist(pkv)
{
    read_kpts_coarse(in_dir + FILE_COARSE, ucell, this->klist, out_dir);
    this->klist_coarse = *this->klist;
    this->period = RI_Util::get_Born_vonKarmen_period(*klist);
    this->Rlist = RI_Util::get_Born_von_Karmen_cells(period);
    // std::cout << "Rlist:" << std::endl;
    // int count = 0;
    // for (const auto& iR: Rlist)
    // {
    //     count++;
    //     std::cout << "iR=" << count <<": "<< iR[0] << " " << iR[1] << " " << iR[2] << std::endl;
    // }
    if (use_fine_kgrid==1)
    {
        read_kpts_fine(in_dir + FILE_FINE_UNIFORM, ucell, this->klist, false, out_dir);
    }
    else if (use_fine_kgrid==2)
    {
        read_kpts_fine(in_dir + FILE_FINE_NONUNIFORM, ucell, this->klist, true, out_dir);
    }
    else if (use_fine_kgrid!=0)
        ModuleBase::WARNING_QUIT("LR_IO", "use_fine_kgrid must be 0, 1 or 2");
};

void RI_kRlist::read_kpts_coarse(const std::string& file, const UnitCell& ucell,
                                 K_Vectors* const klist, const std::string& out_dir)
{
    std::ifstream ifs;
    ifs.open(file);
    if (!ifs) throw std::runtime_error(file + " not found");
    std::string tmp;
    for (int i = 0; i < 7; ++i) { std::getline(ifs, tmp); } // get the 7th line(number of atoms)
    int nat = std::stoi(tmp);
    for (int i = 0; i != nat; ++i) { std::getline(ifs, tmp); }
    int nks_original = klist->get_nks();
    // std::cout << "Origianl klist (Cartesian|Direct)" << std::endl;
    // for (int ik = 0;ik < nks_original;++ik)
    // {
    //     std::cout << "ik=" << std::setw(5) << ik << std::setw(11) << klist->kvec_c[ik].x << std::setw(11) 
    //     << klist->kvec_c[ik].y << std::setw(11) << klist->kvec_c[ik].z << " | " << std::setw(11)
    //     << klist->kvec_d[ik].x << std::setw(11) << klist->kvec_d[ik].y << std::setw(11) << klist->kvec_d[ik].z << std::endl;
    // }

    ifs >> klist->nmp[0] >> klist->nmp[1] >> klist->nmp[2];
    int nk = klist->nmp[0] * klist->nmp[1] * klist->nmp[2];
    int nks = (PARAM.inp.nspin == 2) ? 2 * nk : nk;
    assert(nks == nks_original);

    for (int ik = 0; ik < nk; ++ik)
    {
        ifs >> klist->kvec_c[ik].x >> klist->kvec_c[ik].y >> klist->kvec_c[ik].z;
        klist->kvec_c[ik] /= ModuleBase::TWO_PI * ModuleBase::BOHR_TO_A; // in unit of 2pi/angstrom
        klist->kvec_d[ik] = klist->kvec_c[ik] * ucell.latvec.Transpose();
        set_zero_if_close(klist->kvec_d[ik]);
        klist->wk[ik] = 1.0 / double(nk);
    }
    if (PARAM.inp.nspin == 2)
    {
        for (int ik = 0; ik < nk; ++ik)
        {
            klist->kvec_c[ik + nk] = klist->kvec_c[ik];
            klist->kvec_d[ik + nk] = klist->kvec_d[ik];
            klist->wk[ik + nk] = klist->wk[ik];
        }
    }

    std::ofstream ofs_kpts_coarse(out_dir + "kpts_coarse.dat");
    ofs_kpts_coarse << "kpts_coarse:" << nk << std::setw(16) << "( Cartesian" << std::setw(36) 
        << "|                Direct )" << std::setw(15) << "| wk (normalized as sum = nk)" << std::endl;
    for (int ik = 0; ik < nks; ++ik)
    {
        ofs_kpts_coarse << std::setw(5) << ik << std::setw(12) << klist->kvec_c[ik].x << std::setw(12) 
        << klist->kvec_c[ik].y << std::setw(12) << klist->kvec_c[ik].z << " | " << std::setw(12)
        << klist->kvec_d[ik].x << std::setw(12) << klist->kvec_d[ik].y << std::setw(12) << klist->kvec_d[ik].z 
        << " | " << klist->wk[ik]*nk << std::endl;
    }
    ofs_kpts_coarse.close();
}

void RI_kRlist::read_kpts_fine(const std::string& file, const UnitCell& ucell,
                               K_Vectors* const klist, const bool is_weighted,
                               const std::string& out_dir)
{
    // band_kpath_info format: first line: nband nbasis nspin nk, then kx ky kz per line (direct coords)
    // KPT_bse format: first line = nk, then kx ky kz wk per line (direct coords, BSE weight sum=nk)
    std::ifstream ifs;
    ifs.open(file);
    if (!ifs) throw std::runtime_error(file + " not found");
    int nk;
    if (is_weighted) {ifs >> nk; ifs.ignore(2048, '\n');}
    else {ifs >> nk >> nk >> nk >> nk;}

    int nks = (PARAM.inp.nspin == 2) ? 2 * nk : nk;
    klist->set_nks(nks);
    klist->set_nkstot(nks);
    klist->set_nkstot_full(nk);

    auto klist_reset = [&klist](int kpoint_number){
        klist->kvec_c.resize(0);    klist->kvec_c.resize(kpoint_number);
        klist->kvec_d.resize(0);    klist->kvec_d.resize(kpoint_number);
        klist->wk.resize(0);        klist->wk.resize(kpoint_number);
        klist->isk.resize(0);
        klist->ngk.resize(0);
    };
    klist_reset(nks);

    for (int ik = 0; ik < nk; ++ik)
    {
        ifs >> klist->kvec_d[ik].x >> klist->kvec_d[ik].y >> klist->kvec_d[ik].z;
        if (is_weighted) {
            ifs >> klist->wk[ik];
            klist->wk[ik] /= double(nk);
        }
        else {klist->wk[ik] = 1.0 / double(nk);}
        klist->kvec_c[ik] = klist->kvec_d[ik] * ucell.G;
        set_zero_if_close(klist->kvec_c[ik]);
    }
    std::cout << "Read " << nk << " k-points and weights from " << file << std::endl;
    if (PARAM.inp.nspin == 2)
    {
        for (int ik = 0; ik < nk; ++ik)
        {
            klist->kvec_c[ik + nk] = klist->kvec_c[ik];
            klist->kvec_d[ik + nk] = klist->kvec_d[ik];
            klist->wk[ik + nk] = klist->wk[ik];
        }
    }
    std::ofstream ofs_kpts_fine(out_dir + "kpts_fine.dat");
    ofs_kpts_fine << "kpts_fine:" << nk << std::setw(18) << "( Cartesian" << std::setw(36) 
        << "|                Direct )" << std::setw(15) << "| wk (normalized as sum = nk)" << std::endl;
    for (int ik = 0; ik < nk; ++ik)
    {
        ofs_kpts_fine << std::setw(5) << ik << std::setw(12) << klist->kvec_c[ik].x << std::setw(12) 
        << klist->kvec_c[ik].y << std::setw(12) << klist->kvec_c[ik].z << " | " << std::setw(12)
        << klist->kvec_d[ik].x << std::setw(12) << klist->kvec_d[ik].y << std::setw(12) << klist->kvec_d[ik].z
        << " | " << klist->wk[ik]*nk << std::endl;
    }
    ofs_kpts_fine.close();
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

#endif // __EXX

} // namespace LR_IO
