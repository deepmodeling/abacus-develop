//=======================
// AUTHOR : Rong Shi
// DATE :   2022-12-09
//=======================

#ifndef RPA_LRI_HPP
#define RPA_LRI_HPP
#include "RPA_LRI.h"

#include <cstring>
#include <fstream>
#include <iostream>

template <typename Tdata> void RPA_LRI<Tdata>::init(const MPI_Comm &mpi_comm_in, const K_Vectors &kv_in)
{
    ModuleBase::TITLE("RPA_LRI", "init");
    ModuleBase::timer::tick("RPA_LRI", "init");
    this->mpi_comm = mpi_comm_in;
    this->lcaos = exx_lri_rpa.lcaos;
    this->abfs = exx_lri_rpa.abfs;
    this->abfs_ccp = exx_lri_rpa.abfs_ccp;
    this->p_kv = &kv_in;

    //	this->cv = std::move(exx_lri_rpa.cv);
    //    exx_lri_rpa.cv = exx_lri_rpa.cv;
}

template <typename Tdata> void RPA_LRI<Tdata>::cal_rpa_cv()
{
    std::vector<TA> atoms(GlobalC::ucell.nat);
    for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};

    const std::array<Tcell, Ndim> period_Vs = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs
        = exx_lri_rpa.cv.cal_Vs(list_As_Vs.first, list_As_Vs.second[0], {{"writable_Vws", true}});
    this->Vs_period = RI::RI_Tools::cal_period(Vs, period);

    const std::array<Tcell, Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Cs
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);

    std::pair<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>,
              std::array<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>, 3>>
        Cs_dCs = exx_lri_rpa.cv.cal_Cs_dCs(list_As_Cs.first,
                                           list_As_Cs.second[0],
                                           {{"cal_dC", false},
                                            {"writable_Cws", true},
                                            {"writable_dCws", true},
                                            {"writable_Vws", false},
                                            {"writable_dVws", false}});
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> &Cs = std::get<0>(Cs_dCs);
    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
}

template <typename Tdata>
void RPA_LRI<Tdata>::cal_postSCF_exx(const Local_Orbital_Charge& loc,
                const MPI_Comm& mpi_comm_in,
                const K_Vectors& kv,
                const Parallel_Orbitals& pv)
{
    exx_lri_rpa.mix_DMk_2D.set_nks(kv.nks, GlobalV::GAMMA_ONLY_LOCAL);
    exx_lri_rpa.mix_DMk_2D.set_mixing_mode(Mixing_Mode::No);
    if(GlobalV::GAMMA_ONLY_LOCAL)
        exx_lri_rpa.mix_DMk_2D.mix(loc.dm_gamma, true);
    else
        exx_lri_rpa.mix_DMk_2D.mix(loc.dm_k, true);
    exx_lri_rpa.init(mpi_comm_in, kv);
    exx_lri_rpa.cal_exx_ions();
    exx_lri_rpa.cal_exx_elec(pv);
    // cout<<"postSCF_Eexx: "<<exx_lri_rpa.Eexx<<endl;
}

template <typename Tdata>
void RPA_LRI<Tdata>::out_for_RPA(const Parallel_Orbitals &parav,
                                 const psi::Psi<std::complex<double>> &psi,
                                 const elecstate::ElecState *pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_for_RPA");
    this->out_bands(pelec);
    this->out_eigen_vector(parav, psi);
    this->out_struc();

    this->cal_rpa_cv();
    std::cout << "rpa_pca_threshold: " << this->info.pca_threshold << std::endl;
    std::cout << "rpa_ccp_rmesh_times: " << this->info.ccp_rmesh_times << std::endl;
    std::cout << "rpa_lcao_exx(Ha): " << std::fixed << std::setprecision(15) << exx_lri_rpa.Eexx / 2.0 << std::endl;
    this->out_Cs();
    this->out_coulomb_k();

    std::cout << "etxc(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etxc / 2.0 << std::endl;
    std::cout << "etot(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etot / 2.0 << std::endl;
    std::cout << "Etot_without_rpa(Ha): " << std::fixed << std::setprecision(15)
              << (pelec->f_en.etot - pelec->f_en.etxc + exx_lri_rpa.Eexx) / 2.0 << std::endl;

    return;
}

template <typename Tdata>
void RPA_LRI<Tdata>::out_eigen_vector(const Parallel_Orbitals &parav, const psi::Psi<std::complex<double>> &psi)
{

    ModuleBase::TITLE("DFT_RPA_interface", "out_eigen_vector");

    const int nks_tot = GlobalV::NSPIN == 2 ? p_kv->nks / 2 : p_kv->nks;
    const int npsin_tmp = GlobalV::NSPIN == 2 ? 2 : 1;
    const std::complex<double> zero(0.0, 0.0);

    for (int ik = 0; ik < nks_tot; ik++)
    {
        std::stringstream ss;
        ss << "KS_eigenvector_" << ik << ".dat";

        std::ofstream ofs;
        if (GlobalV::MY_RANK == 0)
            ofs.open(ss.str().c_str(), std::ios::out);
        std::vector<ModuleBase::ComplexMatrix> is_wfc_ib_iw(npsin_tmp);
        for (int is = 0; is < npsin_tmp; is++)
        {
            is_wfc_ib_iw[is].create(GlobalV::NBANDS, GlobalV::NLOCAL);
            for (int ib_global = 0; ib_global < GlobalV::NBANDS; ++ib_global)
            {
                std::vector<std::complex<double>> wfc_iks(GlobalV::NLOCAL, zero);

                const int ib_local = parav.global2local_col[ib_global];

                if (ib_local >= 0)
                    for (int ir = 0; ir < psi.get_nbasis(); ir++)
                        wfc_iks[parav.local2global_row[ir]] = psi(ik + nks_tot * is, ib_local, ir);

                std::vector<std::complex<double>> tmp = wfc_iks;
#ifdef __MPI
                MPI_Allreduce(&tmp[0], &wfc_iks[0], GlobalV::NLOCAL, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
                for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
                    is_wfc_ib_iw[is](ib_global, iw) = wfc_iks[iw];
            } // ib
        } // is
        ofs << ik + 1 << std::endl;
        for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                for (int is = 0; is < npsin_tmp; is++)
                    ofs << std::setw(21) << std::fixed << std::setprecision(15) << is_wfc_ib_iw[is](ib, iw).real()
                        << std::setw(21) << std::fixed << std::setprecision(15) << is_wfc_ib_iw[is](ib, iw).imag()
                        << std::endl;
            }
        }
        ofs.close();
    } // ik
    return;
}

template <typename Tdata> void RPA_LRI<Tdata>::out_struc()
{
    if (GlobalV::MY_RANK != 0)
        return;
    ModuleBase::TITLE("DFT_RPA_interface", "out_struc");
    double TWOPI_Bohr2A = ModuleBase::TWO_PI * ModuleBase::BOHR_TO_A;
    const int nks_tot = GlobalV::NSPIN == 2 ? (int)p_kv->nks / 2 : p_kv->nks;
    ModuleBase::Matrix3 lat = GlobalC::ucell.latvec / ModuleBase::BOHR_TO_A;
    ModuleBase::Matrix3 G = GlobalC::ucell.G * TWOPI_Bohr2A;
    std::stringstream ss;
    ss << "stru_out";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << lat.e11 << std::setw(15) << lat.e12 << std::setw(15) << lat.e13 << std::endl;
    ofs << lat.e21 << std::setw(15) << lat.e22 << std::setw(15) << lat.e23 << std::endl;
    ofs << lat.e31 << std::setw(15) << lat.e32 << std::setw(15) << lat.e33 << std::endl;

    ofs << G.e11 << std::setw(15) << G.e12 << std::setw(15) << G.e13 << std::endl;
    ofs << G.e21 << std::setw(15) << G.e22 << std::setw(15) << G.e23 << std::endl;
    ofs << G.e31 << std::setw(15) << G.e32 << std::setw(15) << G.e33 << std::endl;

    ofs << p_kv->nmp[0] << std::setw(6) << p_kv->nmp[1] << std::setw(6) << p_kv->nmp[2]
        << std::setw(6) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
        ofs << std::setw(15) << std::fixed << std::setprecision(9) << p_kv->kvec_c[ik].x * TWOPI_Bohr2A
            << std::setw(15) << std::fixed << std::setprecision(9) << p_kv->kvec_c[ik].y * TWOPI_Bohr2A
            << std::setw(15) << std::fixed << std::setprecision(9) << p_kv->kvec_c[ik].z * TWOPI_Bohr2A
            << std::endl;
    ofs.close();
    return;
}

template <typename Tdata> void RPA_LRI<Tdata>::out_bands(const elecstate::ElecState *pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_bands");
    if (GlobalV::MY_RANK != 0)
        return;
    const int nks_tot = GlobalV::NSPIN == 2 ? (int)p_kv->nks / 2 : p_kv->nks;
    const int nspin_tmp = GlobalV::NSPIN == 2 ? 2 : 1;
    std::stringstream ss;
    ss << "band_out";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << nks_tot << std::endl;
    ofs << GlobalV::NSPIN << std::endl;
    ofs << GlobalV::NBANDS << std::endl;
    ofs << GlobalV::NLOCAL << std::endl;
    ofs << (pelec->eferm.ef / 2.0) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (int is = 0; is != nspin_tmp; is++)
        {
            ofs << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << std::endl;
            for (int ib = 0; ib != GlobalV::NBANDS; ib++)
                ofs << std::setw(5) << ib + 1 << "   " << std::setw(8) << pelec->wg(ik + is * nks_tot, ib) * nks_tot
                    << std::setw(18) << std::fixed << std::setprecision(8) << pelec->ekb(ik + is * nks_tot, ib) / 2.0
                    << std::setw(18) << std::fixed << std::setprecision(8)
                    << pelec->ekb(ik + is * nks_tot, ib) * ModuleBase::Ry_to_eV << std::endl;
        }
    }
    ofs.close();
    return;
}

template <typename Tdata> void RPA_LRI<Tdata>::out_Cs()
{
    std::stringstream ss;
    ss << "Cs_data_" << GlobalV::MY_RANK << ".txt";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << GlobalC::ucell.nat << "    " << 0 << std::endl;
    for (auto &Ip: this->Cs_period)
    {
        size_t I = Ip.first;
        size_t i_num = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[I]].nw;
        for (auto &JPp: Ip.second)
        {
            size_t J = JPp.first.first;
            auto R = JPp.first.second;
            auto &tmp_Cs = JPp.second;
            size_t j_num = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[J]].nw;

            ofs << I + 1 << "   " << J + 1 << "   " << R[0] << "   " << R[1] << "   " << R[2] << "   " << i_num
                << std::endl;
            ofs << j_num << "   " << tmp_Cs.shape[0] << std::endl;
            for (int i = 0; i != i_num; i++)
                for (int j = 0; j != j_num; j++)
                    for (int mu = 0; mu != tmp_Cs.shape[0]; mu++)
                        ofs << std::setw(15) << std::fixed << std::setprecision(9) << tmp_Cs(mu, i, j) << std::endl;
        }
    }
    ofs.close();
    return;
}

template <typename Tdata> void RPA_LRI<Tdata>::out_coulomb_k()
{
    int all_mu = 0;
    vector<int> mu_shift(GlobalC::ucell.nat);
    for (int I = 0; I != GlobalC::ucell.nat; I++)
    {
        mu_shift[I] = all_mu;
        all_mu += exx_lri_rpa.cv.get_index_abfs_size(GlobalC::ucell.iat2it[I]);
    }
    const int nks_tot = GlobalV::NSPIN == 2 ? (int)p_kv->nks / 2 : p_kv->nks;
    std::stringstream ss;
    ss << "coulomb_mat_" << GlobalV::MY_RANK << ".txt";

    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);

    ofs << nks_tot << std::endl;
    for (auto &Ip: this->Vs_period)
    {
        auto I = Ip.first;
        size_t mu_num = exx_lri_rpa.cv.get_index_abfs_size(GlobalC::ucell.iat2it[I]);

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> Vq_k_IJ;
            for (auto &JPp: Ip.second)
            {
                auto J = JPp.first.first;

                auto R = JPp.first.second;
                if (J < I)
                    continue;
                RI::Tensor<std::complex<double>> tmp_VR = RI::Global_Func::convert<std::complex<double>>(JPp.second);

                const double arg = 1
                                   * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * GlobalC::ucell.latvec))
                                   * ModuleBase::TWO_PI; // latvec
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
                if (Vq_k_IJ[J].empty())
                    Vq_k_IJ[J] = RI::Tensor<std::complex<double>>({tmp_VR.shape[0], tmp_VR.shape[1]});
                Vq_k_IJ[J] = Vq_k_IJ[J] + tmp_VR * kphase;
            }
            for (auto &vq_Jp: Vq_k_IJ)
            {
                auto iJ = vq_Jp.first;
                auto &vq_J = vq_Jp.second;
                size_t nu_num = exx_lri_rpa.cv.get_index_abfs_size(GlobalC::ucell.iat2it[iJ]);
                ofs << all_mu << "   " << mu_shift[I] + 1 << "   " << mu_shift[I] + mu_num << "  " << mu_shift[iJ] + 1
                    << "   " << mu_shift[iJ] + nu_num << std::endl;
                ofs << ik + 1 << "  " << p_kv->wk[ik] / 2.0 * GlobalV::NSPIN << std::endl;
                for (int i = 0; i != vq_J.data->size(); i++)
                {
                    ofs << std::setw(21) << std::fixed << std::setprecision(12) << (*vq_J.data)[i].real()
                        << std::setw(21) << std::fixed << std::setprecision(12) << (*vq_J.data)[i].imag() << std::endl;
                }
            }
        }
    }
    ofs.close();
}

// template<typename Tdata>
// void RPA_LRI<Tdata>::init(const MPI_Comm &mpi_comm_in)
// {
// 	if(this->info == this->exx.info)
// 	{
// 		this->lcaos = this->exx.lcaos;
// 		this->abfs = this->exx.abfs;
// 		this->abfs_ccp = this->exx.abfs_ccp;

// 		exx_lri_rpa.cv = std::move(this->exx.cv);
// 	}
// 	else
// 	{
// 		this->lcaos = ...
// 		this->abfs = ...
// 		this->abfs_ccp = ...

// 		exx_lri_rpa.cv.set_orbitals(
// 			this->lcaos, this->abfs, this->abfs_ccp,
// 			this->info.kmesh_times, this->info.ccp_rmesh_times );
// 	}

// //	for( size_t T=0; T!=this->abfs.size(); ++T )
// //		GlobalC::exx_info.info_ri.abfs_Lmax = std::max( GlobalC::exx_info.info_ri.abfs_Lmax,
// static_cast<int>(this->abfs[T].size())-1 );

// }

// template<typename Tdata>
// void RPA_LRI<Tdata>::cal_rpa_ions()
// {
// 	// this->rpa_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

// 	if(this->info == this->exx.info)
// 		exx_lri_rpa.cv.Vws = std::move(this->exx.cv.Vws);

// 	const std::array<Tcell,Ndim> period_Vs = LRI_CV_Tools::cal_latvec_range<Tcell>(1+this->info.ccp_rmesh_times);
// 	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
// 		list_As_Vs = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);

// 	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>
// 		Vs = exx_lri_rpa.cv.cal_Vs(
// 			list_As_Vs.first, list_As_Vs.second[0],
// 			{{"writable_Vws",true}});

// 	// Vs[iat0][{iat1,cell1}]	按 (iat0,iat1) 分进程，每个进程有所有 cell1
// 	Vqs = FFT(Vs);
// 	out_Vs(Vqs);

// 	if(this->info == this->exx.info)
// 		exx_lri_rpa.cv.Cws = std::move(this->exx.cv.Cws);

// 	const std::array<Tcell,Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2);
// 	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
// 		list_As_Cs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);

// 	std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>, std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>>
// 		Cs_dCs = exx_lri_rpa.cv.cal_Cs_dCs(
// 			list_As_Cs.first, list_As_Cs.second[0],
// 			{{"cal_dC",false},
// 			 {"writable_Cws",true}, {"writable_dCws",true}, {"writable_Vws",false}, {"writable_dVws",false}});
// 	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> &Cs = std::get<0>(Cs_dCs);

// 	out_Cs(Cs);

// 	// rpa_lri.set_Cs(Cs);
// }

#endif