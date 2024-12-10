//=======================
// AUTHOR : Peize Lin
#include "module_parameter/parameter.h"
// DATE :   2022-08-17
//=======================

#ifndef EXX_LRI_HPP
#define EXX_LRI_HPP

#include "Exx_LRI.h"
#include "RI_2D_Comm.h"
#include "RI_Util.h"
#include "module_ri/exx_abfs-construct_orbs.h"
#include "module_ri/exx_abfs-io.h"
#include "module_ri/conv_coulomb_pot_k.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"
#include "module_ri/serialization_cereal.h"
#include "module_ri/Mix_DMk_2D.h"
#include "module_basis/module_ao/parallel_orbitals.h"

#include <RI/distribute/Distribute_Equally.h>
#include <RI/global/Map_Operator-3.h>

#include <fstream>
#include <string>

template<typename Tdata>
void Exx_LRI<Tdata>::init(const MPI_Comm &mpi_comm_in, const K_Vectors &kv_in, const LCAO_Orbitals& orb)
{
	ModuleBase::TITLE("Exx_LRI","init");
	ModuleBase::timer::tick("Exx_LRI", "init");

	// if(GlobalC::exx_info.info_global.separate_loop)
	// {
	// 	Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::No;
	// 	Hexx_para.mixing_beta = 0;
	// }
	// else
	// {
	// 	if("plain"==GlobalC::CHR.mixing_mode)
	// 		Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Plain;
	// 	else if("pulay"==GlobalC::CHR.mixing_mode)
	// 		Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Pulay;
	// 	else
	// 		throw std::invalid_argument("exx mixing error. exx_separate_loop==false, mixing_mode!=plain or pulay");
	// 	Hexx_para.mixing_beta = GlobalC::CHR.mixing_beta;
	// }

	this->mpi_comm = mpi_comm_in;
	this->p_kv = &kv_in;
	this->orb_cutoff_ = orb.cutoffs();

	this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs( orb, this->info.kmesh_times );

	// #ifdef __MPI
	// Exx_Abfs::Util::bcast( this->info.files_abfs, 0, this->mpi_comm );
	// #endif

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs_same_atom = Exx_Abfs::Construct_Orbs::abfs_same_atom( orb, this->lcaos, this->info.kmesh_times, this->info.pca_threshold );
	if(this->info.files_abfs.empty())
		{ this->abfs = abfs_same_atom;}
	else
		{ this->abfs = Exx_Abfs::IO::construct_abfs( abfs_same_atom, orb, this->info.files_abfs, this->info.kmesh_times ); }
	Exx_Abfs::Construct_Orbs::print_orbs_size(this->abfs, GlobalV::ofs_running);

	auto get_ccp_parameter = [this]() -> std::map<std::string,double>
	{
		switch(this->info.ccp_type)
		{
			case Conv_Coulomb_Pot_K::Ccp_Type::Ccp:
				return {};
			case Conv_Coulomb_Pot_K::Ccp_Type::Hf:
			{
				// 4/3 * pi * Rcut^3 = V_{supercell} = V_{unitcell} * Nk
				const int nspin0 = (PARAM.inp.nspin==2) ? 2 : 1;
				const double hf_Rcut = std::pow(0.75 * this->p_kv->get_nkstot_full()/nspin0 * GlobalC::ucell.omega / (ModuleBase::PI), 1.0/3.0);
				return {{"hf_Rcut", hf_Rcut}};
			}
			case Conv_Coulomb_Pot_K::Ccp_Type::Hse:
				return {{"hse_omega", this->info.hse_omega}};
			default:
				throw std::domain_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));	break;
		}
	};
	this->abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(this->abfs, this->info.ccp_type, get_ccp_parameter(), this->info.ccp_rmesh_times);


	for( size_t T=0; T!=this->abfs.size(); ++T )
		{ GlobalC::exx_info.info_ri.abfs_Lmax = std::max( GlobalC::exx_info.info_ri.abfs_Lmax, static_cast<int>(this->abfs[T].size())-1 ); }

	this->cv.set_orbitals(
		orb,
		this->lcaos, this->abfs, this->abfs_ccp,
		this->info.kmesh_times, this->info.ccp_rmesh_times );

	ModuleBase::timer::tick("Exx_LRI", "init");
}

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_ions(const bool write_cv)
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_ions");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_ions");

	// init_radial_table_ions( cal_atom_centres_core(atom_pairs_core_origin), atom_pairs_core_origin );

	// this->m_abfsabfs.init_radial_table(Rradial);
	// this->m_abfslcaos_lcaos.init_radial_table(Rradial);

	std::vector<TA> atoms(GlobalC::ucell.nat);
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		{ atoms[iat] = iat; }
	std::map<TA,TatomR> atoms_pos;
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		{ atoms_pos[iat] = RI_Util::Vector3_to_array3( GlobalC::ucell.atoms[ GlobalC::ucell.iat2it[iat] ].tau[ GlobalC::ucell.iat2ia[iat] ] ); }
	const std::array<TatomR,Ndim> latvec
		= {RI_Util::Vector3_to_array3(GlobalC::ucell.a1),
		   RI_Util::Vector3_to_array3(GlobalC::ucell.a2),
		   RI_Util::Vector3_to_array3(GlobalC::ucell.a3)};
	const std::array<Tcell,Ndim> period = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};

	this->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

	// std::max(3) for gamma_only, list_A2 should contain cell {-1,0,1}. In the future distribute will be neighbour.
	const std::array<Tcell,Ndim> period_Vs = LRI_CV_Tools::cal_latvec_range<Tcell>(1+this->info.ccp_rmesh_times, orb_cutoff_);	
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As_Vs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);

	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>
		Vs = this->cv.cal_Vs(
			list_As_Vs.first, list_As_Vs.second[0],
			{{"writable_Vws",true}});
	this->cv.Vws = LRI_CV_Tools::get_CVws(Vs);
	if (write_cv && GlobalV::MY_RANK == 0)
		{ LRI_CV_Tools::write_Vs_abf(Vs, PARAM.globalv.global_out_dir + "Vs"); }
	this->exx_lri.set_Vs(std::move(Vs), this->info.V_threshold);

	if(PARAM.inp.cal_force || PARAM.inp.cal_stress)
	{
		std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>
			dVs = this->cv.cal_dVs(
				list_As_Vs.first, list_As_Vs.second[0],
				{{"writable_dVws",true}});
		this->cv.dVws = LRI_CV_Tools::get_dCVws(dVs);
		this->exx_lri.set_dVs(std::move(dVs), this->info.V_grad_threshold);
		if(PARAM.inp.cal_stress)
		{
			std::array<std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>,3> dVRs = LRI_CV_Tools::cal_dMRs(dVs);
			this->exx_lri.set_dVRs(std::move(dVRs), this->info.V_grad_R_threshold);
		}
	}

	const std::array<Tcell,Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As_Cs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);

	std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>, std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>>
		Cs_dCs = this->cv.cal_Cs_dCs(
			list_As_Cs.first, list_As_Cs.second[0],
			{{"cal_dC",PARAM.inp.cal_force||PARAM.inp.cal_stress},
			 {"writable_Cws",true}, {"writable_dCws",true}, {"writable_Vws",false}, {"writable_dVws",false}});
	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> &Cs = std::get<0>(Cs_dCs);
	this->cv.Cws = LRI_CV_Tools::get_CVws(Cs);
	if (write_cv && GlobalV::MY_RANK == 0)
		{ LRI_CV_Tools::write_Cs_ao(Cs, PARAM.globalv.global_out_dir + "Cs"); }
	this->exx_lri.set_Cs(std::move(Cs), this->info.C_threshold);

	if(PARAM.inp.cal_force || PARAM.inp.cal_stress)
	{
		std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3> &dCs = std::get<1>(Cs_dCs);
		this->cv.dCws = LRI_CV_Tools::get_dCVws(dCs);
		this->exx_lri.set_dCs(std::move(dCs), this->info.C_grad_threshold);
		if(PARAM.inp.cal_stress)
		{
			std::array<std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>,3> dCRs = LRI_CV_Tools::cal_dMRs(dCs);
			this->exx_lri.set_dCRs(std::move(dCRs), this->info.C_grad_R_threshold);
		}
	}
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_ions");
}

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_elec(const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Ds,
	const Parallel_Orbitals& pv,
	const ModuleSymmetry::Symmetry_rotation* p_symrot)
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_elec");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_elec");

	const std::vector<std::tuple<std::set<TA>, std::set<TA>>> judge = RI_2D_Comm::get_2D_judge(pv);

	if(p_symrot)
		{ this->exx_lri.set_symmetry(true, p_symrot->get_irreducible_sector()); }
	else
		{ this->exx_lri.set_symmetry(false, {}); }

	this->Hexxs.resize(PARAM.inp.nspin);
	this->Eexx = 0;
	for(int is=0; is<PARAM.inp.nspin; ++is)
	{
		const std::string suffix = ((PARAM.inp.cal_force || PARAM.inp.cal_stress) ? std::to_string(is) : "");

		this->exx_lri.set_Ds(Ds[is], this->info.dm_threshold, suffix);
		this->exx_lri.cal_Hs({ "","",suffix });

		if (!p_symrot)
		{
			this->Hexxs[is] = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
				this->mpi_comm, std::move(this->exx_lri.Hs), std::get<0>(judge[is]), std::get<1>(judge[is]));
		}
		else
		{
			// reduce but not repeat
			auto Hs_a2D = this->exx_lri.post_2D.set_tensors_map2(this->exx_lri.Hs);
			// rotate locally without repeat
			Hs_a2D = p_symrot->restore_HR(GlobalC::ucell.symm, GlobalC::ucell.atoms, GlobalC::ucell.st, 'H', Hs_a2D);
			// cal energy using full Hs without repeat
			this->exx_lri.energy = this->exx_lri.post_2D.cal_energy(
				this->exx_lri.post_2D.saves["Ds_" + suffix],
				this->exx_lri.post_2D.set_tensors_map2(Hs_a2D));
			// get repeated full Hs for abacus
			this->Hexxs[is] = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
				this->mpi_comm, std::move(Hs_a2D), std::get<0>(judge[is]), std::get<1>(judge[is]));
		}
		this->Eexx += std::real(this->exx_lri.energy);
		post_process_Hexx(this->Hexxs[is]);
	}
	this->Eexx = post_process_Eexx(this->Eexx);
	this->exx_lri.set_symmetry(false, {});
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_elec");	
}

template<typename Tdata>
void Exx_LRI<Tdata>::post_process_Hexx( std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> &Hexxs_io ) const
{
	ModuleBase::TITLE("Exx_LRI","post_process_Hexx");
	constexpr Tdata frac = -1 * 2;								// why?	Hartree to Ry?
	const std::function<void(RI::Tensor<Tdata>&)>
		multiply_frac = [&frac](RI::Tensor<Tdata> &t)
		{ t = t*frac; };
	RI::Map_Operator::for_each( Hexxs_io, multiply_frac );
}

template<typename Tdata>
double Exx_LRI<Tdata>::post_process_Eexx(const double& Eexx_in) const
{
	ModuleBase::TITLE("Exx_LRI","post_process_Eexx");
	const double SPIN_multiple = std::map<int, double>{ {1,2}, {2,1}, {4,1} }.at(PARAM.inp.nspin);				// why?
	const double frac = -SPIN_multiple;
	return frac * Eexx_in;
}

/*
post_process_old
{
	// D
	const std::map<int,double> SPIN_multiple = {{1,0.5}, {2,1}, {4,1}};							// ???
	DR *= SPIN_multiple.at(NSPIN);

	// H
	HR *= -2;

	// E
	const std::map<int,double> SPIN_multiple = {{1,2}, {2,1}, {4,1}};							// ???
	energy *= SPIN_multiple.at(PARAM.inp.nspin);			// ?
	energy /= 2;					// /2 for Ry
}
*/

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_force()
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_force");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_force");

	this->force_exx.create(GlobalC::ucell.nat, Ndim);
	for(int is=0; is<PARAM.inp.nspin; ++is)
	{
		this->exx_lri.cal_force({"","",std::to_string(is),"",""});
		for(std::size_t idim=0; idim<Ndim; ++idim) {
			for(const auto &force_item : this->exx_lri.force[idim]) {
				this->force_exx(force_item.first, idim) += std::real(force_item.second);
		} }
	}

	const double SPIN_multiple = std::map<int,double>{{1,2}, {2,1}, {4,1}}.at(PARAM.inp.nspin);				// why?
	const double frac = -2 * SPIN_multiple;		// why?
	this->force_exx *= frac;
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_force");
}


template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_stress()
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_stress");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_stress");

	this->stress_exx.create(Ndim, Ndim);
	for(int is=0; is<PARAM.inp.nspin; ++is)
	{
		this->exx_lri.cal_stress({"","",std::to_string(is),"",""});
		for(std::size_t idim0=0; idim0<Ndim; ++idim0) {
			for(std::size_t idim1=0; idim1<Ndim; ++idim1) {
				this->stress_exx(idim0,idim1) += std::real(this->exx_lri.stress(idim0,idim1));
		} }
	}

	const double SPIN_multiple = std::map<int,double>{{1,2}, {2,1}, {4,1}}.at(PARAM.inp.nspin);				// why?
	const double frac = 2 * SPIN_multiple / GlobalC::ucell.omega * GlobalC::ucell.lat0;		// why?
	this->stress_exx *= frac;

	ModuleBase::timer::tick("Exx_LRI", "cal_exx_stress");
}

template<typename Tdata>
std::vector<std::vector<int>> Exx_LRI<Tdata>::get_abfs_nchis() const
{
	std::vector<std::vector<int>> abfs_nchis;
	for (const auto& abfs_T : this->abfs)
	{
		std::vector<int> abfs_nchi_T;
		for (const auto& abfs_L : abfs_T)
			{ abfs_nchi_T.push_back(abfs_L.size()); }
		abfs_nchis.push_back(abfs_nchi_T);
	}
	return abfs_nchis;
}

#endif
