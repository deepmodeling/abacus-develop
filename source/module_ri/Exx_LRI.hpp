//=======================
// AUTHOR : Peize Lin
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
#include "module_ri/conv_coulomb_pot_k-template.h"
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
void Exx_LRI<Tdata>::init(const MPI_Comm &mpi_comm_in, const K_Vectors &kv_in)
{
	ModuleBase::TITLE("Exx_LRI","init");
	ModuleBase::timer::tick("Exx_LRI", "init");

//	if(GlobalC::exx_info.info_global.separate_loop)
//	{
//		Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::No;
//		Hexx_para.mixing_beta = 0;
//	}
//	else
//	{
//		if("plain"==GlobalC::CHR.mixing_mode)
//			Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Plain;
//		else if("pulay"==GlobalC::CHR.mixing_mode)
//			Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Pulay;
//		else
//			throw std::invalid_argument("exx mixing error. exx_separate_loop==false, mixing_mode!=plain or pulay");
//		Hexx_para.mixing_beta = GlobalC::CHR.mixing_beta;
//	}

    this->mpi_comm = mpi_comm_in;
    this->p_kv = &kv_in;

	this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs( GlobalC::ORB, this->info.kmesh_times );

//	#ifdef __MPI
//	Exx_Abfs::Util::bcast( this->info.files_abfs, 0, this->mpi_comm );
//	#endif

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs_same_atom = Exx_Abfs::Construct_Orbs::abfs_same_atom( this->lcaos, this->info.kmesh_times, this->info.pca_threshold );
	if(this->info.files_abfs.empty())
		this->abfs = abfs_same_atom;
	else
		this->abfs = Exx_Abfs::IO::construct_abfs( abfs_same_atom, GlobalC::ORB, this->info.files_abfs, this->info.kmesh_times );
	Exx_Abfs::Construct_Orbs::print_orbs_size(this->abfs, GlobalV::ofs_running);

	auto get_ccp_parameter = [this]() -> std::map<std::string,double>
	{
		switch(this->info.ccp_type)
		{
			case Conv_Coulomb_Pot_K::Ccp_Type::Ccp:
				return {};
			case Conv_Coulomb_Pot_K::Ccp_Type::Hf:
				return {};
			case Conv_Coulomb_Pot_K::Ccp_Type::Hse:
				return {{"hse_omega", this->info.hse_omega}};
			default:
				throw std::domain_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));	break;
		}
	};
    this->abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(this->abfs, this->info.ccp_type, get_ccp_parameter(), this->info.ccp_rmesh_times, this->p_kv->nkstot_full);


	for( size_t T=0; T!=this->abfs.size(); ++T )
		GlobalC::exx_info.info_ri.abfs_Lmax = std::max( GlobalC::exx_info.info_ri.abfs_Lmax, static_cast<int>(this->abfs[T].size())-1 );

	this->cv.set_orbitals(
		this->lcaos, this->abfs, this->abfs_ccp,
		this->info.kmesh_times, this->info.ccp_rmesh_times );

	ModuleBase::timer::tick("Exx_LRI", "init");
}

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_ions()
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_ions");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_ions");

//	init_radial_table_ions( cal_atom_centres_core(atom_pairs_core_origin), atom_pairs_core_origin );

//	this->m_abfsabfs.init_radial_table(Rradial);
//	this->m_abfslcaos_lcaos.init_radial_table(Rradial);

	std::vector<TA> atoms(GlobalC::ucell.nat);
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		atoms[iat] = iat;
	std::map<TA,TatomR> atoms_pos;
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		atoms_pos[iat] = RI_Util::Vector3_to_array3( GlobalC::ucell.atoms[ GlobalC::ucell.iat2it[iat] ].tau[ GlobalC::ucell.iat2ia[iat] ] );
	const std::array<TatomR,Ndim> latvec
		= {RI_Util::Vector3_to_array3(GlobalC::ucell.a1),
		   RI_Util::Vector3_to_array3(GlobalC::ucell.a2),
		   RI_Util::Vector3_to_array3(GlobalC::ucell.a3)};
	const std::array<Tcell,Ndim> period = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};

	this->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

	// std::max(3) for gamma_only, list_A2 should contain cell {-1,0,1}. In the future distribute will be neighbour.
	const std::array<Tcell,Ndim> period_Vs = LRI_CV_Tools::cal_latvec_range<Tcell>(1+this->info.ccp_rmesh_times);	
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As_Vs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);

	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>
		Vs = this->cv.cal_Vs(
			list_As_Vs.first, list_As_Vs.second[0],
			{{"writable_Vws",true}});
	this->cv.Vws = LRI_CV_Tools::get_CVws(Vs);
	this->exx_lri.set_Vs(std::move(Vs), this->info.V_threshold);

	if(GlobalV::CAL_FORCE || GlobalV::CAL_STRESS)
	{
		std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>
			dVs = this->cv.cal_dVs(
				list_As_Vs.first, list_As_Vs.second[0],
				{{"writable_dVws",true}});
		this->cv.dVws = LRI_CV_Tools::get_dCVws(dVs);
		this->exx_lri.set_dVs(std::move(dVs), this->info.V_grad_threshold);
	}

	const std::array<Tcell,Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As_Cs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);

	std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>, std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>>
		Cs_dCs = this->cv.cal_Cs_dCs(
			list_As_Cs.first, list_As_Cs.second[0],
			{{"cal_dC",GlobalV::CAL_FORCE||GlobalV::CAL_STRESS},
			 {"writable_Cws",true}, {"writable_dCws",true}, {"writable_Vws",false}, {"writable_dVws",false}});
	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> &Cs = std::get<0>(Cs_dCs);
	this->cv.Cws = LRI_CV_Tools::get_CVws(Cs);
	this->exx_lri.set_Cs(std::move(Cs), this->info.C_threshold);

	if(GlobalV::CAL_FORCE || GlobalV::CAL_STRESS)
	{
		std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3> &dCs = std::get<1>(Cs_dCs);
		this->cv.dCws = LRI_CV_Tools::get_dCVws(dCs);
		this->exx_lri.set_dCs(std::move(dCs), this->info.C_grad_threshold);
	}
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_ions");
}

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_elec(const std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>> &Ds, const Parallel_Orbitals &pv)
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_elec");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_elec");

	const std::vector<std::tuple<std::set<TA>, std::set<TA>>> judge = RI_2D_Comm::get_2D_judge(pv);

	this->exx_lri.set_csm_threshold(this->info.cauchy_threshold);

	this->Hexxs.resize(GlobalV::NSPIN);
	this->Eexx = 0;
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		if(!(GlobalV::CAL_FORCE || GlobalV::CAL_STRESS))
		{
			this->exx_lri.set_Ds(Ds[is], this->info.dm_threshold);
			this->exx_lri.cal_Hs();
		}
		else
		{
			this->exx_lri.set_Ds(Ds[is], this->info.dm_threshold, std::to_string(is));
			this->exx_lri.cal_Hs({"","",std::to_string(is)});
		}
		this->Hexxs[is] = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
			this->mpi_comm, std::move(this->exx_lri.Hs), std::get<0>(judge[is]), std::get<1>(judge[is]));
		this->Eexx += this->exx_lri.energy;
		post_process_Hexx(this->Hexxs[is]);
	}
	this->Eexx = post_process_Eexx(this->Eexx);
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
Tdata Exx_LRI<Tdata>::post_process_Eexx( const Tdata &Eexx_in ) const
{
	ModuleBase::TITLE("Exx_LRI","post_process_Eexx");
	const Tdata SPIN_multiple = std::map<int,Tdata>{{1,2}, {2,1}, {4,1}}.at(GlobalV::NSPIN);				// why?
	const Tdata frac = - SPIN_multiple;
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
	energy *= SPIN_multiple.at(GlobalV::NSPIN);			// ?
	energy /= 2;					// /2 for Ry
}
*/

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_force()
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_force");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_force");
		
	this->exx_lri.set_csm_threshold(this->info.cauchy_force_threshold);

	this->force_exx.create(GlobalC::ucell.nat, Ndim);
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		this->exx_lri.cal_force({"","",std::to_string(is),"",""});
		for(std::size_t idim=0; idim<Ndim; ++idim)
			for(const auto &force_item : this->exx_lri.force[idim])
				this->force_exx(force_item.first, idim) += std::real(force_item.second);
	}

	const double SPIN_multiple = std::map<int,double>{{1,2}, {2,1}, {4,1}}.at(GlobalV::NSPIN);				// why?
	const double frac = -2 * SPIN_multiple;		// why?
	this->force_exx *= frac;
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_force");
}


template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_stress()
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_stress");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_stress");
		
	this->exx_lri.set_csm_threshold(this->info.cauchy_stress_threshold);

	this->stress_exx.create(Ndim, Ndim);
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		this->exx_lri.cal_stress({"","",std::to_string(is),"",""});
		for(std::size_t idim0=0; idim0<Ndim; ++idim0)
			for(std::size_t idim1=0; idim1<Ndim; ++idim1)
				this->stress_exx(idim0,idim1) += std::real(this->exx_lri.stress(idim0,idim1));
	}

	const double SPIN_multiple = std::map<int,double>{{1,2}, {2,1}, {4,1}}.at(GlobalV::NSPIN);				// why?
	const double frac = 2 * SPIN_multiple / GlobalC::ucell.omega * GlobalC::ucell.lat0;		// why?
	this->stress_exx *= frac;

	ModuleBase::timer::tick("Exx_LRI", "cal_exx_stress");
}


#endif
