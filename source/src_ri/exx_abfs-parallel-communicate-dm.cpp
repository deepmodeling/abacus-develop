#include "exx_abfs-parallel-communicate-dm.h"
#include "../src_lcao/global_fp.h"
#include "../src_pw/global.h"
#include "../src_lcao/record_adj.h"

#include "../src_external/src_test/src_ri/exx_lcao-test.h"
#include "lcao_nnr.h"
//#include <gperftools/profiler.h>

/*
void Exx_Abfs::Parallel::Communicate::DM::set_atom_in_exx( const set<std::pair<size_t,size_t>> &H_atom_pairs_core )
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::DM::set_atom_in_exx");
	
	atom_in_exx.row.resize(GlobalC::ucell.nat);
	atom_in_exx.col.resize(GlobalC::ucell.nat);
	for( const auto std::pair : H_atom_pairs_core )
	{
		atom_in_exx.row[std::pair.first]  = true;
		atom_in_exx.col[std::pair.second] = true;
	}
}
*/

/*
std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>
Exx_Abfs::Parallel::Communicate::DM::m2D_to_a2Dexx( 
	const ModuleBase::matrix &DM_m2D, 
	const int iwt1_index_begin, 
	const int iwt2_index_begin) const
{
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>> DM_a2Dexx;
	for( int iwt1_m2D=0; iwt1_m2D<DM_m2D.nr; ++iwt1_m2D )
	{
		const int iwt1 = GlobalC::GridT.trace_lo[iwt1_index_begin+iwt1_m2D];
		const size_t iat1 = GlobalC::ucell.iwt2iat[iwt1];
		if(!atom_in_exx.row[iat1])	continue;
		const size_t iw1 = GlobalC::ucell.iwt2iw[iwt1];
		for( int iwt2_m2D=0; iwt2_m2D<DM_m2D.nc; ++iwt2_m2D )
		{
			const int iwt1 = GlobalC::GridT.trace_lo[iwt2_index_begin+iwt2_m2D];
			const size_t iat2 = GlobalC::ucell.iwt2iat[iwt2];
			if(!atom_in_exx.col[iat2])	continue;
			const size_t iw2 = GlobalC::ucell.iwt2iw[iwt2];
			
			try
			{ 
				DM_a2Dexx.at(iat1).at(iat2).at(Abfs::Vector3_Order<int>(0,0,0))(iw1,iw2) = DM_m2D(iwt1_m2D,iwt2_m2D); 
			}
			catch(const std::out_of_range&)
			{ 
				DM_a2Dexx[iat1][iat2][Abfs::Vector3_Order<int>(0,0,0)].creat( 
					GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw, 
					GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw); 
				DM_a2Dexx[iat1][iat2][Abfs::Vector3_Order<int>(0,0,0)](iw1,iw2) = DM_m2D(iwt1_m2D,iwt2_m2D);
			}
		}
	}
	return DM_a2Dexx;
}


Exx_Abfs::Parallel::Communicate::DM::a2Dexx_to_exx( 
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>> & m_exx, 
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>> & m_a2Dexx ) const
{
	for( size_t is=0; is!=GlobalV::NSPIN; ++is )
	{
		for( auto & m_a2DexxA : m_a2Dexx[is] )
		{
			const size_t iat1 = m_a2DexxA.first;
			for( auto & m_a2DexxB : m_a2DexxA.second )
			{
				const size_t iat2 = m_a2DexxB.first;
				for( auto & m_a2DexxC : m_a2DexxB.second )
				{
					const Abfs::Vector3_Order<int> & box2 = m_a2DexxC.first;
					if( m_a2DexxC.second.absmax() < threshold )	continue;
					try                       { m_exx[is].at(iat1).at(iat2).at(box2) +=           m_a2DexxC.second ; }
					catch(std::out_of_range&) { m_exx[is]   [iat1]   [iat2]   [box2]  = std::move(m_a2DexxC.second); }
				}
			}
		}
		m_a2Dexx[is].clear();
	}
}


Exx_Abfs::Parallel::Communicate::DM::f( 
	const ModuleBase::matrix &DM_m2D, 
	const int iwt1_index_begin, 
	const int iwt2_index_begin)
{
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>> & DM_a2Dexx = 
		m2D_to_a2Dexx( DM_m2D, iwt1_index_begin, iwt2_index_begin );
	a2Dexx_to_exx( this->DMr, DM_a2Dexx );
}
*/

void Exx_Abfs::Parallel::Communicate::DM::cal_DM( 
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
	const set<std::pair<size_t,size_t>> &H_atom_pairs_core,
    const double threshold,
    complex<double>*** wfc_k_grid,
    double*** DM,
    double** DM_R)
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::DM::cal_DM");
	
std::ofstream ofs_time("time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start;
//gettimeofday( &t_start, NULL);
#if false
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DM_grid = LOC_to_grid( Born_von_Karman_period, threshold, DM, DM_R );
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::DM::LOC_to_grid\t"<<time_during(t_start)<<std::endl;
ofs_matrixes( "DM_grid_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK), DM_grid );

//gettimeofday( &t_start, NULL);	
	MPI_Barrier( MPI_COMM_WORLD );
//ofs_time<<"TIME@ MPI_Barrier\t"<<time_during(t_start)<<std::endl;

//gettimeofday( &t_start, NULL);	
	Allreduce allreduce( MPI_COMM_WORLD, DM_grid, Born_von_Karman_period, H_atom_pairs_core );
	this->DMr = allreduce.grid_to_exx();
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::DM::Allreduce::grid_to_exx\t"<<time_during(t_start)<<std::endl;
#else
	auto cal_dm_my = [&]() -> std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>>
	{
		std::vector<Abfs::Vector3_Order<int>> Born_von_Karman_boxes;
		for( int ix=0; ix<Born_von_Karman_period.x; ++ix )
			for( int iy=0; iy<Born_von_Karman_period.y; ++iy )
				for( int iz=0; iz<Born_von_Karman_period.z; ++iz )
					Born_von_Karman_boxes.push_back({ix,iy,iz});
				
		Exx_Abfs::DM dm_my;
		dm_my.flag_mix = false;
		dm_my.cal_DM( H_atom_pairs_core, Born_von_Karman_boxes, wfc_k_grid );
		
		std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DM_grid(GlobalV::NSPIN);
		for( const auto &DMrA : dm_my.DMr )
		{
			const size_t iat1 = DMrA.first;
			for( const auto &DMrB : DMrA.second )
			{
				const size_t iat2 = DMrB.first;
				for( const auto &DMrC : DMrB.second )
				{
					const auto box2 = DMrC.first;
					for( size_t is=0; is!=DMrC.second.size(); ++is )
						DM_grid[is][iat1][iat2][box2] = std::move(DMrC.second[is]);
				}
			}
		}
		return DM_grid;
	};
//gettimeofday( &t_start, NULL);	
	this->DMr = cal_dm_my();
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::DM::Allreduce::cal_dm_my\t"<<time_during(t_start)<<std::endl;
#endif	
//ofs_time.close();
ofs_matrixes( "DMr_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK), DMr );
}

std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>>
Exx_Abfs::Parallel::Communicate::DM::LOC_to_grid( 
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
    const double threshold,
    double*** DM,
    double** DM_R) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::DM::LOC_to_grid");
	
	const double SPIN_multiple = 0.5*GlobalV::NSPIN;
	
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DM_grid(GlobalV::NSPIN);
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
{
	std::ofstream ofs("LOC.DM_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		for( int i1=0; i1!=GlobalC::GridT.lgd; ++i1 )
		{
			for( int i2=0; i2!=GlobalC::GridT.lgd; ++i2 )
				ofs<<DM[is][i1][i2]<<"\t";
			ofs<<std::endl;
		}
		ofs<<std::endl;
	}
	ofs.close();
}		
		for( size_t is=0; is!=GlobalV::NSPIN; ++is )
		{			
			for( int iat1=0, iwt1_index=0; iat1!=GlobalC::ucell.nat; ++iat1 )
			{
				if(!GlobalC::GridT.in_this_processor[iat1])	continue;
				const int nw1 = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw;
				for( int iat2=0, iwt2_index=0; iat2!=GlobalC::ucell.nat; ++iat2 )
				{
					if(!GlobalC::GridT.in_this_processor[iat2])	continue;
					const int nw2 = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw;
					
					ModuleBase::matrix DM_grid_2D(nw1,nw2);
					for( int iw1=0; iw1!=nw1; ++iw1 )
					{
						for( int iw2=0; iw2!=nw2; ++iw2 )
						{
							DM_grid_2D(iw1,iw2) = DM[is][iwt1_index+iw1][iwt2_index+iw2];
						}
					}
					if( DM_grid_2D.absmax() * SPIN_multiple >= threshold )
						DM_grid[is][iat1][iat2][{0,0,0}] = DM_grid_2D * SPIN_multiple;
					else
						DM_grid[is][iat1][iat2][{0,0,0}];					// to erase atom_unset
					iwt2_index += nw2;
				}
				iwt1_index += nw1;
			}
	
			/*
			for( int iwt1_grid=0; iwt1_grid<GlobalC::GridT.lgd; ++iwt1_grid )
			{				
				const int iwt1 = GlobalC::ParaO.MatrixInfo.row_set[iwt1_grid];
				const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
				const int iw1  = GlobalC::ucell.iwt2iw[iwt1];
std::cout<<iwt1_grid<<"\t"<<iwt1<<"\t"<<iat1<<"\t"<<iw1<<std::endl;
				for( int iwt2_grid=0; iwt2_grid<GlobalC::GridT.lgd; ++iwt2_grid )
				{				
					const int iwt2 = GlobalC::ParaO.MatrixInfo.col_set[iwt2_grid];
					const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
					const int iw2  = GlobalC::ucell.iwt2iw[iwt2];
std::cout<<"\t"<<iwt2_grid<<"\t"<<iwt2<<"\t"<<iat2<<"\t"<<iw2<<std::endl;
					try
					{
						DM_grid[is].at(iat1).at(iat2).at({0,0,0})(iw1,iw2) = DM[is][iwt1_grid][iwt2_grid];
					}
					catch(const std::out_of_range&)
					{
						DM_grid[is][iat1][iat2][{0,0,0}].create(
							GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw,
							GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw);
						DM_grid[is][iat1][iat2][{0,0,0}](iw1,iw2) = DM[is][iwt1_grid][iwt2_grid];
					}
				}
			}*/
		}
	}
	else
	{	
std::ofstream ofs_LOC_DM("LOC.DM_R_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
for( int i=0; i<100; ++i )
	ofs_LOC_DM<<DM_R[0][i]<<"\t";
ofs_LOC_DM<<std::endl<<std::endl;

		Record_adj RA;
		RA.for_grid(GlobalC::GridT);
//std::cout<<__FILE__<<__LINE__<<std::endl;

		for( size_t is=0; is!=GlobalV::NSPIN; ++is )
		{
			for( size_t iat1=0; iat1<GlobalC::ucell.nat; ++iat1 )
			{
				if(!GlobalC::GridT.in_this_processor[iat1])	continue;
				int iw_index = 0;
				for( size_t iat2_2D=0; iat2_2D<RA.na_each[iat1]; ++iat2_2D )
				{
					const size_t iat2 = GlobalC::ucell.itia2iat( RA.info[iat1][iat2_2D][3], RA.info[iat1][iat2_2D][4] );
					const Abfs::Vector3_Order<int> box2( RA.info[iat1][iat2_2D][0], RA.info[iat1][iat2_2D][1], RA.info[iat1][iat2_2D][2] );
					const Abfs::Vector3_Order<int> boxp2 = box2%Born_von_Karman_period;
					const int nw1 = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw;
					const int nw2 = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw;	
{
	ofs_LOC_DM<<"@\t"<<iat1<<"\t"<<iat2<<"\t"<<box2<<std::endl;
	for( int iw1=0; iw1!=nw1; ++iw1 )
	{
		for( int iw2=0; iw2!=nw2; ++iw2 )
			ofs_LOC_DM<<DM_R[is][GlobalC::GridT.nlocstartg[iat1]+iw_index+iw1*nw2+iw2]<<"\t";
		ofs_LOC_DM<<std::endl;
	}
	ofs_LOC_DM<<std::endl;
}			
					if( !ModuleBase::GlobalFunc::MAP_EXIST( DM_grid[is], iat1, iat2, boxp2 ) )
					{					
						ModuleBase::matrix DM_grid_2D(nw1,nw2,false);
						memcpy( DM_grid_2D.c, DM_R[is]+GlobalC::GridT.nlocstartg[iat1]+iw_index, sizeof(double)*(nw1*nw2) );
						if( DM_grid_2D.absmax() * SPIN_multiple >= threshold )
							DM_grid[is][iat1][iat2][boxp2] = DM_grid_2D * SPIN_multiple;
						else
							DM_grid[is][iat1][iat2][boxp2];					// to erase atom_unset
					}
					iw_index += nw1*nw2;															   
				}
			}
		}
ofs_LOC_DM.close();
	}
	return DM_grid;
}

