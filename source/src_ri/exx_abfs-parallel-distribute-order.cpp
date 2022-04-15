#include "exx_abfs-parallel-distribute-order.h"

#include "../src_pw/global.h"
#include "abfs.h"
#include <cmath>

#ifdef __MPI
std::vector<std::pair<size_t,size_t>> Exx_Abfs::Parallel::Distribute::Order::distribute(
	const double rmesh_times )
{
	const std::vector<Abfs::Vector3_Order<int>> Coulomb_potential_boxes = Abfs::get_Coulomb_potential_boxes(rmesh_times);
	auto neighbour = [&](const size_t iat1, const size_t iat2) -> bool
	{
		const int it1 = GlobalC::ucell.iat2it[iat1];
		const int it2 = GlobalC::ucell.iat2it[iat2];
		const ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[it1].tau[GlobalC::ucell.iat2ia[iat1]];
		const ModuleBase::Vector3<double> tau2 = GlobalC::ucell.atoms[it2].tau[GlobalC::ucell.iat2ia[iat2]];
		const double Rcut = std::min( GlobalC::ORB.Phi[it1].getRcut()*rmesh_times+GlobalC::ORB.Phi[it2].getRcut(), GlobalC::ORB.Phi[it1].getRcut()+GlobalC::ORB.Phi[it2].getRcut()*rmesh_times );
		for(const ModuleBase::Vector3<int> box2 : Coulomb_potential_boxes)
		{
			const double R = (-tau1 + tau2 + box2 * GlobalC::ucell.latvec).norm();
			if(R*GlobalC::ucell.lat0 < Rcut)
				return true;
		}
		return false;
	};
	
	std::vector<std::pair<size_t,size_t>> pairs;
	for( size_t iat1=0; iat1<GlobalC::ucell.nat; ++iat1 )
		for( size_t iat2=iat1; iat2<GlobalC::ucell.nat; ++iat2 )
			if(neighbour(iat1,iat2))
				pairs.push_back( {iat1,iat2} );
			
	const size_t work = std::ceil(static_cast<double>(pairs.size())/GlobalV::NPROC);
	return std::vector<std::pair<size_t,size_t>>{ 
		std::min( pairs.begin()+ GlobalV::MY_RANK   *work, pairs.end() ),
		std::min( pairs.begin()+(GlobalV::MY_RANK+1)*work, pairs.end() )};
}
#endif