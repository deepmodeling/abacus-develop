#ifndef EXX_ABFS_CONSTRUCT_ORBS_TEST_H
#define EXX_ABFS_CONSTRUCT_ORBS_TEST_H

#include <fstream>
#include <vector>
#include <string>

static void print_orbs( const std::vector<std::vector<std::vector<std::vector<double>>>> & orbs, const std::string & file_name )
{
	std::ofstream ofsN(file_name);
	for( size_t T=0; T!=orbs.size(); ++T )
		for( size_t L=0; L!=orbs[T].size(); ++L )
			for( size_t N=0; N!=orbs[T][L].size(); ++N )
			{
				ofsN<<T<<"\t"<<L<<"\t"<<N<<std::endl;
				for( size_t ir=0; ir!=orbs[T][L][N].size(); ++ir )
					ofsN<<orbs[T][L][N][ir]<<"\t";
				ofsN<<std::endl;
			}	
	ofsN.close();
}

static void print_orbs( const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> & orbs, const std::string & file_name )
{
	std::ofstream ofsN(file_name);
	for( size_t T=0; T!=orbs.size(); ++T )
		for( size_t L=0; L!=orbs[T].size(); ++L )
			for( size_t N=0; N!=orbs[T][L].size(); ++N )
			{
				ofsN<<T<<"\t"<<L<<"\t"<<N<<std::endl;
				for( size_t ir=0; ir!=orbs[T][L][N].getNr(); ++ir )
					ofsN<<orbs[T][L][N].getPsi(ir)<<"\t";
				ofsN<<std::endl;
			}	
	ofsN.close();
}

#endif	// EXX_ABFS_CONSTRUCT_ORBS_TEST_H