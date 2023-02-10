#include <map>
#include <vector>
#include <set>
#include <complex>
#include "module_base/global_function.h"
#include "module_base/matrix.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

class Exx_Abfs_DM_Test
{
public:
	static std::map<size_t,std::map<size_t,std::vector<ModuleBase::ComplexMatrix>>> cal_DMk_raw_readfile( const set<std::pair<size_t,size_t>> &atom_pairs );
	static std::vector<std::vector<std::complex<double>>> read_wfc( const std::string &file_name );
	static ModuleBase::matrix read_wg( const std::string &file_name );
};

std::map<size_t,std::map<size_t,std::vector<ModuleBase::ComplexMatrix>>> 
Exx_Abfs_DM_Test::cal_DMk_raw_readfile( const set<std::pair<size_t,size_t>> &atom_pairs )
{
	static int istep=-1;	++istep;	
	ModuleBase::TITLE("cal_DMk_raw_readfile_"+ModuleBase::GlobalFunc::TO_STRING(istep));
	ModuleBase::matrix wf_wg = read_wg( "wf.wg/wf.wg_"+ModuleBase::GlobalFunc::TO_STRING(istep) );
	std::vector<std::vector<std::vector<std::complex<double>>>> wfc(GlobalC::kv.nks);
	for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		wfc[ik] = read_wfc( "hvec/hvec_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(ik) );	
	
	
	const double SPIN_multiple = 0.5*GlobalV::NSPIN;
	
	std::map<size_t,std::map<size_t,std::vector<ModuleBase::ComplexMatrix>>> DMk_raw;
	for( const std::pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t it2 = GlobalC::ucell.iat2it[iat2];
		const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
		const size_t ia2 = GlobalC::ucell.iat2ia[iat2];

		DMk_raw[iat1][iat2] = std::vector<ModuleBase::ComplexMatrix>( GlobalC::kv.nks, {GlobalC::ucell.atoms[it1].nw,GlobalC::ucell.atoms[it2].nw} );
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		{
			for( size_t iw1=0; iw1!=GlobalC::ucell.atoms[it1].nw; ++iw1 )
			{
				for( size_t iw2=0; iw2!=GlobalC::ucell.atoms[it2].nw; ++iw2 )
				{
					for( size_t ib=0; ib!=GlobalV::NBANDS; ++ib )
					{
							DMk_raw[iat1][iat2][ik](iw1,iw2) += wf_wg(ik,ib) * conj(wfc[ik][ib][GlobalC::ucell.itiaiw2iwt(it1,ia1,iw1)]) * wfc[ik][ib][GlobalC::ucell.itiaiw2iwt(it2,ia2,iw2)];
					}
				}
			}
			DMk_raw[iat1][iat2][ik] *= SPIN_multiple;
		}
	}
	return DMk_raw;
}



std::vector<std::vector<std::complex<double>>> Exx_Abfs_DM_Test::read_wfc( const std::string &file_name )
{
	std::vector<std::vector<std::complex<double>>> wfc(GlobalV::NBANDS,std::vector<std::complex<double>>(GlobalV::NLOCAL));
	std::ifstream ifs(file_name);
	for( size_t iw=0; iw!=GlobalV::NLOCAL; ++iw )
	{
		for( size_t ib=0; ib!=GlobalV::NBANDS; ++ib )
		{
			double a,b;
			ifs>>a>>b;
			wfc[ib][iw] = {a,b};
		}
	}
	ifs.close();
	return wfc;
}



ModuleBase::matrix Exx_Abfs_DM_Test::read_wg( const std::string &file_name )
{
	ModuleBase::matrix wf_wg(GlobalC::kv.nks,GlobalV::NBANDS);
	std::ifstream ifs(file_name);
	for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
	{
		for( size_t ib=0; ib!=GlobalV::NBANDS; ++ib )
		{
			ifs>>wf_wg(ik,ib);
		}
	}
	ifs.close();
	return wf_wg;
}
