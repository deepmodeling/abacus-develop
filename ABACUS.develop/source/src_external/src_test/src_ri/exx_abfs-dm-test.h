#include <map>
#include <vector>
#include <set>
#include <complex>
#include "src_global/global_function.h"
#include "module_base/matrix.h"
#include "src_pw/global.h"

class Exx_Abfs_DM_Test
{
public:
	static map<size_t,map<size_t,vector<ComplexMatrix>>> cal_DMk_raw_readfile( const set<pair<size_t,size_t>> &atom_pairs );
	static vector<vector<complex<double>>> read_wfc( const string &file_name );
	static matrix read_wg( const string &file_name );
};

map<size_t,map<size_t,vector<ComplexMatrix>>> 
Exx_Abfs_DM_Test::cal_DMk_raw_readfile( const set<pair<size_t,size_t>> &atom_pairs )
{
	static int istep=-1;	++istep;	
	TITLE("cal_DMk_raw_readfile_"+TO_STRING(istep));
	matrix wf_wg = read_wg( "wf.wg/wf.wg_"+TO_STRING(istep) );
	vector<vector<vector<complex<double>>>> wfc(kv.nks);
	for( size_t ik=0; ik!=kv.nks; ++ik )
		wfc[ik] = read_wfc( "hvec/hvec_"+TO_STRING(istep)+"_"+TO_STRING(ik) );	
	
	
	const double SPIN_multiple = 0.5*NSPIN;
	
	map<size_t,map<size_t,vector<ComplexMatrix>>> DMk_raw;
	for( const pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = ucell.iat2it[iat1];
		const size_t it2 = ucell.iat2it[iat2];
		const size_t ia1 = ucell.iat2ia[iat1];
		const size_t ia2 = ucell.iat2ia[iat2];

		DMk_raw[iat1][iat2] = vector<ComplexMatrix>( kv.nks, {ucell.atoms[it1].nw,ucell.atoms[it2].nw} );
		for( size_t ik=0; ik!=kv.nks; ++ik )
		{
			for( size_t iw1=0; iw1!=ucell.atoms[it1].nw; ++iw1 )
			{
				for( size_t iw2=0; iw2!=ucell.atoms[it2].nw; ++iw2 )
				{
					for( size_t ib=0; ib!=NBANDS; ++ib )
					{
							DMk_raw[iat1][iat2][ik](iw1,iw2) += wf_wg(ik,ib) * conj(wfc[ik][ib][ucell.itiaiw2iwt(it1,ia1,iw1)]) * wfc[ik][ib][ucell.itiaiw2iwt(it2,ia2,iw2)];
					}
				}
			}
			DMk_raw[iat1][iat2][ik] *= SPIN_multiple;
		}
	}
	return DMk_raw;
}



vector<vector<complex<double>>> Exx_Abfs_DM_Test::read_wfc( const string &file_name )
{
	vector<vector<complex<double>>> wfc(NBANDS,vector<complex<double>>(NLOCAL));
	ifstream ifs(file_name);
	for( size_t iw=0; iw!=NLOCAL; ++iw )
	{
		for( size_t ib=0; ib!=NBANDS; ++ib )
		{
			double a,b;
			ifs>>a>>b;
			wfc[ib][iw] = {a,b};
		}
	}
	ifs.close();
	return wfc;
}



matrix Exx_Abfs_DM_Test::read_wg( const string &file_name )
{
	matrix wf_wg(kv.nks,NBANDS);
	ifstream ifs(file_name);
	for( size_t ik=0; ik!=kv.nks; ++ik )
	{
		for( size_t ib=0; ib!=NBANDS; ++ib )
		{
			ifs>>wf_wg(ik,ib);
		}
	}
	ifs.close();
	return wf_wg;
}
