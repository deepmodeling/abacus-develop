#include "optical.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

bool Optical::opt_epsilon2 = false;
int  Optical::opt_nbands = 0;

Optical::Optical(){}
Optical::~Optical(){}

void Optical::cal_epsilon2(const int &nbands)
{
	if(!Optical::opt_epsilon2) return;
	ModuleBase::TITLE("Optical","cal_epsilon2");
	ModuleBase::timer::tick("Optical","cal_epsilon2");

	if(Optical::opt_nbands > GlobalV::NBANDS)
	{
		opt_nbands = GlobalV::NBANDS;
	}

	assert( GlobalC::wf.ekb!=0 );
	assert( GlobalC::wf.evc!=0 );

	std::cout << " begin to calculate the epsilon2." << std::endl;
	
	std::ofstream ofs;

	if(GlobalV::MY_RANK==0)
	{
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "EPSILON2.dat";
		ofs.open( ss.str().c_str() );
	}

	const double de = 0.01; // unit: eV

	double maxe = GlobalC::wf.ekb[0][0];
	double mine = GlobalC::wf.ekb[0][0];

	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		for(int ib=0; ib<nbands; ib++)
		{
			maxe = std::max( GlobalC::wf.ekb[ik][ib], maxe);
			mine = std::min( GlobalC::wf.ekb[ik][ib], mine);
		}
	}

#ifdef __MPI
	Parallel_Reduce::gather_max_double_all( maxe );
	Parallel_Reduce::gather_min_double_all( mine );
#endif

	maxe *= ModuleBase::Ry_to_eV;
	mine *= ModuleBase::Ry_to_eV;

	double range = maxe - mine;
	int np = int(range / de) + 1; 
	int n_occ = static_cast<int>( (GlobalC::CHR.nelec+1)/2 + 1.0e-8 );

	std::cout << " n_occ = " << n_occ << std::endl;
	std::cout << " nbands = " << opt_nbands << std::endl;	
	std::cout << " unit is eV" << std::endl;
	std::cout << " maxe = " << maxe << std::endl;
	std::cout << " mine = " << mine << std::endl;
	std::cout << " energy range = " << maxe - mine << std::endl;
	std::cout << " de = " << de << " points = " << np << std::endl;

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"n_occ",n_occ);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nbands for optical",opt_nbands);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max energy",maxe);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"min energy",mine);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"energy range",maxe-mine);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"de",de);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"points",np);

	double *epsilon2 = new double[np];
	ModuleBase::GlobalFunc::ZEROS(epsilon2, np);

	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		for(int iv=0; iv<n_occ; iv++)
		{
			const double ev = GlobalC::wf.ekb[ik][iv];
			for(int ic=n_occ; ic<opt_nbands; ic++)
			{
				const double ec = GlobalC::wf.ekb[ik][ic];
				const int ie = int((ec - ev)*ModuleBase::Ry_to_eV/de);
				assert(ie < np);
				epsilon2[ie] += GlobalC::kv.wk[ik] * this->element_cvk(ik, iv, ic);
			}
		}
	}

#ifdef __MPI
	Parallel_Reduce::reduce_double_allpool( epsilon2, np);	
#endif

	if(GlobalV::MY_RANK==0)
	{
		ofs << np << std::endl;
		ofs << GlobalC::kv.nkstot << std::endl;
		for(int ie=0; ie<np; ie++)
		{
			ofs << ie*de << " " << epsilon2[ie] << std::endl;
		}
		delete[] epsilon2;
		ofs.close();
	}

	ModuleBase::timer::tick("Optical","cal_epsilon2");
	return;
}


double Optical::element_cvk(const int &ik, const int &iv, const int &ic)
{
	double v=0.0;

	std::complex<double> tmp[3];
	ModuleBase::GlobalFunc::ZEROS(tmp, 3);
	for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
	{
		const std::complex<double> uvc = conj( GlobalC::wf.evc[ik](ic,ig) ) * GlobalC::wf.evc[ik](iv, ig);
		tmp[0] += uvc * GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, ig), 0);
		tmp[1] += uvc * GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, ig), 1);
		tmp[2] += uvc * GlobalC::pw.get_GPlusK_cartesian_projection(ik, GlobalC::wf.igk(ik, ig), 2);
	}

	Parallel_Reduce::reduce_complex_double_pool( tmp,3 );

	v = norm(tmp[0]) + norm(tmp[1]) + norm(tmp[2]);

	return v;
}
