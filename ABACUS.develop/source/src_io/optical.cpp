#include "optical.h"
#include "../src_pw/tools.h"
#include "../src_pw/global.h"

bool Optical::opt_epsilon2 = false;
int  Optical::opt_nbands = 0;

Optical::Optical(){}
Optical::~Optical(){}

void Optical::cal_epsilon2(const int &nbands)
{
	if(!Optical::opt_epsilon2) return;
	TITLE("Optical","cal_epsilon2");
	timer::tick("Optical","cal_epsilon2");

	if(Optical::opt_nbands > NBANDS)
	{
		opt_nbands = NBANDS;
	}

	assert( wf.ekb!=0 );
	assert( wf.evc!=0 );

	cout << " begin to calculate the epsilon2." << endl;
	
	ofstream ofs;

	if(MY_RANK==0)
	{
		stringstream ss;
		ss << global_out_dir << "EPSILON2.dat";
		ofs.open( ss.str().c_str() );
	}

	const double de = 0.01; // unit: eV

	double maxe = wf.ekb[0][0];
	double mine = wf.ekb[0][0];

	for(int ik=0; ik<kv.nks; ik++)
	{
		for(int ib=0; ib<nbands; ib++)
		{
			maxe = std::max( wf.ekb[ik][ib], maxe);
			mine = std::min( wf.ekb[ik][ib], mine);
		}
	}

#ifdef __MPI
	Parallel_Reduce::gather_max_double_all( maxe );
	Parallel_Reduce::gather_min_double_all( mine );
#endif

	maxe *= Ry_to_eV;
	mine *= Ry_to_eV;

	double range = maxe - mine;
	int np = int(range / de) + 1; 
	int n_occ = static_cast<int>( (ucell.nelec+1)/2 + 1.0e-8 );

	cout << " n_occ = " << n_occ << endl;
	cout << " nbands = " << opt_nbands << endl;	
	cout << " unit is eV" << endl;
	cout << " maxe = " << maxe << endl;
	cout << " mine = " << mine << endl;
	cout << " energy range = " << maxe - mine << endl;
	cout << " de = " << de << " points = " << np << endl;

	OUT(ofs_running,"n_occ",n_occ);
	OUT(ofs_running,"nbands for optical",opt_nbands);
	OUT(ofs_running,"max energy",maxe);
	OUT(ofs_running,"min energy",mine);
	OUT(ofs_running,"energy range",maxe-mine);
	OUT(ofs_running,"de",de);
	OUT(ofs_running,"points",np);

	double *epsilon2 = new double[np];
	ZEROS(epsilon2, np);

	for(int ik=0; ik<kv.nks; ik++)
	{
		for(int iv=0; iv<n_occ; iv++)
		{
			const double ev = wf.ekb[ik][iv];
			for(int ic=n_occ; ic<opt_nbands; ic++)
			{
				const double ec = wf.ekb[ik][ic];
				const int ie = int((ec - ev)*Ry_to_eV/de);
				assert(ie < np);
				epsilon2[ie] += kv.wk[ik] * this->element_cvk(ik, iv, ic);
			}
		}
	}

#ifdef __MPI
	Parallel_Reduce::reduce_double_allpool( epsilon2, np);	
#endif

	if(MY_RANK==0)
	{
		ofs << np << endl;
		ofs << kv.nkstot << endl;
		for(int ie=0; ie<np; ie++)
		{
			ofs << ie*de << " " << epsilon2[ie] << endl;
		}
		delete[] epsilon2;
		ofs.close();
	}

	timer::tick("Optical","cal_epsilon2");
	return;
}


double Optical::element_cvk(const int &ik, const int &iv, const int &ic)
{
	double v=0.0;

	complex<double> tmp[3];
	ZEROS(tmp, 3);
	for(int ig=0; ig<kv.ngk[ik]; ig++)
	{
		const complex<double> uvc = conj( wf.evc[ik](ic,ig) ) * wf.evc[ik](iv, ig);
		tmp[0] += uvc * pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, ig), 0);
		tmp[1] += uvc * pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, ig), 1);
		tmp[2] += uvc * pw.get_GPlusK_cartesian_projection(ik, wf.igk(ik, ig), 2);
	}

	Parallel_Reduce::reduce_complex_double_pool( tmp,3 );

	v = norm(tmp[0]) + norm(tmp[1]) + norm(tmp[2]);

	return v;
}
