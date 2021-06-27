//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================

#include"vdwd2.h"
#include"src_global/global_function.h"
#include"module_base/constants.h"
#include<cmath>

Vdwd2::Vdwd2(const UnitCell_pseudo &unit_in, Vdwd2_Parameters &para_in):
	ucell(unit_in),
	para(para_in){}

void Vdwd2::cal_energy()
{
    TITLE("Vdwd2","energy");
	para.initset(ucell);

	energy = 0;
	for( int it1=0; it1!=ucell.ntype; ++it1 )
	{
		for( int it2=0; it2!=ucell.ntype; ++it2 )
		{
			const double C6_product = sqrt( para.C6.at(ucell.atoms[it1].psd) * para.C6.at(ucell.atoms[it2].psd) )/pow(ucell.lat0,6) ;
			const double R0_sum = ( para.R0.at(ucell.atoms[it1].psd) + para.R0.at(ucell.atoms[it2].psd) )/ucell.lat0;
			if(!R0_sum)
				WARNING_QUIT("Input", "R0_sum can not be 0");		
			for( int ia1=0; ia1!=ucell.atoms[it1].na; ++ia1 )
			{
				for( int ia2=0; ia2!=ucell.atoms[it2].na; ++ia2 )
				{
					Vector3<int> ilat_loop;
					for( ilat_loop.x = -para.period.x/2; ilat_loop.x <= (para.period.x-1)/2; ++ilat_loop.x )
						for( ilat_loop.y = -para.period.y/2; ilat_loop.y <= (para.period.y-1)/2; ++ilat_loop.y )
							for( ilat_loop.z = -para.period.z/2; ilat_loop.z <= (para.period.z-1)/2; ++ilat_loop.z )
							{
								if( (!( ilat_loop.x || ilat_loop.y || ilat_loop.z )) && (it1==it2) && (ia1==ia2) )
									continue;
								const Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
								const Vector3<double> tau2 = ucell.atoms[it2].tau[ia2] + ilat_loop * ucell.latvec;
								const double r_sqr = (tau1 - tau2).norm2();
								const double r = sqrt(r_sqr);
								const double tmp_damp_recip = 1+ exp( -para.damping* (r/R0_sum-1) );
								energy -= C6_product/ pow(r_sqr,3)/ tmp_damp_recip/ 2;
							} // end for ilat_loop
				} // end for ia2
			} // end for ia1
		} // end for it2
	} // end for it1
	energy *= para.scaling;
}

void Vdwd2::cal_force()
{
    TITLE("Vdwd2","force");
	para.initset(ucell);

	force.clear();
	force.resize(ucell.nat);
	
	for( int it1=0; it1!=ucell.ntype; ++it1 )
	{
		for( int it2=0; it2!=ucell.ntype; ++it2 )
		{
			const double C6_product = sqrt( para.C6.at(ucell.atoms[it1].psd) * para.C6.at(ucell.atoms[it2].psd) )/pow(ucell.lat0,6);
			const double R0_sum = ( para.R0.at(ucell.atoms[it1].psd) + para.R0.at(ucell.atoms[it2].psd) )/ucell.lat0;
			if(!R0_sum)
				WARNING_QUIT("Input", "R0_sum can not be 0");
			for( int ia1=0; ia1!=ucell.atoms[it1].na; ++ia1 )
			{
				for( int ia2=0; ia2!=ucell.atoms[it2].na; ++ia2 )
				{
					Vector3<int> ilat_loop;
					for( ilat_loop.x = -para.period.x/2; ilat_loop.x <= (para.period.x-1)/2; ++ilat_loop.x )
						for( ilat_loop.y = -para.period.y/2; ilat_loop.y <= (para.period.y-1)/2; ++ilat_loop.y )
							for( ilat_loop.z = -para.period.z/2; ilat_loop.z <= (para.period.z-1)/2; ++ilat_loop.z )
							{
								if( (!( ilat_loop.x || ilat_loop.y || ilat_loop.z )) && (it1==it2) && (ia1==ia2) )
									continue;
								const Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
								const Vector3<double> tau2 = ucell.atoms[it2].tau[ia2] + ilat_loop * ucell.latvec;
								const double r_sqr = (tau1 - tau2).norm2();
								const double r = sqrt(r_sqr);
								const double tmp_exp = exp( -para.damping* (r/R0_sum-1) );
								const double tmp_factor = C6_product/ pow(r_sqr,3)/ r/ (1+tmp_exp)* ( -6/r + tmp_exp/(1+tmp_exp)*para.damping/R0_sum);
								force[ucell.itia2iat(it1,ia1)] += tmp_factor*(tau1-tau2);
							} // end for ilat_loop
				} // end for ia2
			} // end for ia1
		} // end for it2
	} // end for it1
	for( int iat=0; iat!=ucell.nat; ++iat )
	{
		force[iat] *= para.scaling/ucell.lat0;
	}
}


void Vdwd2::cal_stress()
{
    TITLE("Vdwd2","stress");
	para.initset(ucell);

	stress.Zero();
	
	for( int it1=0; it1!=ucell.ntype; ++it1 )
	{
		for( int it2=0; it2!=ucell.ntype; ++it2 )
		{
			const double C6_product = sqrt( para.C6.at(ucell.atoms[it1].psd) * para.C6.at(ucell.atoms[it2].psd) )/pow(ucell.lat0,6);
			const double R0_sum = ( para.R0.at(ucell.atoms[it1].psd) + para.R0.at(ucell.atoms[it2].psd) )/ucell.lat0;
			if(!R0_sum)
				WARNING_QUIT("Input", "R0_sum can not be 0");
			for( int ia1=0; ia1!=ucell.atoms[it1].na; ++ia1 )
			{
				for( int ia2=0; ia2!=ucell.atoms[it2].na; ++ia2 )
				{
					Vector3<int> ilat_loop;
					for( ilat_loop.x = -para.period.x/2; ilat_loop.x <= (para.period.x-1)/2; ++ilat_loop.x )
						for( ilat_loop.y = -para.period.y/2; ilat_loop.y <= (para.period.y-1)/2; ++ilat_loop.y )
							for( ilat_loop.z = -para.period.z/2; ilat_loop.z <= (para.period.z-1)/2; ++ilat_loop.z )
							{
								if( (!( ilat_loop.x || ilat_loop.y || ilat_loop.z )) && (it1==it2) && (ia1==ia2) )
									continue;
								const Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
								const Vector3<double> tau2 = ucell.atoms[it2].tau[ia2] + ilat_loop * ucell.latvec;
								const Vector3<double> dr = tau2 - tau1;
								const double r_sqr = (tau1 - tau2).norm2();
								const double r = sqrt(r_sqr);
								const double tmp_exp = exp( -para.damping* (r/R0_sum-1) );
								const double tmp_factor = C6_product/ pow(r_sqr,3)/ r/ (1+tmp_exp)* ( -6/r + tmp_exp/(1+tmp_exp)*para.damping/R0_sum);
								stress += tmp_factor / 2 * Matrix3(
									dr.x*dr.x, dr.x*dr.y, dr.x*dr.z,
									dr.y*dr.x, dr.y*dr.y, dr.y*dr.z,
									dr.z*dr.x, dr.z*dr.y, dr.z*dr.z);
							} // end for ilat_loop
				} // end for ia2
			} // end for ia1
		} // end for it2
	} // end for it1
	stress *= para.scaling / ucell.omega;
}
