#include "conv_coulomb_pot.h"
#include "conv_coulomb_pot-inl.h"
#include "src_global/mathzone.h"
#include "src_global/constants.h"
#include "src_global/global_function.h"

#include "src_external/src_test/test_function.h"

Conv_Coulomb_Pot::Conv_Coulomb_Pot(const Numerical_Orbital_Lm &orb_in)
	:orb(orb_in)
{
	conv_coulomb_pot.resize(orb.getNr());
}


double Conv_Coulomb_Pot::get_conv_coulomb_pot(const size_t &ir) const
{
	if( ir < conv_coulomb_pot.size() )
	{
		return conv_coulomb_pot[ir];
	}
	else
	{
		double r_radial = get_r_radial_outofrange( ir );
		return conv_coulomb_pot_extra / pow(r_radial,orb.getL()+1);		
	}
}


void Conv_Coulomb_Pot::cal_conv_coulomb_pot()
{
	vector<double> tmp_func( orb.getNr() );
	vector<double> tmp_integral( orb.getNr() );

	// \int_{r'=0}^{r} dr' f(r') * r'^{L+2} / r^{L+1}
	for( size_t ir=0; ir!=orb.getNr(); ++ir )
	{
		tmp_func[ir] = orb.getPsi(ir) * pow( orb.getRadial(ir), orb.getL()+2 );
	}	
	Mathzone::Simpson_Integral_0toall( orb.getNr(), VECTOR_TO_PTR(tmp_func), orb.getRab(), VECTOR_TO_PTR(tmp_integral) );
	conv_coulomb_pot[0]=0;
	for( size_t ir=1; ir!=orb.getNr(); ++ir )
	{
		conv_coulomb_pot[ir] = tmp_integral[ir] / pow( orb.getRadial(ir), orb.getL()+1 );
	}

	// \int_{r'=0}^{Rcut} dr' f(r') * r'^{L+2}
	conv_coulomb_pot_extra = tmp_integral.back();

	// \int_{r'=r}^{+\infty} dr' f(r') * r^{L} / r'^{L-1}
	tmp_func[0]=0;
	for( size_t ir=1; ir!=orb.getNr(); ++ir )
	{
		tmp_func[ir] = orb.getPsi(ir) / pow( orb.getRadial(ir), orb.getL()-1 );
	}
	Mathzone::Simpson_Integral_alltoinf( orb.getNr(), VECTOR_TO_PTR(tmp_func), orb.getRab(), VECTOR_TO_PTR(tmp_integral) );
	for( size_t ir=0; ir!=orb.getNr(); ++ir )
	{
		conv_coulomb_pot[ir] += tmp_integral[ir] * pow( orb.getRadial(ir), orb.getL() );
	}

	// 4pi/(2L+1)
	const double coefficient = 4*PI/(2*orb.getL()+1);
	for( size_t ir=0; ir!=orb.getNr(); ++ir )
	{
		conv_coulomb_pot[ir] *= coefficient;
	}
	conv_coulomb_pot_extra *= coefficient;
}


template<>
void Conv_Coulomb_Pot::cal_orbs_ccp<Numerical_Orbital_Lm>( 
	const Numerical_Orbital_Lm & orbs, 
	Numerical_Orbital_Lm & orbs_ccp, 
	size_t rmesh_times, 
	size_t kmesh_times )
{
	Conv_Coulomb_Pot ccp( orbs );
	ccp.cal_conv_coulomb_pot();

	const size_t mesh = (rmesh_times * orbs.getNr()) | 1;		// mesh must be odd for simspon integral
	vector<double> psi( mesh );
	for( size_t ir=0; ir!=mesh ; ++ir)
	{
		psi[ir] = ccp.get_conv_coulomb_pot(ir);
	}

	vector<double> rab( mesh );
	for( size_t ir=0; ir!=mesh ; ++ir)
	{
		if( ir<orbs.getNr() )
			rab[ir] = orbs.getRab(ir);
		else
			rab[ir] = orbs.getRab(orbs.getNr()-1);
	}

	vector<double> r_radial( mesh );
	for( size_t ir=0; ir!=mesh ; ++ir)
	{
		if( ir<orbs.getNr() )
			r_radial[ir] = orbs.getRadial(ir);
		else
			r_radial[ir] = ccp.get_r_radial_outofrange(ir);
	}
	
	orbs_ccp.set_orbital_info(
 		orbs.getLabel(),
	 	orbs.getType(),
		orbs.getL(),
		orbs.getChi(),
	    mesh,
		VECTOR_TO_PTR(rab),
		VECTOR_TO_PTR(r_radial),
		Numerical_Orbital_Lm::Psi_Type::Psi,
		VECTOR_TO_PTR(psi),
		orbs.getNk() * kmesh_times | 1,								// Nk must be odd
		orbs.getDk(),												// Peize Lin change 2017-04-16
//		orbs.getDk() / kmesh_times,
		orbs.getDruniform(),
		false,
		true);
}
