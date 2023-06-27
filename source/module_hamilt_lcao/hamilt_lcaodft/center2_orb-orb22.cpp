//=========================================================
//AUTHOR : Peize Lin
//DATE : 2016-01-24
//=========================================================

#include "center2_orb-orb22.h"

Center2_Orb::Orb22::Orb22(
	const Numerical_Orbital_Lm &nA1_in,
	const Numerical_Orbital_Lm &nA2_in,
	const Numerical_Orbital_Lm &nB1_in,
	const Numerical_Orbital_Lm &nB2_in,
	const ORB_table_phi &MOT_in,
	const ORB_gaunt_table &MGT_in
)
	:nA1(nA1_in),
	 nA2(nA2_in),
	 nB1(nB1_in),
	 nB2(nB2_in),
	 MOT(MOT_in),
	 MGT(MGT_in){}

void Center2_Orb::Orb22::init_radial_table()
{
	const Numerical_Orbital_Lm & nB_short = (nB1.getNr()<=nB2.getNr()) ? nB1 : nB2;

	std::vector<double> nB_tmp(nB_short.getNr());
	for( size_t ir=0; ir!=nB_tmp.size(); ++ir)
	{
		nB_tmp[ir] = nB1.getPsi(ir) * nB2.getPsi(ir);
	}

	const int LB1 = nB1.getL();
	const int LB2 = nB2.getL();
	for( int LB = abs(LB1-LB2); LB<=LB1+LB2; ++LB)
	{
		if( (LB-std::abs(LB1-LB2))%2==1 )			// if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
			continue;

		this->nB[LB].set_orbital_info(
			nB_short.getLabel(),
			nB_short.getType(),
			LB,
			1,						// N?
			nB_short.getNr(),
			nB_short.getRab(),
			nB_short.getRadial(),
			Numerical_Orbital_Lm::Psi_Type::Psi,
			nB_tmp.data(),
			nB_short.getNk(),
			nB_short.getDk(),
			nB_short.getDruniform(),
			false,
			true,
			GlobalV::CAL_FORCE); // mohan add 2021-05-07

		this->orb21s.insert( std::make_pair( LB, Center2_Orb::Orb21( nA1, nA2, this->nB[LB], this->MOT, this->MGT ) ) );

		this->orb21s.at(LB).init_radial_table();
	}
}

void Center2_Orb::Orb22::init_radial_table( const std::set<size_t> &radials )
{
	const Numerical_Orbital_Lm & nB_short = (nB1.getNr()<=nB2.getNr()) ? nB1 : nB2;

	std::vector<double> nB_tmp(nB_short.getNr());
	for( size_t ir=0; ir!=nB_tmp.size(); ++ir)
	{
		nB_tmp[ir] = nB1.getPsi(ir) * nB2.getPsi(ir);
	}

	const int LB1 = nB1.getL();
	const int LB2 = nB2.getL();
	for( int LB = abs(LB1-LB2); LB<=LB1+LB2; ++LB)
	{
		if( (LB-std::abs(LB1-LB2))%2==1 )			// if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
			continue;

		this->nB[LB].set_orbital_info(
			nB_short.getLabel(),
			nB_short.getType(),
			LB,
			1,						// N?
			nB_short.getNr(),
			nB_short.getRab(),
			nB_short.getRadial(),
			Numerical_Orbital_Lm::Psi_Type::Psi,
			nB_tmp.data(),
			nB_short.getNk(),
			nB_short.getDk(),
			nB_short.getDruniform(),
			false,
			true, GlobalV::CAL_FORCE);

		this->orb21s.insert( std::make_pair( LB, Center2_Orb::Orb21( nA1, nA2, this->nB[LB], this->MOT, this->MGT ) ) );

		this->orb21s.at(LB).init_radial_table(radials);
	}
}

double Center2_Orb::Orb22::cal_overlap(
	const ModuleBase::Vector3<double> &RA, const ModuleBase::Vector3<double> &RB,
	const int &mA1, const int &mA2, const int &mB1, const int &mB2) const
{
	const int LB1 = nB1.getL();
	const int LB2 = nB2.getL();

	double overlap = 0.0;

	for( const auto& orb21 : this->orb21s )
	{
		const int LB = orb21.first;

		for( int mB=0; mB!=2*LB+1; ++mB )
		// const int mB=mB1+mB2;
		{
			const double Gaunt_real_B1_B2_B12 =
				this->MGT.Gaunt_Coefficients (
					this->MGT.get_lm_index(LB1,mB1),
					this->MGT.get_lm_index(LB2,mB2),
					this->MGT.get_lm_index(LB,mB));
			if( 0==Gaunt_real_B1_B2_B12 )	continue;

			overlap += Gaunt_real_B1_B2_B12 * orb21.second.cal_overlap(RA, RB, mA1, mA2, mB);
		}
	}

	return overlap;
}

ModuleBase::Vector3<double> Center2_Orb::Orb22::cal_grad_overlap(
	const ModuleBase::Vector3<double> &RA, const ModuleBase::Vector3<double> &RB,
	const int &mA1, const int &mA2, const int &mB1, const int &mB2) const
{
	const int LB1 = nB1.getL();
	const int LB2 = nB2.getL();

	ModuleBase::Vector3<double> grad_overlap(0.0, 0.0, 0.0);

	for( const auto& orb21 : this->orb21s )
	{
		const int LB = orb21.first;

		for( int mB=0; mB!=2*LB+1; ++mB )
		// const int mB=mB1+mB2;
		{
			const double Gaunt_real_B1_B2_B12 =
				this->MGT.Gaunt_Coefficients (
					this->MGT.get_lm_index(LB1,mB1),
					this->MGT.get_lm_index(LB2,mB2),
					this->MGT.get_lm_index(LB,mB));
			if( 0==Gaunt_real_B1_B2_B12 )	continue;

			grad_overlap += Gaunt_real_B1_B2_B12 * orb21.second.cal_grad_overlap(RA, RB, mA1, mA2, mB);
		}
	}

	return grad_overlap;
}
