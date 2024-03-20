//=========================================================
//AUTHOR : Peize Lin
//DATE : 2016-01-24
//=========================================================

#include "center2_orb-orb11.h"

#include "module_base/constants.h"
#include "module_base/ylm.h"
#include "module_base/math_polyint.h"

#include <cmath>

Center2_Orb::Orb11::Orb11(
	const Numerical_Orbital_Lm &nA_in,
	const Numerical_Orbital_Lm &nB_in,
	const ORB_table_phi &MOT_in,
	const ORB_gaunt_table &MGT_in
)
	:nA(nA_in),
	 nB(nB_in),
	 MOT(MOT_in),
	 MGT(MGT_in)
{}

void Center2_Orb::Orb11::init_radial_table(void)
{
	const int LA = this->nA.getL();
	const int LB = this->nB.getL();
	for( int LAB = std::abs(LA-LB); LAB<=LA+LB; ++LAB)
	{
		if( (LAB-std::abs(LA-LB))%2==1 )			// if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
		{
			continue;
		}

		const int rmesh = this->MOT.get_rmesh(this->nA.getRcut(),this->nB.getRcut()) ;

		this->Table_r[LAB].resize(rmesh,0);
		this->Table_dr[LAB].resize(rmesh,0);

		this->MOT.cal_ST_Phi12_R(
			1,
			LAB,
			this->nA,
			this->nB,
			rmesh,
			this->Table_r[LAB].data(),
			this->Table_dr[LAB].data());
	}
	return;
}

void Center2_Orb::Orb11::init_radial_table( const std::set<size_t> &radials )
{
	const int LA = this->nA.getL();
	const int LB = this->nB.getL();

	const size_t rmesh = this->MOT.get_rmesh(this->nA.getRcut(),this->nB.getRcut());

	std::set<size_t> radials_used;
	for( const size_t &ir : radials )
		if( ir<rmesh )
			radials_used.insert(ir);

	for( int LAB = std::abs(LA-LB); LAB<=LA+LB; ++LAB)
	{
		if( (LAB-std::abs(LA-LB))%2==1 )			// if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
			continue;

		this->Table_r[LAB].resize(rmesh,0);
		this->Table_dr[LAB].resize(rmesh,0);

		this->MOT.cal_ST_Phi12_R(
			1,
			LAB,
			this->nA,
			this->nB,
			radials_used,
			this->Table_r[LAB].data(),
			this->Table_dr[LAB].data());
	}
}

double Center2_Orb::Orb11::cal_overlap(
	const ModuleBase::Vector3<double> &RA, const ModuleBase::Vector3<double> &RB,
	const int &mA, const int &mB) const
{
	const double tiny1 = 1e-12;		// why 1e-12?
	const double tiny2 = 1e-10;		// why 1e-10?

	const ModuleBase::Vector3<double> delta_R = RB-RA;
	const double distance_true = delta_R.norm();
	const double distance = (distance_true>=tiny1) ? distance_true : distance_true+tiny1;
	const double RcutA = this->nA.getRcut();
	const double RcutB = this->nB.getRcut();
	if( distance > (RcutA + RcutB) ) return 0.0;

	const int LA = this->nA.getL();
	const int LB = this->nB.getL();

	std::vector<double> rly;
	ModuleBase::Ylm::rl_sph_harm (
		LA+LB,											// max LAB
		delta_R.x, delta_R.y, delta_R.z,
		rly);

	double overlap = 0.0;

	for( const auto &tb_r : this->Table_r )
	{
		const int LAB = tb_r.first;

		for( int mAB=0; mAB!=2*LAB+1; ++mAB )
		// const int mAB = mA + mB;
		{
			const double Gaunt_real_A_B_AB =
				this->MGT.Gaunt_Coefficients (
					this->MGT.get_lm_index(LA,mA),
					this->MGT.get_lm_index(LB,mB),
					this->MGT.get_lm_index(LAB,mAB));
			if( 0==Gaunt_real_A_B_AB )	continue;

			const double ylm_solid = rly[ this->MGT.get_lm_index(LAB, mAB) ];
			if( 0==ylm_solid ) continue;
			const double ylm_real =
				(distance > tiny2) ?
				ylm_solid / pow(distance,LAB) :
				ylm_solid;

			const double i_exp = std::pow(-1.0, (LA - LB - LAB) / 2);

			const double Interp_Tlm =
				(distance > tiny2)
				? ModuleBase::PolyInt::Polynomial_Interpolation(
					tb_r.second.data(),
					this->MOT.get_rmesh(RcutA, RcutB),
					this->MOT.dr,
					distance )
				: tb_r.second.at(0);

			overlap +=
	//			pow(2*PI,1.5)
				i_exp
				* Gaunt_real_A_B_AB
				* Interp_Tlm
				* ylm_real;

		}
	}

	return overlap;
}

ModuleBase::Vector3<double> Center2_Orb::Orb11::cal_grad_overlap(    //caoyu add 2021-11-19
	const ModuleBase::Vector3<double> &RA, const ModuleBase::Vector3<double> &RB,
	const int& mA, const int& mB) const
{
	const double tiny1 = 1e-12;		// same as `cal_overlap`
	const double tiny2 = 1e-10;		//same as `cal_overlap`

	const ModuleBase::Vector3<double> delta_R = RB-RA;
	const double distance_true = delta_R.norm();
	const double distance = (distance_true>=tiny1) ? distance_true : distance_true+tiny1;
	const double RcutA = this->nA.getRcut();
	const double RcutB = this->nB.getRcut();
	if( distance > (RcutA + RcutB) )
		return ModuleBase::Vector3<double>(0.0, 0.0, 0.0);

	const int LA = this->nA.getL();
	const int LB = this->nB.getL();

	std::vector<double> rly;
    std::vector<vector<double>> tmp_grly;
    std::vector<ModuleBase::Vector3<double>> grly;
    ModuleBase::Ylm::grad_rl_sph_harm(
        LA + LB,
        delta_R.x, delta_R.y, delta_R.z,
        rly, tmp_grly);
    for (const auto& tmp_ele : tmp_grly)
    {
        ModuleBase::Vector3<double> ele(tmp_ele[0], tmp_ele[1], tmp_ele[2]);
        grly.push_back(ele);
    }

    ModuleBase::Vector3<double> grad_overlap(0.0, 0.0, 0.0);

    for (const auto& tb_r : this->Table_r)
	{
        const int LAB = tb_r.first;
		for( int mAB=0; mAB!=2*LAB+1; ++mAB )
		// const int mAB = mA + mB;
		{
			const double Gaunt_real_A_B_AB =
				this->MGT.Gaunt_Coefficients (
					this->MGT.get_lm_index(LA,mA),
					this->MGT.get_lm_index(LB,mB),
					this->MGT.get_lm_index(LAB,mAB));
			if( 0==Gaunt_real_A_B_AB )	continue;

			const double ylm_solid = rly[ this->MGT.get_lm_index(LAB, mAB) ];
			const double ylm_real =
				(distance > tiny2) ?
				ylm_solid / pow(distance,LAB) :
                ylm_solid;

            const ModuleBase::Vector3<double> gylm_solid = grly[ this->MGT.get_lm_index(LAB, mAB) ];
			const ModuleBase::Vector3<double> gylm_real =
				(distance > tiny2)
				? gylm_solid / pow(distance,LAB)
                : gylm_solid;

            const double i_exp = std::pow(-1.0, (LA - LB - LAB) / 2);

			const double Interp_Tlm =
				(distance > tiny2)
				? ModuleBase::PolyInt::Polynomial_Interpolation(
					tb_r.second.data(),
					this->MOT.get_rmesh(RcutA, RcutB),
					this->MOT.dr,
					distance )
                : tb_r.second.at(0);

            const double grad_Interp_Tlm =
                (distance > tiny2)
                ? ModuleBase::PolyInt::Polynomial_Interpolation(
                    this->Table_dr.at(LAB).data(),
                    this->MOT.get_rmesh(RcutA, RcutB),
                    this->MOT.dr,
                    distance)   //Interp(Table_dr)
                  - Interp_Tlm * LAB / distance
				: 0.0;

            grad_overlap +=
                i_exp   //			pow(2*PI,1.5)
                * Gaunt_real_A_B_AB
                * (Interp_Tlm * gylm_real
                + grad_Interp_Tlm * ylm_real *delta_R / distance);
		}
	}

	return grad_overlap;
}