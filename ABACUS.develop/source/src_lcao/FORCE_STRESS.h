#ifndef FORCE_STRESS_LCAO_H
#define FORCE_STRESS_LCAO_H

#include "../src_pw/tools.h"
#include "FORCE_k.h"
//#include "./force_lcao_gamma.h"
#include "../src_pw/stress_func.h"
#include "../input_conv.h"
#include "../src_pw/forces.h"

class Force_Stress_LCAO  
{
	// mohan add 2021-02-09
	friend class md;
	friend class Run_MD;
	friend void Input_Conv::Convert();
	friend class Update_input;
	friend class LOOP_ions;

	public :
	
	Force_Stress_LCAO ();
	~Force_Stress_LCAO ();

	private:

	Force_LCAO_k flk;
//	Force_LCAO_gamma flg;
	Stress_Func sc_pw;
	Forces f_pw;

	void allocate (void);
	void destroy (void);
	
	void print_force(const string &name, matrix& f, const bool screen, bool ry)const;
	void printforce_total (const bool ry, const bool istestf, matrix& fcs);
	
	void getForceStress(
		const bool isforce, 
		const bool isstress, 
		const bool istestf, 
		const bool istests, 
		matrix &fcs, 
		matrix &scs);

	void forceSymmetry(matrix &fcs);

	void calForcePwPart(
		matrix &fvl_dvl, 
		matrix &fewalds, 
		matrix &fcc, 
		matrix &fscc);

	void calForceStressIntegralPart(
		const bool isGammaOnly,
		const bool isforce,
		const bool isstress,
		matrix &foverlap,
		matrix &ftvnl_dphi,
		matrix &fvnl_dbeta,	
		matrix &fvl_dphi,
		matrix &soverlap,
		matrix &stvnl_dphi,
		matrix &svnl_dbeta,
		matrix &svl_dphi);

	void calStressPwPart(
		matrix &sigmadvl,
		matrix &sigmahar,
		matrix &sigmaewa,
		matrix &sigmacc,
		matrix &sigmaxc);
	
	static double force_invalid_threshold_ev;
	static double output_acc; // control the accuracy
};
#endif
