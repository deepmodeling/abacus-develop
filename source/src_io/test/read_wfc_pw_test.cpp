#include"../wf_io.h"
#include"gtest/gtest.h"

/************************************************
*  unit test of read_wfc
***********************************************/

/**
 *  - Tested functions:
 *  	read_wfc()
 * 
 */

TEST(ReadWfcPWTest,txt)
{
	// input parameters
	std::string fn1 = "WAVEFUNC1.txt";
	std::string fn2 = "WAVEFUNC2.txt";
	int out_wfc_pw = 1;
	psi::Psi<std::complex<double>> psi;
	K_Vectors * p_kv = new K_Vectors;
	ModulePW::PW_Basis_K * wfc_basis = new ModulePW::PW_Basis_K;
	int ikstot = 1;
	int nkstot =2;
	int npw;
	int nbands;
	double ecut;
    	p_kv->kvec_c.resize(nkstot);
    	p_kv->wk.resize(nkstot);
	// call the function
	WF_io::read_wfc(fn1,out_wfc_pw,psi,p_kv,wfc_basis,ikstot,nbands,ecut,nkstot);

	int ik = ikstot-1;
	EXPECT_EQ(ikstot,1);
	EXPECT_EQ(nkstot,2);
	EXPECT_NEAR(p_kv->kvec_c[ik].x,0,1e-12);
	EXPECT_NEAR(p_kv->kvec_c[ik].y,0,1e-12);
	EXPECT_NEAR(p_kv->kvec_c[ik].z,0,1e-12);
	EXPECT_NEAR(p_kv->wk[ik],1,1e-12);
	EXPECT_EQ(wfc_basis->npw,10131);
	EXPECT_EQ(nbands,2);
	EXPECT_NEAR(ecut,20,1e-12);
	EXPECT_NEAR(wfc_basis->lat0,1.88973,1e-12);
	EXPECT_NEAR(wfc_basis->tpiba,3.32492,1e-12);
	EXPECT_NEAR(wfc_basis->getgcar(ik,wfc_basis->npw-1).z,1.3,1e-12);
	EXPECT_NEAR(psi(0,0,wfc_basis->npw-1).real(),-3.514392e-04,1e-12);
	EXPECT_NEAR(psi(0,0,wfc_basis->npw-1).imag(),-2.902740e-05,1e-12);
	EXPECT_NEAR(psi(0,1,wfc_basis->npw-1).real(),-6.761676e-05,1e-12);
	EXPECT_NEAR(psi(0,1,wfc_basis->npw-1).imag(),-5.420499e-06,1e-12);

	//
	ikstot = 2;

	WF_io::read_wfc(fn2,out_wfc_pw,psi,p_kv,wfc_basis,ikstot,nbands,ecut,nkstot);

	ik = ikstot-1;
	EXPECT_EQ(ikstot,2);
	EXPECT_EQ(nkstot,2);
	EXPECT_NEAR(p_kv->kvec_c[ik].x,0,1e-12);
	EXPECT_NEAR(p_kv->kvec_c[ik].y,0,1e-12);
	EXPECT_NEAR(p_kv->kvec_c[ik].z,0.05,1e-12);
	EXPECT_NEAR(p_kv->wk[ik],1,1e-12);
	EXPECT_EQ(wfc_basis->npw,10154);
	EXPECT_EQ(nbands,2);
	EXPECT_NEAR(ecut,20,1e-12);
	EXPECT_NEAR(wfc_basis->lat0,1.88973,1e-12);
	EXPECT_NEAR(wfc_basis->tpiba,3.32492,1e-12);
	EXPECT_NEAR(wfc_basis->getgcar(ik,wfc_basis->npw-1).z,1.2,1e-12);
	EXPECT_NEAR(psi(1,0,wfc_basis->npw-1).real(),9.199584e-05,1e-12);
	EXPECT_NEAR(psi(1,0,wfc_basis->npw-1).imag(),-4.052787e-04,1e-12);
	EXPECT_NEAR(psi(1,1,wfc_basis->npw-1).real(),5.699829e-06,1e-12);
	EXPECT_NEAR(psi(1,1,wfc_basis->npw-1).imag(),4.346193e-06,1e-12);

	// delete pointers
	delete p_kv;
	delete wfc_basis;
}

