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
	int nbands;
	double ecut;
	int ikstot = 1;
	int nkstot =2;
	int ik = ikstot-1;
	psi::Psi<std::complex<double>> psi;
	K_Vectors * p_kv = new K_Vectors;
	ModulePW::PW_Basis_K * wfc_basis[nkstot];
	for (int ik=0;ik<nkstot;ik++)
	{
		wfc_basis[ik] = new ModulePW::PW_Basis_K;
	}
    	p_kv->kvec_c.resize(nkstot);
    	p_kv->wk.resize(nkstot);
	// call the function
	WF_io::read_wfc(fn1,out_wfc_pw,psi,p_kv,wfc_basis[ik],ikstot,nbands,ecut,nkstot);

	EXPECT_EQ(ikstot,1);
	EXPECT_EQ(nkstot,2);
	EXPECT_NEAR(p_kv->kvec_c[ik].x,0,1e-12);
	EXPECT_NEAR(p_kv->kvec_c[ik].y,0,1e-12);
	EXPECT_NEAR(p_kv->kvec_c[ik].z,0,1e-12);
	EXPECT_NEAR(p_kv->wk[ik],1,1e-12);
	EXPECT_EQ(wfc_basis[ik]->npw,10131);
	EXPECT_EQ(nbands,2);
	EXPECT_NEAR(ecut,20,1e-12);
	EXPECT_NEAR(wfc_basis[ik]->lat0,1.88973,1e-12);
	EXPECT_NEAR(wfc_basis[ik]->tpiba,3.32492,1e-12);
	EXPECT_NEAR(wfc_basis[ik]->gcar[wfc_basis[ik]->npw-1].z,1.3,1e-12);
	EXPECT_NEAR(psi(0,0,wfc_basis[ik]->npw-1).real(),-3.514392e-04,1e-12);
	EXPECT_NEAR(psi(0,0,wfc_basis[ik]->npw-1).imag(),-2.902740e-05,1e-12);
	EXPECT_NEAR(psi(0,1,wfc_basis[ik]->npw-1).real(),-6.761676e-05,1e-12);
	EXPECT_NEAR(psi(0,1,wfc_basis[ik]->npw-1).imag(),-5.420499e-06,1e-12);

	// read the second WAVEFUNC file
	ikstot = 2;
	ik = ikstot-1;
	WF_io::read_wfc(fn2,out_wfc_pw,psi,p_kv,wfc_basis[ik],ikstot,nbands,ecut,nkstot);

	EXPECT_EQ(ikstot,2);
	EXPECT_EQ(nkstot,2);
	EXPECT_NEAR(p_kv->kvec_c[ik].x,0,1e-12);
	EXPECT_NEAR(p_kv->kvec_c[ik].y,0,1e-12);
	EXPECT_NEAR(p_kv->kvec_c[ik].z,0.05,1e-12);
	EXPECT_NEAR(p_kv->wk[ik],1,1e-12);
	EXPECT_EQ(wfc_basis[ik]->npw,10154);
	EXPECT_EQ(nbands,2);
	EXPECT_NEAR(ecut,20,1e-12);
	EXPECT_NEAR(wfc_basis[ik]->lat0,1.88973,1e-12);
	EXPECT_NEAR(wfc_basis[ik]->tpiba,3.32492,1e-12);
	EXPECT_NEAR(wfc_basis[0]->gcar[wfc_basis[0]->npw-1].z,1.3,1e-12);
	EXPECT_NEAR(wfc_basis[ik]->gcar[wfc_basis[ik]->npw-1].z,1.2,1e-12);
	EXPECT_NEAR(psi(1,0,wfc_basis[ik]->npw-1).real(),9.199584e-05,1e-12);
	EXPECT_NEAR(psi(1,0,wfc_basis[ik]->npw-1).imag(),-4.052787e-04,1e-12);
	EXPECT_NEAR(psi(1,1,wfc_basis[ik]->npw-1).real(),5.699829e-06,1e-12);
	EXPECT_NEAR(psi(1,1,wfc_basis[ik]->npw-1).imag(),4.346193e-06,1e-12);

	// delete pointers
	delete p_kv;
	for (int ik=0;ik<nkstot;ik++)
	{
		delete wfc_basis[ik];
	}
}

