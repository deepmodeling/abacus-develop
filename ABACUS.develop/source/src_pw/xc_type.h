//==========================================================
// AUTHOR : Lixin He, mohan
// DATE : 2008-11-08
//==========================================================
#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

#include "tools.h"
class xcfunc
{

	public:

	//module funct
	// char dft[20];
	// dft is the exchange-correlation functional, described by
	// any nonconflicting combination of the following keywords
	// (case-insensitive):

	// Exchange:    "nox"    none                           iexch=0
	//              "sla"    Slater (alpha=2/3)             iexch=1 (default)
	//              "sl1"    Slater (alpha=1.0)             iexch=2
	//              "rxc"    Relativistic Slater            iexch=3
	//				"oep"    Optimized Effective Potential  iexch=4
	//				"hf"     Hartree-Fock                   iexch=5
	//				"pb0x"   PBE0                           iexch=6
	//				"b3lpx"  B3LYP                          iexch=7
	//				"kzk"    Finite-size corrections        iexch=8
	//				"hsex"   HSE06                          iexch=9
	//
	// Correlation: "noc"    none                           icorr=0
	//              "pz"     Perdew-Zunger                  icorr=1 (default)
	//              "vwn"    Vosko-Wilk-Nusair              icorr=2
	//              "lyp"    Lee-Yang-Parr                  icorr=3
	//              "pw"     Perdew-Wang                    icorr=4
	//              "wig"    Wigner                         icorr=5
	//              "hl"     Hedin-Lunqvist                 icorr=6
	//              "obz"    Ortiz-Ballone form for PZ      icorr=7
	//              "obw"    Ortiz-Ballone form for PW      icorr=8
	//              "gl"     Gunnarson-Lunqvist             icorr=9
	//
	// Gradient Correction on Exchange:
	//              "nogx"   none                           igcx =0 (default)
	//              "b88"    Becke88 (beta=0.0042)          igcx =1
	//              "ggx"    Perdew-Wang 91                 igcx =2
	//              "pbx"    Perdew-Burke-Ernzenhof exch    igcx =3
	//              "rpb"    revised PBE by Zhang-Yang      igcx =4
	//              "hcth"   Cambridge exch, Handy et al    igcx =5
	//              "optx"   Handy's exchange functional    igcx =6
	//				"meta"   meta-gga                       igcx =7
	//				"pb0x"   PBE0                           igcx =8
	//				"b3lp"   B3LYP                          igcx =9
	//				"psx"    PBEsol exchange                igcx =10
	//				"wcx"    Wu-Cohen                       igcx =11
	//				"hsex"   HSE06                          igcx =12
	//
	// Gradient Correction on Correlation:
	//              "nogc"   none                           igcc =0 (default)
	//              "p86"    Perdew86                       igcc =1
	//              "ggc"    Perdew-Wang 91 corr.           igcc =2
	//              "blyp"   Lee-Yang-Parr                  igcc =3
	//              "pbc"    Perdew-Burke-Ernzenhof corr    igcc =4
	//              "hcth"   Cambridge corr, Handy et al    igcc =5
	//				"meta"   meta-gga                       igcc =6				
	//
	// Special cases:
	//              "bp"    = "b88+p86"           = Becke-Perdew grad.corr.
	//              "pw91"  = "pw +ggx+ggc"       = PW91 (aka GGA)
	//              "blyp"  = "sla+b88+lyp+blyp"  = BLYP
	//              "pbe"   = "sla+pw+pbx+pbc"    = PBE
	//              "revpbe"= "sla+pw+rpb+pbc"    = revPBE (Zhang-Yang)
	//              "pbesol"= "sla+pw+psx+psc"    = PBEsol
	//              "hcth"  = "nox+noc+hcth+hcth" = HCTH/120
	//              "olyp"  = "nox+lyp+optx+blyp"!!! UNTESTED !!!
	//              "tpss"  = "sla+pw+meta+meta"  = TPSS Meta-GGA
	//              "wc"    = "sla+pw+wcx+pbc"    = Wu-Cohen
	//              "pbe0"  = "pb0x+pw+pb0x+pbc"  = PBE0
	//              "b3lyp" = "b3lp+vwn+b3lp+b3lp"= B3LYP
	//              "hse"   = "hsex+pw+hsex+pbc"  = HSE	
	
	
	// References:
	//              pz      J.P.Perdew and A.Zunger, PRB 23, 5048 (1981)
	//              vwn     S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
	//				wig		E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)
	//				hl		L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)
	//				gl		O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976)
	//              pw      J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)
	//              obpz    G.Ortiz and P.Ballone, PRB 50, 1391 (1994)
	//              obpw    as above
	//              b88     A.D.Becke, PRA 38, 3098 (1988)
	//              p86     J.P.Perdew, PRB 33, 8822 (1986)
	//              pbe     J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
	//              pw91    J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)
	//              blyp    C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)
	//              hcth    Handy et al, JCP 109, 6264 (1998)
	//              olyp    Handy et al, JCP 116, 5411 (2002)
	//              revPBE  Zhang and Yang, PRL 80, 890 (1998)
	//				meta    J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria,PRL 91, 146401 (2003)
	//				kzk     H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)
	//				pbe0    J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)
	//				b3lyp   P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch J.Phys.Chem 98, 11623 (1994)
	//				pbesol	J.P. Perdew et al., PRL 100, 136406 (2008)
	//				wc		Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
	//              hse     Paier J, Marsman M, Hummer K, et al, JPC 124(15): 154709 (2006)

	int iexch; 
	int icorr;
	int igcx;
	int igcc;
	// internal indices for exchange-correlation
	//    iexch: type of exchange
	//    icorr: type of correlation
	//    igcx:  type of gradient correction on exchange
	//    igcc:  type of gradient correction on correlations
	// see comments above and routine "which_dft" below
	
	// Peize Lin add 2016-12-03
	int iexch_now; 
	int icorr_now;
	int igcx_now;
	int igcc_now;

	xcfunc();
	~xcfunc();

	// there are four values of dft.
	void which_dft(const string *dft);
	void printdft(ofstream &ofs);
	void ostreamdft(ostream &ofs);  // zws add 20150108
private:
	void set_dft_value(int &m,const int i);
	bool match_one(const string* dft, const string &name)const; 
	void copy_to_now();	// Peize Lin add 2016-12-03
	void hybrid_first();	// Peize Lin add 2016-12-03
};

#endif //FUNCTION_H
