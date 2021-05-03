//==========================================================
// Author: Xin Qu
// DATE : 2019-12-10
//==========================================================
#ifndef DFTU_RELAX_H
#define DFTU_RELAX_H

#include <vector>
#include "dftu_yukawa.h"
#include "../src_global/complexmatrix.h"

using namespace std;

//==========================================================
// CLASS :
// NAME : DFTU_RELAX
//==========================================================
class DFTU_RELAX : public DFTU_Yukawa
{

public:

    DFTU_RELAX();
    ~DFTU_RELAX();

    void force_stress();
    void folding_dSm_soverlap();
    void allocate_force_stress();
    void erase_force_stress();
    void cal_force_k(const vector<vector<complex<double>>> &VU);
    void cal_force_gamma(const vector<vector<double>> &VU);
    void cal_stress_k(const vector<vector<complex<double>>> &VU);
    void cal_stress_gamma(const vector<vector<double>> &VU);

    double get_onebody_eff_pot
    (
        const int T, const int iat,
	    const int L, const int N, const int spin, 
	    const int m0,
        const int m1,
        const int type, const bool newlocale
    );

    //forces and stress
    vector<vector<double>> force_dftu;      //force_dftu[iat][dim] 
    vector<vector<double>> stress_dftu;

    //vector<vector<vector<complex<double>>>> dSm_k;            //dSm_k[ik][dim][irc]
    complex<double> ***dSm_k;                                   //dSm_k[ik][dim][irc]
    //vector<vector<vector<complex<double>>>> soverlap_k;       //soverlap_k[ik][xy][irc]
    complex<double> ***soverlap_k;                              //soverlap_k[ik][xy][irc]
    //vector<vector<double>> soverlap_gamma;                    //soverlap_gamma[xy][irc]
    double **soverlap_gamma;                                    //soverlap_gamma[xy][irc]

    //transform between iwt index and it, ia, L, N and m index
    vector<vector<vector<vector<vector<int>>>>> iatlnmipol2iwt;   //iatlnm2iwt[iat][l][n][m][ipol]
    vector<int> iwt2it;                               //iwt2it[iwt]
    vector<int> iwt2iat;                              //iwt2iat[iwt]
    vector<int> iwt2l;                                //iwt2l[iwt]
    vector<int> iwt2n;                                //iwt2n[iwt]
    vector<int> iwt2m;                                //iwt2m[iwt]
    vector<int> iwt2ipol;                             //iwt2ipol[iwt]
    vector<int> iat2it; 

    //local occupancy matrix of the correlated subspace
    //locale: the out put local occupation number matrix of correlated electrons in the current electronic step
    //locale_save: the input local occupation number matrix of correlated electrons in the current electronic step
    vector<vector<vector<vector<matrix>>>> locale;            // locale[iat][l][n][spin](m1,m2)
    vector<vector<vector<vector<matrix>>>> locale_save;       // locale_save[iat][l][n][spin](m1,m2)

};

#endif
