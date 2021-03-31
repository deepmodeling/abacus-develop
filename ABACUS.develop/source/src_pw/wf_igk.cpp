#include "global.h"
#include "wf_igk.h"

WF_igk::WF_igk()
{
    g2kin = new double[1];
}

WF_igk::~WF_igk()
{ 
	if(test_deconstructor)
	{
		cout << " ~WF_igk()" << endl;
	}
    delete[] g2kin;
}

//=======================================================
//  find out all G vectors that |k+G|^2<=GPsi2 for each
//  k point, and create indices ind_Gk(ngk,nk)
//  between the G vectors in the
//  k point array and the 1d Gvector array G1d[:].
//  set npwx
//  set igk
//========================================================
int WF_igk::setupIndGk(const PW_Basis &pwb,const int nks)
{
    TITLE("WF_igk","setupIndGk");
    timer::tick("WF_igk","setupIndGk");

    //============================================
    // Find out for each k point,
    // how many planewave within the radius ggpsi
    // and the index ind_Gk array
    //============================================
    int npw_max = 0;

    //=============================================
    // Done this in each cpu ( use nks,not nkstot )
    // notice : using cartesian coordinate
    //=============================================
    for (int ik = 0; ik < nks; ik++)
    {
        int ng = 0;
        for (int ig = 0; ig < pwb.ngmw ; ig++)
        {
            Vector3<double> f =  pwb.gcar[ig] + kv.kvec_c[ik];
            const double gk2 = f * f;
            const double k2 = kv.kvec_c[ik] * kv.kvec_c[ik];
            if (sqrt(pwb.gg[ig]) > sqrt(pwb.ggpsi) + sqrt(k2))
            {
                break;
            }
            if (gk2 <= pwb.ggpsi)
            {
                ng++;
            }
        }
        kv.ngk[ik] = ng;
        if ( npw_max < ng)
        {
            npw_max = ng;
        }
//		ofs_running << " " << setw(8) << ik << setw(10) << kv.ngk[ik] << endl;
    }

    if (test_wf > 1) OUT("npw_max",npw_max);
//----------------------------------------------------------
// EXPLAIN : correspondence K + G <- -> G
//----------------------------------------------------------
    igk.create(nks, npw_max);
    Memory::record("WF_igk","igk",nks*npw_max,"int");

//----------------------------------------------------------
// EXPLAIN : Calculate again ! (Not a smart job)
//----------------------------------------------------------
    for (int ik = 0; ik < nks; ik++)
    {
        int ng = 0;
        for (int ig = 0; ig < pwb.ngmw; ig++)
        {
            Vector3<double> f =  pwb.gcar[ig] + kv.kvec_c[ik];
            const double gk2 = f * f;
            const double k2 = kv.kvec_c[ik] * kv.kvec_c[ik];
            if (sqrt(pwb.gg[ig]) > sqrt(pwb.ggpsi) + sqrt(k2))
            {
                break;
            }
            if (gk2 <= pwb.ggpsi)
            {
                this->igk(ik, ng) = ig;
                ng++;
            }
        }
    }

	bool out_gk =0; //DIY! mohan 2011-10-03
	if(out_gk)
	{
		stringstream ss;
		ss << global_out_dir << "PW_GK" << MY_RANK+1 << ".dat";
		ofstream ofs( ss.str().c_str() );
		ofs << pw.ggpsi << " (ggpsi, Ry)" << endl;
		ofs << pw.ggwfc << " (ggwfc, Ry)" << endl;
		ofs << pw.ggwfc2 << " (ggwfc2, Ry)" << endl;
		Vector3<double> f;
		for(int ik=0; ik < nks; ++ik)
		{
			ofs << ik+1 << " (Index of k)" << endl;
			ofs << pwb.ngmw << " (Number of plane waves)" << endl;	
			for(int ig=0; ig < pwb.ngmw; ++ig)
			{
				f =  pwb.gcar[ig] + kv.kvec_c[ik];
				ofs << f.x << " " << f.y << " " << f.z << " " << f.norm() << endl;
			}
		}
		ofs.close();
	}

    timer::tick("WF_igk","setupIndGk");
    return npw_max;
} // end setupIndGk()


//--------------------------------------------------------
// Compute kinetic energy for each k-point
//--------------------------------------------------------
void WF_igk::ekin(const int ik)
{
    timer::tick("WF_igk","ekin");
    zeros( g2kin, this->npwx);

    for (int ig = 0;ig < kv.ngk[ik];ig++)
    {
        Vector3<double> v3 = kv.kvec_c[ik] + pw.gcar[ this->igk(ik, ig) ];
        this->g2kin[ig] = v3 * v3;
//--------------------------------------------------------
// EXPLAIN : Put the correct units on the kinetic energy
//--------------------------------------------------------
        this->g2kin[ig] = g2kin[ig] * ucell.tpiba2;

    }
    timer::tick("WF_igk","ekin");
    return ;
}


Vector3<double> WF_igk::get_1qvec_cartesian(const int ik,const int ig)const
{
    Vector3<double> qvec = kv.kvec_c[ik] + pw.gcar[ this->igk(ik, ig) ];

	/*
	if(igk(ik,ig)==0)
	{
		cout << " g add = " << pw.gcar << endl;
		cout << "\n igk = " << igk(ik,ig);
		cout << " k=" << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z;
		cout << " g=" << pw.gcar[ this->igk(ik, ig) ].x << " " << pw.gcar[ this->igk(ik, ig) ].y << " " << pw.gcar[ this->igk(ik, ig) ].z;
	}
	*/
	
    return qvec;
}


double* WF_igk::get_qvec_cartesian(const int &ik)
{
    double *qmod = new double[ kv.ngk[ik] ];
    for (int ig=0; ig< kv.ngk[ik]; ig++)
    {
        // cartesian coordinate
        Vector3<double> qvec = kv.kvec_c[ik] + pw.gcar[ this->igk(ik, ig) ];
        // modulus, in ucell.tpiba unit.
        const double q2 = qvec*qvec;
        qmod[ig] = ucell.tpiba * sqrt(q2);//sqrt(q2) * ucell.tpiba;
    }
    return qmod;
}


complex<double>* WF_igk::get_sk(const int ik, const int it, const int ia)const
{
    timer::tick("WF_igk","get_sk");
    const double arg = (kv.kvec_c[ik] * ucell.atoms[it].tau[ia]) * TWO_PI;
    const complex<double> kphase = complex <double> ( cos(arg),  -sin(arg) );
    complex<double> *sk = new complex<double>[ kv.ngk[ik] ];
    for (int ig=0; ig<kv.ngk[ik]; ig++)
    {
        const int iig = this->igk(ik, ig);
        const int iat = ucell.itia2iat(it, ia);
        sk[ig] = kphase
                 * pw.eigts1(iat, pw.ig1[iig])
                 * pw.eigts2(iat, pw.ig2[iig])
                 * pw.eigts3(iat, pw.ig3[iig]);
    }
    timer::tick("WF_igk","get_sk");
    return sk;
}


complex<double>* WF_igk::get_skq(int ik, const int it, const int ia, Vector3<double> q)   //pengfei 2016-11-23 
{
    complex<double> *skq = new complex<double>[ kv.ngk[ik] ];

    for (int ig=0; ig<kv.ngk[ik]; ig++)
    {
         Vector3<double> qkq = kv.kvec_c[ik] + pw.gcar[this->igk(ik, ig)] + q;
         double arg = (qkq * ucell.atoms[it].tau[ia]) * TWO_PI;
         skq[ig] = complex <double> ( cos(arg),  -sin(arg) );
    }

    return skq;
}

