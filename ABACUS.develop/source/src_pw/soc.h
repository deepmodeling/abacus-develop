#ifndef SOC_H
#define SOC_H

#include "tools.h"

class Fcoef
{

public:

     Fcoef();
     ~Fcoef();

     int ind1,ind4,ind5;

     complex<double> *p;

     inline complex<double> &operator()(const int &i1,const int &i2,const int &i3,const int &i4,const int &i5)
     { return p[ind2*ind3*ind4*ind5*i1 + ind3*ind4*ind5*i2 + ind4*ind5*i3 + ind5*i4 + i5]; }

     inline const complex<double> &operator()(const int &i1,const int &i2,const int &i3,const int &i4,const int &i5)const
     { return p[ind2*ind3*ind4*ind5*i1 + ind3*ind4*ind5*i2 + ind4*ind5*i3 + ind5*i4 + i5]; }

     void create(const int i1, const int i2, const int i3);

//   void free();

private:

     int ind2,ind3;
};

//-----------------------
// spin-orbital coupling 
//-----------------------
class Soc
{

	public:

    Soc();
    ~Soc();

	void init();

	double spinor(const int l, const double j, const int m, const int spin);

	int sph_ind(const int l, const double j, const int m, const int spin);

	void rot_ylm(const int lmax);
	// complex<double> **rotylm;

	complex<double> *p_rot;

	complex<double> & rotylm(const int &i1,const int &i2)
	{ return p_rot[l_max*i1 + i2]; }

	Fcoef fcoef;

	Vector3<double> *m_loc;   //magnetization for each atom axis
	double *angle1;
	double *angle2;

	double ux[3];
	void cal_ux(const int ntype);
	bool lsign ;

	//int npol;
	private:

	bool judge_parallel(double a[3],Vector3<double> b);

	int l_max;
};
#endif
