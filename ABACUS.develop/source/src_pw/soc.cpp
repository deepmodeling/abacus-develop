#include"soc.h"

Fcoef::Fcoef()
{
     ind2 = 2;
     ind3 = 2;
     ind1 = 1;
     ind4 = 1;
     ind5 = 1;
     p = new complex<double> [1];
}

Fcoef::~Fcoef()
{
     delete[] p;
}

void Fcoef::create(const int i1, const int i2, const int i3)
{
     if(i1*i2*i3>0&&i1>0&&i2>0&&i3>0)
     {
         ind1 = i1;
         ind4 = i2;
         ind5 = i3;
         delete[] p;
         int tot = ind1*ind2*ind3*ind4*ind5;
         p = new complex<double> [tot];
         for(int i=0;i<tot;i++) p[i] = complex<double>(0.0,0.0);
     }
     else
     {
         cout<<"not allowed!"<<endl;
     }

     return;
}

Soc::Soc()
{
	m_loc = new Vector3<double> [1];
	angle1 = new double[1];
	angle2 = new double[1];
	p_rot = new complex<double> [1];
}

Soc::~Soc()
{
	delete[] m_loc;
	delete[] angle1;
	delete[] angle2;
	delete[] p_rot;
}


double Soc::spinor(const int l, const double j, const int m, const int spin)
{
	double den;    // denominator
	double spinor0;

	if (spin!=0&&spin!=1) WARNING_QUIT("spinor","spin direction unknown");
	if (m<-l-1||m>l) WARNING_QUIT("spinor","m not allowed");

	den = 1.0/(2.0*l+1.0);
	if (fabs(j-l-0.50)<1e-8) 
	{
		//j == l+0.5
		if (spin==0) 
		{
			spinor0= sqrt((l+m+1.0)*den);
		}
		if (spin==1) 
		{
			spinor0= sqrt((l-m)*den);
		}
	}
	else if (fabs(j-l+0.50)<1e-8){
		//j == l-0.5
		if (m<-l+1)
		{
			spinor0=0.0;
		}
		else
		{
			if (spin==0) spinor0= sqrt((l-m+1.0)*den);
			if (spin==1) spinor0= -sqrt((l+m)*den);
		}
	}
	else
	{
		WARNING_QUIT("soc::spinor","j and l not compatible");
	}

	return spinor0;
}

void Soc::rot_ylm(const int lmax)
{
	int m,n,l;
	l_max_ = lmax;
	l2plus1_ = 2 * l_max_ + 1;
	delete[] p_rot;
	this->p_rot = new complex<double> [l2plus1_ * l2plus1_];
	ZEROS(p_rot, l2plus1_ * l2plus1_);
	this->p_rot[l_max_] = complex<double>(1.0, 0.0);

	l = l_max_;
	for(int i=1;i<2*l+1;i += 2)
	{
		m = (i+1)/2;
		n = l-m;
		this->p_rot[l2plus1_*i + n] = complex<double>(pow(-1.0 , m) / sqrt(2) , 0.0);
		this->p_rot[l2plus1_*(i+1) + n] = complex<double>(0.0,-pow(-1.0 , m)/sqrt(2));
		n = l+m;
		this->p_rot[l2plus1_*i + n] = complex<double>(1.0/sqrt(2), 0.0);
		this->p_rot[l2plus1_*(i+1) + n] = complex<double>(0.0, 1.0/sqrt(2));
	}
	return;
}

int Soc::sph_ind(const int l, const double j, const int m, const int spin)
{
    // This function calculates the m index of the spherical harmonic
    // in a spinor with orbital angular momentum l, total angular 
    // momentum j, projection along z of the total angular momentum m+-1/2. 
    // Spin selects the up (spin=1) or down (spin=2) coefficient.
    int sph_ind0;
    if (spin != 0 && spin!= 1)
	{
        WARNING_QUIT("sph_ind","spin must be 0 1");
	}
    if (m< -l-1 || m> l)
	{
        WARNING_QUIT("sph_ind","m not allowed");
	}
    if (fabs(j-l-0.5)< 1e-8)
	{
        if(spin == 0) sph_ind0 = m;
        if(spin == 1) sph_ind0 = m+1;
    }
    else if (fabs(j-l+0.5)<1e-8 )
	{
        if(m<-l+1) 
		{
			sph_ind0 = 0;
		}
        else 
		{
			if(spin ==0) 
			{
				sph_ind0 = m-1;
			}
			if(spin ==1) 
			{
				sph_ind0 = m;
			}
        }
    }
    else 
	{
		cout<<"l= "<<l<<" j= "<<j<<endl;
		WARNING_QUIT("sph_ind","l and j not suitable");
	}
    if(sph_ind0 < -l || sph_ind0 >l )
	{
		sph_ind0 =0;
	}
    return sph_ind0;
}


void Soc::cal_ux(const int ntype)
{
	double amag, uxmod;
	int starting_it;
	bool is_paraller;
	//do not sign feature in teh general case
	lsign = false;
	ZEROS(ux, 3);

	starting_it = 0;
	for(int it = 0;it<ntype;it++)
	{
		amag = pow(m_loc[it].x,2) + pow(m_loc[it].y,2) + pow(m_loc[it].z,2);
		if(amag > 1e-6)
		{
			ux[0] = m_loc[it].x;
			ux[1] = m_loc[it].y;
			ux[2] = m_loc[it].z;
			starting_it = it;
			lsign = true;
			break;
		}
	}
	//initial magnetizations should be parallel
	for(int it = starting_it+1; it<ntype;it++)
	{
		lsign = lsign && judge_parallel(ux, m_loc[it]);
	}
	if(lsign)
	{
		uxmod =  pow(ux[0],2) + pow(ux[1],2) +pow(ux[2],2);
		if(uxmod<1e-6) 
		{
			WARNING_QUIT("cal_ux","wrong uxmod");
		}
		for(int i = 0;i<3;i++)
		{
			ux[i] *= 1/sqrt(uxmod);
		}
		//       cout<<"    Fixed quantization axis for GGA: "
		//<<setw(10)<<ux[0]<<"  "<<setw(10)<<ux[1]<<"  "<<setw(10)<<ux[2]<<endl;
	}
	return;
}


bool Soc::judge_parallel(double a[3], Vector3<double> b)
{
   bool jp=false;
   double cross;
   cross = pow((a[1]*b.z-a[2]*b.y),2) +  pow((a[2]*b.x-a[0]*b.z),2) + pow((a[0]*b.y-a[1]*b.x),2);
   jp = (fabs(cross)<1e-6);
   return jp;
}

void Soc::set_fcoef(
	const int &l1,
	const int &l2,
	const int &is1,
	const int &is2,
	const int &m1,
	const int &m2,
	const double &j1,
	const double &j2,
	const int &it,
	const int &ip1,
	const int &ip2)
{
	complex<double> coeff = complex<double>(0.0,0.0);
	for(int m=-l1-1;m<l1+1;m++)
	{
		const int mi = sph_ind(l1,j1,m,is1) + l_max_ ;
		const int mj = sph_ind(l2,j2,m,is2) + l_max_ ;
		coeff += rotylm(m1,mi) * spinor(l1,j1,m,is1)*
				conj(rotylm(m2,mj)) * spinor(l2,j2,m,is2);
	}
	this->fcoef(it,is1,is2,ip1,ip2) = coeff;
}