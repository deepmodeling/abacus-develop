#include"vdwd3.h"
#include"../src_global/element_name.h"
#include"../src_global/constants.h"
#include"../src_global/global_function.h"
#include<cstring>
#include<algorithm>

Vdwd3::Vdwd3( const UnitCell_pseudo &unitcell):
	    energy_result(0),
		noabc(true),
		alp(14.0),
		max_elem(94),
		maxc(5),
		nline(32385),
        ucell(unitcell),
		init_set(false)
{	
	disp = 0.0;
	init_C6_tmp();
	init_r2r4();
	add_r2r4=VECTOR_TO_PTR(r2r4);
	init_rcov();
	add_rcov=VECTOR_TO_PTR(rcov);
}

void Vdwd3::atomkind (const UnitCell_pseudo &unitcell)
{
	atom_kind.resize(unitcell.ntype) ;
	for ( vector<int>::size_type i=0 ; i!=unitcell.ntype ; i++ )
	{
		for (int j=0; j!=element_name.size(); j++)
		{
			if ( unitcell.atoms[i].psd == element_name[j] )
			{
				atom_kind[i] = j ;
				break;
			}
		}
	}
	return ;
}


void Vdwd3::XYZ( const UnitCell_pseudo &unitcell,double (*xyz)[3],int *iz)
{
	int k = 0;
	for(int it=0 ; it<unitcell.ntype ; it++)
	{
		for(int ia=0 ; ia<unitcell.atoms[it].na ; ia++)
		{
			iz[k] = atom_kind[it]+1;
			//xyz in Cartesian_au
			xyz[k][0] = unitcell.atoms[it].tau[ia].x * unitcell.lat0;
			xyz[k][1] = unitcell.atoms[it].tau[ia].y * unitcell.lat0;
			xyz[k][2] = unitcell.atoms[it].tau[ia].z * unitcell.lat0;
			k += 1;
		}
	}
        return;
}
void Vdwd3::init_mxc()
{
	for(int i = 0 ; i < max_elem; i++)
	{
		mxc[i] = 1;
	}
	return ;
}

int Vdwd3::limit(int &i)
{
	int icn = 1;
	while(i>=100)
	{
		i-=100;
		icn+=1;
	}
	return icn;
}

void Vdwd3::loadc6()
{
	init_mxc();
	int k = 0,iat,jat,iatcn,jatcn;
	for(int n=0;n<nline;n++)
	{
		k = n*5;
		
		iat = static_cast<int>(C6_tmp[k+1])-1;
		jat = static_cast<int>(C6_tmp[k+2])-1;
        iatcn = limit(iat);
        jatcn = limit(jat);

		mxc[iat] = max(mxc[iat],iatcn);
		mxc[jat] = max(mxc[jat],jatcn);

		c6ab[0][jatcn-1][iatcn-1][jat][iat] = C6_tmp[k];
		c6ab[1][jatcn-1][iatcn-1][jat][iat] = C6_tmp[k+3];
		c6ab[2][jatcn-1][iatcn-1][jat][iat] = C6_tmp[k+4];

		c6ab[0][iatcn-1][jatcn-1][iat][jat] = C6_tmp[k];
		c6ab[1][iatcn-1][jatcn-1][iat][jat] = C6_tmp[k+4];
		c6ab[2][iatcn-1][jatcn-1][iat][jat] = C6_tmp[k+3];

	}
	return;
}

void Vdwd3::setlat( const UnitCell_pseudo &unitcell)
{
        lat[0][0] = unitcell.latvec.e11* unitcell.lat0;
        lat[0][1] = unitcell.latvec.e12* unitcell.lat0;
        lat[0][2] = unitcell.latvec.e13* unitcell.lat0;

        lat[1][0] = unitcell.latvec.e21* unitcell.lat0;
        lat[1][1] = unitcell.latvec.e22* unitcell.lat0;
        lat[1][2] = unitcell.latvec.e23* unitcell.lat0;

        lat[2][0] = unitcell.latvec.e31* unitcell.lat0;
        lat[2][1] = unitcell.latvec.e32* unitcell.lat0;
        lat[2][2] = unitcell.latvec.e33* unitcell.lat0;

        return;
}

/* void Vdwd3::setfuncpar(int &ver)
{
	if(ver==1)
	{
	    s6 = 1.0;
        rs6 = 1.217;
        s18 = 0.722;
        rs18 = 1.0;
	}
	else if(ver==2)
    {
	    s6 = 1.0;
        rs6 = 0.4289;
        s18 = 0.7875;
        rs18 = 4.4407;
	}
	return;
} */

void Vdwd3::initset()
{
	if(!init_set)
	{
		if(abc)
		{
			noabc = !abc;
		}
		atomkind(ucell);
	    loadc6();
	    setlat(ucell);
		setr0ab_(&max_elem,&BOHR_TO_A,r0ab);
/*         if((s6==0)&&(rs6==0)&&(s18==0)&&(rs18==0)&&(rthr2==0))
		{
			setfuncpar(version);
			rthr2=95**2;
		}	 */		
		for(int i=0 ; i<3 ; i++)
        {
            tau_max[i] = 0;
        }
		if(model=="radius")
		{
			set_criteria_(&rthr2,lat,tau_max);
		    for(int i=0 ; i<3 ; i++)
		    {
			    rep_vdw[i] = int(tau_max[i])+1;
		    }
		}
		set_criteria_(&cn_thr2,lat,tau_max);
        for(int i=0 ; i<3 ; i++)
        {
            rep_cn[i] = int(tau_max[i])+1;
        }
		for(int i=0 ; i<3 ; i++)
		{
			for(int j=0 ; j<3 ; j++)
			{
				stress[i][i] = 0;
				sigma[i][i] = 0;
			}
		}

		init_set = true;
	}
        return;
}

double Vdwd3::energy()
{
	TITLE("Vdwd3","energy");
    initset();
    double xyz[ucell.nat][3];
	int iz[ucell.nat];
	
	XYZ(ucell,xyz,iz);
	
	double alp2 = alp+2;
	double alp4 = alp+4;
    pbcedisp_(&max_elem,&maxc,&ucell.nat,xyz,iz,c6ab,mxc,add_r2r4,
		r0ab,add_rcov,&rs6,&rs18,&rs18,&alp,&alp2,&alp4,
	    &version,&noabc,&e6,&e8,&e10,&e12,&e63,lat,
		&rthr2,rep_vdw,&cn_thr2,rep_cn);
 
	energy_result = (-s6*e6-s18*e8-e63)*2;
	return energy_result;
}

vector< vector<double> > Vdwd3::force(
	const bool force_for_vdw, 
	const bool stress_for_vdw, 
	matrix &force_vdw, 
	matrix &stress_result
)
{
	TITLE("Vdwd3","force");
    initset();
	double xyz[ucell.nat][3];
	int iz[ucell.nat];
	XYZ(ucell,xyz,iz);
	double g[ucell.nat][3];
    for(int iat=0; iat!=ucell.nat; ++iat)
	{
		for(int j=0; j!=3; ++j)
		{
			g[iat][j] = 0;
		}
	}
	
	force_result.resize(ucell.nat);
	for( int iat=0; iat!=ucell.nat; ++iat )
	{
		force_result[iat].assign(3,0);
	}
	if(stress_for_vdw) stress_result.zero_out();
	
	double alp2 = alp+2;
	double alp4 = alp+4;
	pbcgdisp_(&max_elem,&maxc,&ucell.nat,xyz,iz,c6ab,mxc,add_r2r4,
			r0ab,add_rcov,&s6,&s18,&rs6,&rs18,&rs18,&alp,&alp2,&alp4,
            &noabc,&version,g,&disp,stress,sigma,lat,rep_vdw,rep_cn,
			&rthr2,&cn_thr2);
			
	for(int iat=0 ; iat<ucell.nat ; iat++)
	{
		for(int i=0 ; i<3 ; i++)
		{
			force_result[iat][i] = -2*g[iat][i];			
		}
	}
	for(int ipol=0;ipol<3;ipol++)
	{
		if(stress_for_vdw)
		{
			for(int jpol=0;jpol<3;jpol++)
			{
				stress_result(ipol,jpol) = 2*sigma[ipol][jpol]/ucell.omega;
			}
		}
		if(force_for_vdw)
		{
			for(int iat=0;iat<ucell.nat;iat++)
			{
				force_vdw(iat,ipol) += force_result[iat][ipol];
			}
		}
	}
	return force_result;
}


