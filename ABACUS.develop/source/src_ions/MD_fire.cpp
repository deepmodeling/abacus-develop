#include "src_ions/MD_fire.h"

MD_fire::MD_fire()
{
    dt_max=-1.0;
    alpha_start=0.10;
    alpha = alpha_start;

    finc=1.1;
    fdec=0.5;
    f_alpha=0.99;
    N_min=4;
    negative_count=0;
}

void MD_fire::check_FIRE(
    const int& numIon, 
    const Vector3<double>* force,
    double& deltaT, 
    Vector3<double>* vel)
{
    if(dt_max<0) dt_max = deltaT * 2.5; //initial dt_max
	double P(0.0);
	double sumforce(0.0);
	double normvel(0.0);
	
	for(int iatom =0;iatom<numIon;iatom++)
    {
		P += vel[iatom].x*force[iatom].x+vel[iatom].y*force[iatom].y+vel[iatom].z*force[iatom].z;
		
		sumforce += force[iatom].x*force[iatom].x+force[iatom].y*force[iatom].y+force[iatom].z*force[iatom].z;
		
		normvel += vel[iatom].x*vel[iatom].x +vel[iatom].y*vel[iatom].y +vel[iatom].z*vel[iatom].z ;
	}
	
	sumforce = sqrt(sumforce);
	normvel = sqrt(normvel);
	
	for(int iatom=0; iatom<numIon; iatom ++)
    {
		vel[iatom].x = (1.0-alpha)*vel[iatom].x+alpha*force[iatom].x/sumforce*normvel;
		vel[iatom].y = (1.0-alpha)*vel[iatom].y+alpha*force[iatom].y/sumforce*normvel;
		vel[iatom].z = (1.0-alpha)*vel[iatom].z+alpha*force[iatom].z/sumforce*normvel;
	}
	
	if(P > 0 )
	{
	    negative_count +=1 ;
		if(negative_count >=N_min)
        {
			deltaT=min(deltaT*finc,dt_max);
			alpha= alpha *f_alpha;
		}
	}
	else if( P<=0)
	{
		deltaT=deltaT*fdec;
		negative_count = 0;
		
		for(int iatom=0; iatom<numIon ; iatom++)
        {
			vel[iatom].x = 0.0;
			vel[iatom].y = 0.0;
			vel[iatom].z = 0.0;
		}
		
		alpha=alpha_start;
	}
}
