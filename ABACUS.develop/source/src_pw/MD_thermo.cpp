#include"MD_thermo.h"

MD_thermo::init_NHC(
    const int &MNHC_in, 
    const int &NVT_control, 
    ofstream &ofs, 
    const int &numIon
    )
{
    this->MNHC = MNHC_in;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	init_by_array(init, length);
	ofs<<" ...............Nose-Hoover Chain parameter initialization...............  " << endl;
	ofs<<" Temperature =    "<< temperature << endl;
	ofs<<" Temperature2 =    "<< temperature/K_BOLTZMAN_AU << endl;
	ofs<<" NHC frequency =    "<< 1.0/NVT_tau << endl;
	ofs<<" NHC chain =    "<< MNHC << endl;
	ofs<<" Qmass  =    "<< Qmass << endl;
	ofs<<" ...............................................................  " << endl;
	w[0]=0.7845136105;
	w[6]=0.7845136105;
	w[1]=0.2355732134;
	w[5]=0.2355732134;
	w[2]=-1.177679984;
	w[4]=-1.177679984;
	w[3]=1-w[0]-w[1]-w[2]-w[4]-w[5]-w[6];
	//
	for(int i=0;i<nsy;i++)
    {
		delta[i]=w[i]*dt;
	}
		
    delete[] G;	
	G=new Vector3<double>[MNHC*numIon];
	delete[] NHCeta;	
	NHCeta=new Vector3<double>[MNHC*numIon];
	delete[] NHCpeta;	
	NHCpeta=new Vector3<double>[MNHC*numIon];

	for(int j=0;j<MNHC;j++)
    {
		for(int k=0;k<numIon;k++)
        {
			NHCeta[j*numIon+k].x=genrand_real2()-0.5; 
			NHCeta[j*numIon+k].y=genrand_real2()-0.5; 
			NHCeta[j*numIon+k].z=genrand_real2()-0.5;
			NHCpeta[j*numIon+k].x=gaussrand()*sqrt(Qmass * temperature);
			NHCpeta[j*numIon+k].y=gaussrand()*sqrt(Qmass * temperature);
			NHCpeta[j*numIon+k].z=gaussrand()*sqrt(Qmass * temperature);
		}
	}
		
	for(int j=MNHC-1;j>=1;j--)
    {
		for(int k=0;k<numIon;k++)
        {
			G[j*numIon+k].x=pow(NHCpeta[(j-1)*numIon+k].x,2)/Qmass - 1.0 * temperature;
			G[j*numIon+k].y=pow(NHCpeta[(j-1)*numIon+k].y,2)/Qmass - 1.0 * temperature;
			G[j*numIon+k].z=pow(NHCpeta[(j-1)*numIon+k].z,2)/Qmass - 1.0 * temperature;
		}
	}
		
	for(int k=0;k<numIon;k++)
    {
		G[0*numIon+k].x=pow(vel[k].x,2)*allmass[k]-1.0 * temperature;
		G[0*numIon+k].y=pow(vel[k].y,2)*allmass[k]-1.0 * temperature;
		G[0*numIon+k].z=pow(vel[k].z,2)*allmass[k]-1.0 * temperature;
	}
	ofs << "finish NHC thermostat define" << endl; //zifei
}

double MD_thermo::NHChamiltonian(
    const double &KE,
    const double &PE,
    const int &numIon,
    const double &temperature,
    const int &nfrozen )
{
	double NHChamiltonian0; // The conserved quantity

	NHChamiltonian0 = KE + PE ;
//	cout << "hamiltonian0 =    "<< NHChamiltonian0 <<endl;

	for(int i=0;i<numIon;i++)
    {
		for(int j=0;j<MNHC;j++)
        {
			NHChamiltonian0 += pow(NHCpeta[i*MNHC+j].x,2)/2.0/Qmass+temperature*NHCeta[i*MNHC+j].x;
			NHChamiltonian0 += pow(NHCpeta[i*MNHC+j].y,2)/2.0/Qmass+temperature*NHCeta[i*MNHC+j].y;
			NHChamiltonian0 += pow(NHCpeta[i*MNHC+j].z,2)/2.0/Qmass+temperature*NHCeta[i*MNHC+j].z;
		}
	}

	if (!MY_RANK)
    {
		out<< " "<<endl;
        out<< " --------------------------------------------------"<<endl;
        out<< " SUMMARY OF NVT CALCULATION"<<endl;
        out<<" --------------------------------------------------"<<endl;
        out<<" NVT Conservation     : "<<setw(10)<< NHChamiltonian0*2<<" (Rydberg)"<<endl;
        out<<" NVT Temperature      : "<<setw(10)<< KE*2/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<" (K)"<<endl;
        out<<" NVT Kinetic energy   : "<<setw(10)<< KE*2<<" (Rydberg)"<<endl;
        out<<" NVT Potential energy : "<<setw(10)<< PE*2<<" (Rydberg)"<<endl;
        ofs_running<< " "<<endl;
        ofs_running<< " --------------------------------------------------"<<endl;
        ofs_running<< " SUMMARY OF NVT CALCULATION"<<endl;
        ofs_running<<" --------------------------------------------------"<<endl;
        ofs_running<<" NVT Conservation     : "<<setw(10)<< NHChamiltonian0*2<<" (Rydberg)"<<endl;
        ofs_running<<" NVT Temperature      : "<<setw(10)<< KE*2/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<" (K)"<<endl;
        ofs_running<<" NVT Kinetic energy   : "<<setw(10)<< KE*2<<" (Rydberg)"<<endl;
        ofs_running<<" NVT Potential energy : "<<setw(10)<< PE*2<<" (Rydberg)"<<endl;
    }
//	cout << "hamiltonian1 =    "<< NHChamiltonian0 <<endl;
   
	return NHChamiltonian0;
}

double MD_thermo::gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) 
    {
        do 
        {
            double U1 = genrand_real2();
            double U2 = genrand_real2();
             
            V1 = 2.0 * U1 - 1.0;
            V2 = 2.0 * U2 - 1.0;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2.0 * log(S) / S);
    } 
    else
    {
        X = V2 * sqrt(-2.0 * log(S) / S);
    }

    phase = 1 - phase;
 
    return X;
}

//zifei
/* initializes mt[N] with a seed */
void MD_thermo::init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N_mt19937; mti++) 
    {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}


/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void MD_thermo::init_by_array(unsigned long init_key[], int key_length)
{
    init_genrand(19650218UL);
    int i=1; 
    int j=0;
    int k = (N_mt19937>key_length ? N_mt19937 : key_length);
    for (; k; k--)
    {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N_mt19937) { mt[0] = mt[N_mt19937-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N_mt19937-1; k; k--) 
    {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N_mt19937) { mt[0] = mt[N_mt19937-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long MD_thermo::genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N_mt19937) 
    { /* generate N words at one time */
        int kk;

        if (mti == N_mt19937+1)
        {   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */
        }
        for (kk=0;kk<N_mt19937-M_mt19937;kk++) 
        {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M_mt19937] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N_mt19937-1;kk++) 
        {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M_mt19937-N_mt19937)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N_mt19937-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N_mt19937-1] = mt[M_mt19937-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long MD_thermo::genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double MD_thermo::genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double MD_thermo::genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double MD_thermo::genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double MD_thermo::genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 


void MD_thermo::Integrator(const int control)
{
	if(control == 1) NHCIntegrator();
	else if(control == 2) LGVIntegrator();
	else if(control == 3) ADSIntegrator();
	else WARNING_QUIT("MD_thermo:Integrator", "please choose available reservoir!!!");
	return;
}

void MD_thermo::LGVIntegrator()
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function propagates the Langevin dynamics 
//--------------------------------------------------------------------
	double randomx,randomy,randomz;
	double tempV;
	double c1k = exp(-1/NVT_tau * dt);
	//c1k=e^(-gamma*dt)
	double c2k=sqrt(1.0-c1k*c1k);

    for(int iatom=0;iatom<numIon;iatom++)
    {
        randomx = gaussrand();
        tempV = vel[iatom].x;
        vel[iatom].x = c1k*tempV + c2k*sqrt(allmass[iatom] * temperature)*randomx/allmass[iatom];

        randomy = gaussrand();
        tempV = vel[iatom].y;
        vel[iatom].y = c1k*tempV + c2k*sqrt(allmass[iatom] * temperature)*randomy/allmass[iatom];

        randomz = gaussrand();
        tempV = vel[iatom].z;
        vel[iatom].z = c1k*tempV + c2k*sqrt(allmass[iatom] * temperature)*randomz/allmass[iatom];
    }


    return;
}


//added by zifei
void MD_thermo::ADSIntegrator()
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function propagates the Andersen thermostat
//--------------------------------------------------------------------

    double uniform_random;
    double ranx;
    double rany;
    double ranz;

	double nu = 1/ NVT_tau;

    for(int iatom=0;iatom<numIon;iatom++)
    {
        uniform_random = genrand_real2() ;

        if(uniform_random < 1.0-exp(-nu*dt) ) 
        {
        //if(uniform_random < 0.3 ) {
            ranx=sqrt(allmass[iatom] * temperature);
            ranx *= gaussrand();

            rany=sqrt(allmass[iatom] * temperature);
            rany *= gaussrand();

            ranz=sqrt(allmass[iatom] * temperature);
            ranz *= gaussrand();
            
            vel[iatom].x=ranx/allmass[iatom];
            vel[iatom].y=rany/allmass[iatom];
            vel[iatom].z=ranz/allmass[iatom];
        }
        else
        {//do nothing
            //vel[iatom]=vel[iatom];
        }
    }	

    return;
}


//zifei
void MD_thermo::NHCIntegrator(){
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function propagates the Nose-Hoover Chain
//   NHCpeta :NHC momentum
//   NHCeta : NHC position
//   Q : NHC mass
//--------------------------------------------------------------------

    for(int i=0;i<numIon;i++)
    {
        G[0*numIon+i].x=pow(vel[i].x,2)*allmass[i]-1.0 * temperature;
        G[0*numIon+i].y=pow(vel[i].y,2)*allmass[i]-1.0 * temperature;
        G[0*numIon+i].z=pow(vel[i].z,2)*allmass[i]-1.0 * temperature;
    }
  
    for(int i=0;i<nsy;i++)
    {
        for(int j=0;j<numIon;j++)
        {
            NHCpeta[(MNHC-1)*numIon+j] += G[(MNHC-1)*numIon+j]*delta[i]/2.0;
        }

        for(int j=MNHC-2;j>=0;j--)
        {
            for(int k=0;k<numIon;k++)
            {
                NHCpeta[j*numIon+k].x *=exp(-NHCpeta[(j+1)*numIon+k].x/Qmass*delta[i]/4.0);
                NHCpeta[j*numIon+k].y *=exp(-NHCpeta[(j+1)*numIon+k].y/Qmass*delta[i]/4.0);
                NHCpeta[j*numIon+k].z *=exp(-NHCpeta[(j+1)*numIon+k].z/Qmass*delta[i]/4.0);

                NHCpeta[j*numIon+k]+=G[j*numIon+k]*delta[i]/2.0;

                NHCpeta[j*numIon+k].x *=exp(-NHCpeta[(j+1)*numIon+k].x/Qmass*delta[i]/4.0);
                NHCpeta[j*numIon+k].y *=exp(-NHCpeta[(j+1)*numIon+k].y/Qmass*delta[i]/4.0);
                NHCpeta[j*numIon+k].z *=exp(-NHCpeta[(j+1)*numIon+k].z/Qmass*delta[i]/4.0);
            }
        }//MNHC-2
  	
        for(int j=MNHC-1;j>=0;j--)
        {
            for(int k=0;k<numIon;k++)
            {
                NHCeta[j*numIon+k].x += NHCpeta[j*numIon+k].x/Qmass*delta[i];
                NHCeta[j*numIon+k].y += NHCpeta[j*numIon+k].y/Qmass*delta[i];
                NHCeta[j*numIon+k].z += NHCpeta[j*numIon+k].z/Qmass*delta[i];
            }
        }
        //p=p*dexp(-peta(1,:)/Q(1,:)*delta(i))
        for(int k=0;k<numIon;k++)
        {
            vel[k].x *= exp(-NHCpeta[0*numIon+k].x/Qmass*delta[i]);
            vel[k].y *= exp(-NHCpeta[0*numIon+k].y/Qmass*delta[i]);
            vel[k].z *= exp(-NHCpeta[0*numIon+k].z/Qmass*delta[i]);
        }

  	//G(1,:)=p**2/mass-1/beta
  	    for(int k=0;k<numIon;k++){
            G[0*numIon+k].x=pow(vel[0*numIon+k].x,2)*allmass[k]-1.0 * temperature;
            G[0*numIon+k].y=pow(vel[0*numIon+k].y,2)*allmass[k]-1.0 * temperature;
            G[0*numIon+k].z=pow(vel[0*numIon+k].z,2)*allmass[k]-1.0 * temperature;
  	    }
  		
        for(int j=0;j<MNHC-1;j++)
        {
            for(int k=0;k<numIon;k++)
            {
                NHCpeta[j*numIon+k].x *=exp(-NHCpeta[(j+1)*numIon+k].x/Qmass*delta[i]/4.0);
                NHCpeta[j*numIon+k].y *=exp(-NHCpeta[(j+1)*numIon+k].y/Qmass*delta[i]/4.0);
                NHCpeta[j*numIon+k].z *=exp(-NHCpeta[(j+1)*numIon+k].z/Qmass*delta[i]/4.0);

                NHCpeta[j*numIon+k]+=G[j*numIon+k]*delta[i]/2.0;

                NHCpeta[j*numIon+k].x *=exp(-NHCpeta[(j+1)*numIon+k].x/Qmass*delta[i]/4.0);
                NHCpeta[j*numIon+k].y *=exp(-NHCpeta[(j+1)*numIon+k].y/Qmass*delta[i]/4.0);
                NHCpeta[j*numIon+k].z *=exp(-NHCpeta[(j+1)*numIon+k].z/Qmass*delta[i]/4.0);

                G[(j+1)*numIon+k].x = pow(NHCpeta[j*numIon+k].x,2)/Qmass-1.0 * temperature;
                G[(j+1)*numIon+k].y = pow(NHCpeta[j*numIon+k].y,2)/Qmass-1.0 * temperature;
                G[(j+1)*numIon+k].z = pow(NHCpeta[j*numIon+k].z,2)/Qmass-1.0 * temperature;
            }
        }

        for(int k=0;k<numIon;k++)
        {
            NHCpeta[(MNHC-1)*numIon+k]+=G[(MNHC-1)*numIon+k]*delta[i]/2.0;
        }			
    }//nsy

    return;
}