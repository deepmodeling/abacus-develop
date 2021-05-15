#include"MD_thermo.h"

MD_thermo::MD_thermo()
{
    G = new Vector3<double>[1];
    NHCeta = new Vector3<double>[1];
    NHCpeta = new Vector3<double>[1];
}

MD_thermo::~MD_thermo()
{
    delete[] G;
    delete[] NHCeta;
    delete[] NHCpeta;
}

void MD_thermo::init_NHC(
    const int &MNHC_in, 
    const double &Qmass_in, 
    const double &NVT_tau_in, 
    const double &dt_in,
    const int &NVT_control, 
    ofstream &ofs, 
    const int &numIon,
    const double &temperature,
    const Vector3<double>* vel,
    const double* allmass
    )
{
    this->MNHC_ = MNHC_in;
    this->Qmass_ = Qmass_in;
    this->NVT_tau_ = NVT_tau_in;
    this->dt_ = dt_in;
    this->numIon_ = numIon;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
	init_by_array(init, length);
	ofs<<" ...............Nose-Hoover Chain parameter initialization...............  " << endl;
	ofs<<" Temperature =    "<< temperature << endl;
	ofs<<" Temperature2 =    "<< temperature/K_BOLTZMAN_AU << endl;
	ofs<<" NHC frequency =    "<< 1.0/NVT_tau_ << endl;
	ofs<<" NHC chain =    "<< MNHC_ << endl;
	ofs<<" Qmass  =    "<< Qmass_ << endl;
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
		delta[i]=w[i]*dt_;
	}
		
    delete[] G;	
	G=new Vector3<double>[MNHC_*numIon_];
	delete[] NHCeta;	
	NHCeta=new Vector3<double>[MNHC_*numIon_];
	delete[] NHCpeta;	
	NHCpeta=new Vector3<double>[MNHC_*numIon_];

	for(int j=0;j<MNHC_;j++)
    {
		for(int k=0;k<numIon_;k++)
        {
			NHCeta[j*numIon_+k].x=genrand_real2()-0.5; 
			NHCeta[j*numIon_+k].y=genrand_real2()-0.5; 
			NHCeta[j*numIon_+k].z=genrand_real2()-0.5;
			NHCpeta[j*numIon_+k].x=gaussrand()*sqrt(Qmass_ * temperature);
			NHCpeta[j*numIon_+k].y=gaussrand()*sqrt(Qmass_ * temperature);
			NHCpeta[j*numIon_+k].z=gaussrand()*sqrt(Qmass_ * temperature);
		}
	}
		
	for(int j=MNHC_-1;j>=1;j--)
    {
		for(int k=0;k<numIon_;k++)
        {
			G[j*numIon_+k].x=pow(NHCpeta[(j-1)*numIon_+k].x,2)/Qmass_ - 1.0 * temperature;
			G[j*numIon_+k].y=pow(NHCpeta[(j-1)*numIon_+k].y,2)/Qmass_ - 1.0 * temperature;
			G[j*numIon_+k].z=pow(NHCpeta[(j-1)*numIon_+k].z,2)/Qmass_ - 1.0 * temperature;
		}
	}
		
	for(int k=0;k<numIon_;k++)
    {
		G[0*numIon_+k].x=pow(vel[k].x,2)*allmass[k]-1.0 * temperature;
		G[0*numIon_+k].y=pow(vel[k].y,2)*allmass[k]-1.0 * temperature;
		G[0*numIon_+k].z=pow(vel[k].z,2)*allmass[k]-1.0 * temperature;
	}
	ofs << "finish NHC thermostat define" << endl; //zifei
}

double MD_thermo::NHChamiltonian(
    const double &KE,
    const double &PE,
    const double &temperature,
    const int &nfrozen )
{
	double NHChamiltonian0; // The conserved quantity

	NHChamiltonian0 = KE + PE ;
//	cout << "hamiltonian0 =    "<< NHChamiltonian0 <<endl;

	for(int i=0;i<numIon_;i++)
    {
		for(int j=0;j<MNHC_;j++)
        {
			NHChamiltonian0 += pow(NHCpeta[i*MNHC_+j].x,2)/2.0/Qmass_+temperature*NHCeta[i*MNHC_+j].x;
			NHChamiltonian0 += pow(NHCpeta[i*MNHC_+j].y,2)/2.0/Qmass_+temperature*NHCeta[i*MNHC_+j].y;
			NHChamiltonian0 += pow(NHCpeta[i*MNHC_+j].z,2)/2.0/Qmass_+temperature*NHCeta[i*MNHC_+j].z;
		}
	}

	if (!MY_RANK)
    {
        ofs_running<< " --------------------------------------------------"<<endl;
        ofs_running<< " SUMMARY OF NVT CALCULATION"<<endl;
        ofs_running<<" --------------------------------------------------"<<endl;
        ofs_running<<" NVT Conservation     : "<<setw(10)<< NHChamiltonian0*2<<" (Rydberg)"<<endl;
        ofs_running<<" NVT Temperature      : "<<setw(10)<< KE*2/(3*double(numIon_-nfrozen))/K_BOLTZMAN_AU<<" (K)"<<endl;
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


void MD_thermo::Integrator(
    const int control,
    const double &temperature,
    Vector3<double>* vel,
    const double* allmass)
{
	if(control == 1) NHCIntegrator(temperature, vel, allmass);
	else if(control == 2) LGVIntegrator(temperature, vel, allmass);
	else if(control == 3) ADSIntegrator(temperature, vel, allmass);
	else WARNING_QUIT("MD_thermo:Integrator", "please choose available reservoir!!!");
	return;
}

void MD_thermo::LGVIntegrator(
    const double &temperature,
    Vector3<double>* vel,
    const double* allmass
)
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function propagates the Langevin dynamics 
//--------------------------------------------------------------------
	double randomx,randomy,randomz;
	double tempV;
	double c1k = exp(-1/NVT_tau_ * dt_);
	//c1k=e^(-gamma*dt)
	double c2k=sqrt(1.0-c1k*c1k);

    for(int iatom=0;iatom<numIon_;iatom++)
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
void MD_thermo::ADSIntegrator(
    const double &temperature,
    Vector3<double>* vel,
    const double* allmass
)
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function propagates the Andersen thermostat
//--------------------------------------------------------------------

    double uniform_random;
    double ranx;
    double rany;
    double ranz;

	double nu = 1/ NVT_tau_;

    for(int iatom=0;iatom<numIon_;iatom++)
    {
        uniform_random = genrand_real2() ;

        if(uniform_random < 1.0-exp(-nu*dt_) ) 
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
void MD_thermo::NHCIntegrator(
    const double &temperature,
    Vector3<double>* vel,
    const double* allmass
){
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function propagates the Nose-Hoover Chain
//   NHCpeta :NHC momentum
//   NHCeta : NHC position
//   Q : NHC mass
//--------------------------------------------------------------------

    for(int i=0;i<numIon_;i++)
    {
        G[0*numIon_+i].x=pow(vel[i].x,2)*allmass[i]-1.0 * temperature;
        G[0*numIon_+i].y=pow(vel[i].y,2)*allmass[i]-1.0 * temperature;
        G[0*numIon_+i].z=pow(vel[i].z,2)*allmass[i]-1.0 * temperature;
    }
  
    for(int i=0;i<nsy;i++)
    {
        for(int j=0;j<numIon_;j++)
        {
            NHCpeta[(MNHC_-1)*numIon_+j] += G[(MNHC_-1)*numIon_+j]*delta[i]/2.0;
        }

        for(int j=MNHC_-2;j>=0;j--)
        {
            for(int k=0;k<numIon_;k++)
            {
                NHCpeta[j*numIon_+k].x *=exp(-NHCpeta[(j+1)*numIon_+k].x/Qmass_*delta[i]/4.0);
                NHCpeta[j*numIon_+k].y *=exp(-NHCpeta[(j+1)*numIon_+k].y/Qmass_*delta[i]/4.0);
                NHCpeta[j*numIon_+k].z *=exp(-NHCpeta[(j+1)*numIon_+k].z/Qmass_*delta[i]/4.0);

                NHCpeta[j*numIon_+k]+=G[j*numIon_+k]*delta[i]/2.0;

                NHCpeta[j*numIon_+k].x *=exp(-NHCpeta[(j+1)*numIon_+k].x/Qmass_*delta[i]/4.0);
                NHCpeta[j*numIon_+k].y *=exp(-NHCpeta[(j+1)*numIon_+k].y/Qmass_*delta[i]/4.0);
                NHCpeta[j*numIon_+k].z *=exp(-NHCpeta[(j+1)*numIon_+k].z/Qmass_*delta[i]/4.0);
            }
        }//MNHC-2
  	
        for(int j=MNHC_-1;j>=0;j--)
        {
            for(int k=0;k<numIon_;k++)
            {
                NHCeta[j*numIon_+k].x += NHCpeta[j*numIon_+k].x/Qmass_*delta[i];
                NHCeta[j*numIon_+k].y += NHCpeta[j*numIon_+k].y/Qmass_*delta[i];
                NHCeta[j*numIon_+k].z += NHCpeta[j*numIon_+k].z/Qmass_*delta[i];
            }
        }
        //p=p*dexp(-peta(1,:)/Q(1,:)*delta(i))
        for(int k=0;k<numIon_;k++)
        {
            vel[k].x *= exp(-NHCpeta[0*numIon_+k].x/Qmass_*delta[i]);
            vel[k].y *= exp(-NHCpeta[0*numIon_+k].y/Qmass_*delta[i]);
            vel[k].z *= exp(-NHCpeta[0*numIon_+k].z/Qmass_*delta[i]);
        }

  	//G(1,:)=p**2/mass-1/beta
  	    for(int k=0;k<numIon_;k++){
            G[0*numIon_+k].x=pow(vel[0*numIon_+k].x,2)*allmass[k]-1.0 * temperature;
            G[0*numIon_+k].y=pow(vel[0*numIon_+k].y,2)*allmass[k]-1.0 * temperature;
            G[0*numIon_+k].z=pow(vel[0*numIon_+k].z,2)*allmass[k]-1.0 * temperature;
  	    }
  		
        for(int j=0;j<MNHC_-1;j++)
        {
            for(int k=0;k<numIon_;k++)
            {
                NHCpeta[j*numIon_+k].x *=exp(-NHCpeta[(j+1)*numIon_+k].x/Qmass_*delta[i]/4.0);
                NHCpeta[j*numIon_+k].y *=exp(-NHCpeta[(j+1)*numIon_+k].y/Qmass_*delta[i]/4.0);
                NHCpeta[j*numIon_+k].z *=exp(-NHCpeta[(j+1)*numIon_+k].z/Qmass_*delta[i]/4.0);

                NHCpeta[j*numIon_+k]+=G[j*numIon_+k]*delta[i]/2.0;

                NHCpeta[j*numIon_+k].x *=exp(-NHCpeta[(j+1)*numIon_+k].x/Qmass_*delta[i]/4.0);
                NHCpeta[j*numIon_+k].y *=exp(-NHCpeta[(j+1)*numIon_+k].y/Qmass_*delta[i]/4.0);
                NHCpeta[j*numIon_+k].z *=exp(-NHCpeta[(j+1)*numIon_+k].z/Qmass_*delta[i]/4.0);

                G[(j+1)*numIon_+k].x = pow(NHCpeta[j*numIon_+k].x,2)/Qmass_-1.0 * temperature;
                G[(j+1)*numIon_+k].y = pow(NHCpeta[j*numIon_+k].y,2)/Qmass_-1.0 * temperature;
                G[(j+1)*numIon_+k].z = pow(NHCpeta[j*numIon_+k].z,2)/Qmass_-1.0 * temperature;
            }
        }

        for(int k=0;k<numIon_;k++)
        {
            NHCpeta[(MNHC_-1)*numIon_+k]+=G[(MNHC_-1)*numIon_+k]*delta[i]/2.0;
        }			
    }//nsy

    return;
}

void MD_thermo::NHC_info_out(const int& step, const int& recordFreq, const int& numIon)
{
    bool pass;
	pass = 0;

	if (recordFreq==1||step==1||( recordFreq > 1&& step%recordFreq==0 ) )
		pass =1;
	if (!pass) return;

	if(!MY_RANK){
		stringstream ssc;
		ssc << global_out_dir << "Restart_md.dat";
		ofstream file(ssc.str().c_str(), ios::app);
        file<<'\n';
		file<<"MD_THERMOSTAT"<<endl;
		file<<"MNHC: "<<MNHC_<<endl;
		file<<"G: "<<endl;
		for(int i=0;i<numIon_*MNHC_;i++){
			file<<setprecision (12)<<G[i].x<<" "<<setprecision (12)<<G[i].y<<" "<<setprecision (12)<<G[i].z<<endl;
		}
        file<<"NHCeta: "<<endl;
		for(int i=0;i<numIon_*MNHC_;i++){
			file<<setprecision (12)<<NHCeta[i].x<<" "<<setprecision (12)<<NHCeta[i].y<<" "<<setprecision (12)<<NHCeta[i].z<<endl;
		}
        file<<"NHCpeta: "<<endl;
		for(int i=0;i<numIon_*MNHC_;i++){
			file<<setprecision (12)<<NHCpeta[i].x<<" "<<setprecision (12)<<NHCpeta[i].y<<" "<<setprecision (12)<<NHCpeta[i].z<<endl;
		}                
		file.close();
	}
	return;
    return;
}

#include<cstring>
void MD_thermo::NHC_restart()
{
    int error=0;
    double *nhcg=new double[numIon_*3*MNHC_];
    double *nhcp=new double[numIon_*3*MNHC_];
    double *nhce=new double[numIon_*3*MNHC_];
    if (!MY_RANK)
	{
		stringstream ssc;
		ssc << global_readin_dir << "Restart_md.dat";
		ifstream file(ssc.str().c_str());

        char word[80];
        int mnhc;
        while(file.good())//search and read MNHC
        {
            file>>word;
            if(std::strcmp("MNHC:",word)==0)
            {
                file>>mnhc;
                break;
            }
            else
            {
                file.ignore(150, '\n');
            }
        }

        if(mnhc!=MNHC_)
        {
            cout<<"please ensure whether 'Restart_md.dat' right!"<<endl;
            error = 1;
        }

        if(!error)
        {
            file.get();
            file.ignore(3, '\n');//read G
            for(int i = 0;i<numIon_*3*MNHC_;i++)
            {
                file>>nhcg[i];
            }
            file.get();
            file.ignore(8, '\n');//read eta
            for(int i = 0;i<numIon_*3*MNHC_;i++)
            {
                file>>nhce[i];
            }
            file.get();
            file.ignore(9, '\n');//read peta
            for(int i = 0;i<numIon_*3*MNHC_;i++)
            {
                file>>nhcp[i];
            }
            file.close();
        }

    }
#ifdef __MPI
    MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    if(error)
    {
        delete[] nhcg;
        delete[] nhcp;
        delete[] nhce;
        exit(0);
    }
#ifdef __MPI
	MPI_Bcast(nhcg,numIon_*3*MNHC_,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(nhcp,numIon_*3*MNHC_,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(nhce,numIon_*3*MNHC_,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

    for(int i=0;i<numIon_*MNHC_;i++)
    {
        G[i].x = nhcg[i*3];
        G[i].y = nhcg[i*3+1];
        G[i].z = nhcg[i*3+2];
        NHCpeta[i].x = nhcp[i*3];
        NHCpeta[i].y = nhcp[i*3+1];
        NHCpeta[i].z = nhcp[i*3+2];
        NHCeta[i].x = nhce[i*3];
        NHCeta[i].y = nhce[i*3+1];
        NHCeta[i].z = nhce[i*3+2];
    }
    delete[] nhcg;
    delete[] nhcp;
    delete[] nhce;
    return;
}