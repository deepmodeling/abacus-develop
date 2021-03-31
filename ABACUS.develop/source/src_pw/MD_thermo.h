#ifndef MD_THERMO_H
#define MD_THERMO_H

class MD_thermo 
{
    public:
    MD_thermo();
    ~MD_thrmo();

    void Integrator(const int control);
    void init_NHC(
        const int &MNHC_in, 
        const int &NVT_control, 
        ofstream &ofs, 
        const int &numIon
        );
    double NHChamiltonian(
        const double &KE,
        const double &PE,
        const int &numIon,
        const double &temperature,
        const int &nfrozen );

    private:    
    void ADSIntegrator();
	void LGVIntegrator();
	void NHCIntegrator();

	double gaussrand();
	
	void init_genrand(unsigned long s);
	void init_by_array(unsigned long init_key[], int key_length);
	unsigned long genrand_int32(void);
	long genrand_int31(void);
	double genrand_real1(void);
	double genrand_real2(void);
	double genrand_real3(void);
	double genrand_res53(void) ;

    //NVT thermostat parameters
	const static int N_mt19937=624;
	const static int M_mt19937=397;
	unsigned long const MATRIX_A=0x9908b0dfUL ;  /* constant vector a */
	unsigned long const UPPER_MASK=0x80000000UL; /* most significant w-r bits */
	unsigned long const LOWER_MASK=0x7fffffffUL; /* least significant r bits */
	unsigned long mt[N_mt19937]; /* the array for the state vector  */
	int mti=625; /* mti==N+1 means mt[N] is not initialized */
	long double gamma; //langevin friction coefficient
	long double nu; //Andersen collision frequency
	long long double c1k; //parameter in Langevin
	double c2k; //parameter in Langevin
    const static int nsy=7; //parameter in NHC, constant, no need to modification
    long double w[nsy]; //parameters in NHC
    long double delta[nsy]; //parameters in NHC

    //input parameters
    int MNHC;
    int Qmass;

    //need to be allocated
    Vector3<double> *G; //parameter in NHC
    Vector3<double> *NHCeta; //NHC position
    Vector3<double> *NHCpeta; //NHC momentum

};