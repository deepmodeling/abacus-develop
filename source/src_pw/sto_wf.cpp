#include "sto_wf.h"
#include "global.h"

Stochastic_WF::Stochastic_WF()
{
    chiortho  = new ModuleBase::ComplexMatrix[1];
	chi0  = new ModuleBase::ComplexMatrix[1];
	emax_sto = 1;
	emin_sto = -1;
	stotype = "pw";
}

Stochastic_WF::~Stochastic_WF()
{ 
    delete[] chiortho;
    delete[] chi0;
}

void Stochastic_WF::init(const int & nchi_in,
    const int & nche_sto_in,
    const int & seed_sto_in,
    const double & emax_sto_in,
    const double & emin_sto_in,
    const std::string & stotype_in)
{
    ///initialize parameters

    this->nchi = nchi_in;
    this->nche_sto = nche_sto_in;
    this->seed_sto = seed_sto_in;
    this->emax_sto = emax_sto_in;
    this->emin_sto = emin_sto_in;
    this->stotype = stotype_in;
    //wait for init

    int ndim=0;
    if(stotype == "pw")
    {
        ndim = GlobalC::kv.ngk[0]; // only GAMMA point temporarily
    }
    else
    {
        ndim = GlobalC::pw.nrxx;
    }
    int totnpw=0;
#ifdef __MPI
    int * npwip = new int [GlobalV::NPROC_IN_POOL];
    int *rec = new int [GlobalV::NPROC_IN_POOL];
    int *displ = new int [GlobalV::NPROC_IN_POOL];
    for(int ip = 0 ; ip < GlobalV::NPROC_IN_POOL ;++ip)
    {
        rec[ip] = 1;
        displ[ip] = ip;
    }
    MPI_Allgatherv(&GlobalC::kv.ngk[0], 1, MPI_INT, npwip, rec, displ, MPI_INT, POOL_WORLD);
    for(int ip = 0; ip < GlobalV::NPROC_IN_POOL; ++ip)
    {
        totnpw += npwip[ip];
    }
#else
    totnpw = ndim;
#endif
    
    
    //distribute nchi for each process
    bool allbase = false;
    if(nchi == 0)
    {
        nchi = totnpw-GlobalV::NBANDS;
        std::cout<<"Using all normal bases: "<<totnpw<<std::endl;
        allbase = true;
    }
    nchip = int(nchi/GlobalV::KPAR);
    if(GlobalV::KPAR - GlobalV::MY_POOL - 1 < nchi%GlobalV::KPAR) ++nchip;

    std::complex<double> ui(0,1);

    //We temporarily init one group of orbitals for all k points.
    delete[] chi0;
    chi0 = new ModuleBase::ComplexMatrix[1]; 
    
    
    //srand((unsigned)time(NULL)+GlobalV::MY_RANK*10000);
    srand((unsigned)GlobalV::MY_RANK*10000);
    //srand((unsigned)0);
    
    if(allbase)
    {
        chi0[0].create(nchip,ndim,true);
        int re = GlobalV::KPAR - nchi % GlobalV::KPAR;
        int ip = 0, ig0 = 0;
        int ig;
        for(int i = 0 ; i < nchip ; ++i)
        {
#ifdef __MPI
            if(GlobalV::MY_POOL < re)
			{
                ig = GlobalV::MY_POOL * nchip + i - ig0;
			}
            else
			{
                ig = GlobalV::MY_POOL * nchip - re + i - ig0;
			}
            while(ig >= npwip[ip])
            {
                ig -= npwip[ip];
                ig0 += npwip[ip];
                ip++;
            }
            if(GlobalV::RANK_IN_POOL == ip)
			{
                chi0[0](i , ig) = 1;
			}
#else
                chi0[0](i , ig) = 1;
#endif
        }
    }
    else
    {
        //init with random number
        chi0[0].create(nchip,ndim,false);
        for(int i=0; i<chi0[0].size; ++i)
        {
            chi0[0].c[i]=exp(2*ModuleBase::PI*rand()/double(RAND_MAX)*ui) / sqrt(double(nchi));
        }
    }
    

    //

    

    delete[] chiortho;
#ifdef __MPI
    delete[] npwip;
    delete[] rec;
    delete[] displ;
#endif
    int nkk = 1; // We temporarily use gamma k point.
    chiortho = new ModuleBase::ComplexMatrix[1];
    if(GlobalV::NBANDS > 0)
    {
        chiortho[0].create(nchip,ndim,false);
    }
    
    return;
}


