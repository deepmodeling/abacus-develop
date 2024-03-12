#include "../pw_basis.h"
#include "mpi.h"

using namespace std;


int main(int argc, char **argv)
{
    int nproc = 1, myrank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(myrank == 0) 
        cout<<"compare time between xprime and yprime"<<endl;
    int Nx = 0, Ny = 0, Nz =0;
//----------------------------------------------------
    ModuleBase::Matrix3 latvec(1, 0, 0, 0, 1, 0, 0, 0, 1);
    double lat0=10;
    // Nx = 100;
    // Ny = 50;
    // Nz =100;
    double wfcecut=50;
    double rhoecut=4 * wfcecut;
    int distribution_type = 1;
    int nbands = 100;
//----------------------------------------------------

    ModulePW::PW_Basis pwtest1;
    pwtest1.initmpi(nproc, myrank, MPI_COMM_WORLD);
    if(Nx*Ny*Nz==0)
    {
        pwtest1.initgrids(lat0,latvec,rhoecut);
        pwtest1.initparameters(false,wfcecut,distribution_type,true);
    }
    else
    {
        pwtest1.initgrids(lat0,latvec,Nx,Ny,Nz);
        pwtest1.initparameters(false,pwtest1.gridecut_lat*pwtest1.tpiba2/4,distribution_type,true);
    }
    pwtest1.setuptransform();

    ModulePW::PW_Basis pwtest2;
    pwtest2.initmpi(nproc, myrank, MPI_COMM_WORLD);
    if(Nx*Ny*Nz==0)
    {
        pwtest2.initgrids(lat0,latvec,rhoecut);
        pwtest2.initparameters(false,wfcecut,distribution_type,false);
    }
    else
    {
        pwtest2.initgrids(lat0,latvec,Nx,Ny,Nz);
        pwtest2.initparameters(false,pwtest1.gridecut_lat*pwtest1.tpiba2/4,distribution_type,false);
    }
    pwtest2.setuptransform();

    ModulePW::PW_Basis pwtest3;
    pwtest3.initmpi(nproc, myrank, MPI_COMM_WORLD);
    if(Nx*Ny*Nz==0)
    {
        pwtest3.initgrids(lat0,latvec,rhoecut);
        pwtest3.initparameters(true,wfcecut,distribution_type,true);
    }
    else
    {
        pwtest3.initgrids(lat0,latvec,Nx,Ny,Nz);
        pwtest3.initparameters(true,pwtest1.gridecut_lat*pwtest1.tpiba2/4,distribution_type,true);
    }
    pwtest3.setuptransform();

    ModulePW::PW_Basis pwtest4;
    pwtest4.initmpi(nproc, myrank, MPI_COMM_WORLD);
    if(Nx*Ny*Nz==0)
    {
        pwtest4.initgrids(lat0,latvec,rhoecut);
        pwtest4.initparameters(true,wfcecut,distribution_type,false);
    }
    else
    {
        pwtest4.initgrids(lat0,latvec,Nx,Ny,Nz);
        pwtest4.initparameters(true,pwtest1.gridecut_lat*pwtest1.tpiba2/4,distribution_type,false);
    }
    pwtest4.setuptransform();
    
    
    if(myrank == 0) 
    {
        cout<<"nx: "<<pwtest1.nx<<" ny: "<<pwtest1.ny<<" nz: "<<pwtest1.nz<<endl;
        cout<<"fullnpw   npw1: "<<pwtest1.npw<<" nst1: "<<pwtest1.nst<<";  npw2: "<<pwtest2.npw<<" nst2: "<<pwtest2.nst<<endl;
        cout<<"halfnpw   npw3: "<<pwtest3.npw<<" nst3: "<<pwtest3.nst<<";  npw4: "<<pwtest4.npw<<" nst4: "<<pwtest4.nst<<endl;
    }
    if(pwtest1.nx != pwtest2.nx|| pwtest2.nx != pwtest3.nx || pwtest3.nx != pwtest4.nx || 
       pwtest1.ny != pwtest2.ny|| pwtest2.ny != pwtest3.ny || pwtest3.ny != pwtest4.ny ||
       pwtest1.nz != pwtest2.nz|| pwtest2.nz != pwtest3.nz || pwtest3.nz != pwtest4.nz)
    {
        cout<<"Error"<<endl;
        exit(0);
    }

    int nrxx = pwtest1.nrxx;
    int npw = pwtest1.npw;
    int npw_half1 = pwtest3.npw;
    int npw_half2 = pwtest4.npw;

    double start,end;
    double t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    for(int it = 0 ; it < 10 ; ++it)
    {
        //init
        double * vr1 = new double[nrxx];
        double * vr2 = new double[nrxx];
        double * vr3 = new double[nrxx];
        double * vr4 = new double[nrxx];
        double * vpsi1 = new double[nrxx];
        double * vpsi2 = new double[nrxx];
        complex<double>* vpsic1 = new complex<double>[nrxx];
        complex<double>* vpsic2 = new complex<double>[nrxx];
        for(int i = 0 ; i < nrxx ; ++i) 
        {
            vr4[i] = vr3[i] = vr2[i] = vr1[i] = rand()/double(RAND_MAX);
        }

        complex<double> *psi = new complex<double> [npw*nbands];
        complex<double> *psi2 = new complex<double> [npw*nbands];
        complex<double> *psi3 = new complex<double> [npw_half1*nbands];
        complex<double> *psi4 = new complex<double> [npw_half2*nbands];
        complex<double> *psiout = new complex<double> [npw*nbands];
        complex<double> *psiout2 = new complex<double> [npw*nbands];
        complex<double> *psiout3 = new complex<double> [npw_half1*nbands];
        complex<double> *psiout4 = new complex<double> [npw_half2*nbands];
        for(int i = 0 ; i < npw*nbands ; ++i) 
        {
            psi2[i] = psi[i] = complex<double>(rand()/double(RAND_MAX), rand()/double(RAND_MAX));
        }

        for(int i = 0 ; i < npw_half1 * nbands ; ++i) 
        {
            psi3[i] = complex<double>(rand()/double(RAND_MAX), rand()/double(RAND_MAX));
        }
        for(int i = 0 ; i < npw_half2 * nbands ; ++i) 
        {
            psi4[i] = complex<double>(rand()/double(RAND_MAX), rand()/double(RAND_MAX));
        }


        MPI_Barrier(MPI_COMM_WORLD);
        start=MPI_Wtime();
        for(int i = 0 ; i < nbands ; ++i)
        {
            complex<double> * tmp = psi;
            complex<double> * tmpout = psiout;
            pwtest1.recip2real(tmp, vpsic1);
            for(int j = 0 ; j < nrxx ; ++j) vpsic1[j]*=vr1[j];
            pwtest1.real2recip(vpsic1,tmpout);
            tmp+=npw;
            tmpout+=npw;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        end=MPI_Wtime();
        t1 += end - start;

        MPI_Barrier(MPI_COMM_WORLD);
        start=MPI_Wtime();
        for(int i = 0 ; i < nbands ; ++i)
        {
            complex<double> * tmp = psi2;
            complex<double> * tmpout = psiout2;
            pwtest2.recip2real(tmp, vpsic2);
            for(int j = 0 ; j < nrxx ; ++j) vpsic2[j]*=vr2[j];
            pwtest2.real2recip(vpsic2,tmpout);
            tmp+=npw;
            tmpout+=npw;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        end=MPI_Wtime();
        t2 += end - start;

        MPI_Barrier(MPI_COMM_WORLD);
        start=MPI_Wtime();
        for(int i = 0 ; i < nbands ; ++i)
        {
            complex<double> * tmp = psi3;
            complex<double> * tmpout = psiout3;
            pwtest3.recip2real(tmp, vpsi1);
            for(int j = 0 ; j < nrxx ; ++j) vpsi1[j]*=vr3[j];
            pwtest3.real2recip(vpsi1,tmpout);
            tmp+=npw_half1;
            tmpout+=npw_half1;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        end=MPI_Wtime();
        t3 += end - start;

        MPI_Barrier(MPI_COMM_WORLD);
        start=MPI_Wtime();
        for(int i = 0 ; i < nbands ; ++i)
        {
            complex<double> * tmp = psi4;
            complex<double> * tmpout = psiout4;
            pwtest4.recip2real(tmp, vpsi2);
            for(int j = 0 ; j < nrxx ; ++j) vpsi2[j]*=vr4[j];
            pwtest4.real2recip(vpsi2,tmpout);
            tmp+=npw_half2;
            tmpout+=npw_half2;
        }
        MPI_Barrier(MPI_COMM_WORLD);
        end=MPI_Wtime();
        t4 += end - start;

        for(int i  = 0 ; i < nbands*npw; ++i)
        {
            double error = std::abs(psiout2[i]-psiout[i]);
            if(error > 1e-4) 
            {
                cout<<"Wrong"<<endl;
                exit(0);
            }
        }

        delete[] psi;
        delete[] psi2;
        delete[] psi3;
        delete[] psi4;
        delete[] psiout;
        delete[] psiout2;
        delete[] psiout3;
        delete[] psiout4;
        delete[] vr1;
        delete[] vr2;
        delete[] vr3;
        delete[] vr4;
        delete[] vpsic1;
        delete[] vpsic2;
        delete[] vpsi1;
        delete[] vpsi2;
    }

    if(myrank == 0)
    {
        cout<<setiosflags(ios::left);
	    cout<<setw(10)<<"time"<<setw(15)<<"xprime"<<setw(15)<<"yprime"<<endl;
	    cout<<setw(10)<<"fullnpw"<<setw(15)<<double(t1)<<setw(15)<<double(t2)<<"s"<<endl;
        cout<<setw(10)<<"halfnpw"<<setw(15)<<double(t3)<<setw(15)<<double(t4)<<"s"<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    

}