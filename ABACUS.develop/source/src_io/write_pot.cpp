#include "../src_pw/potential.h"
#include "../src_pw/global.h"

// translate from write_rho in charge.cpp.
void Potential::write_potential(
	const int &is, 
	const int &iter, 
	const string &fn, 
	const matrix &v, 
	const int &precision, 
	const int &hartree)const
{
    TITLE("potential","write_potential");

    if(out_potential==0) 
    {
        return;
    }
    else if(out_potential<0)
    {
        if(hartree==0) return;
    }
    else if(iter % out_potential != 0)
    {
        return;
    }
    timer::tick("potential","write_potential");

    ofstream ofs;

    if(MY_RANK==0)
    {
        ofs.open( fn.c_str() );

        ofs << ucell.latName << endl;//1
        ofs << " " << ucell.lat0 * 0.529177 << endl;
        ofs << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
        ofs << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
        ofs << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

        for(int it=0; it<ucell.ntype; it++)
        {
            ofs << " " << ucell.atoms[it].label;
        }
        ofs << endl;
        for(int it=0; it<ucell.ntype; it++)
        {
            ofs << " " << ucell.atoms[it].na;
        }
        ofs << endl;
        ofs << "Direct" << endl;

        for(int it=0; it<ucell.ntype; it++)
        {
            for(int ia=0; ia<ucell.atoms[it].na; ia++)
            {
                ofs << " " << ucell.atoms[it].taud[ia].x
                    << " " << ucell.atoms[it].taud[ia].y
                    << " " << ucell.atoms[it].taud[ia].z << endl;
            }
        }
        ofs << pw.ncx << " " << pw.ncy << " " << pw.ncz;
        ofs << setprecision(precision);
        ofs << scientific; 
        if(!ofs)
        {
            WARNING("potential::write_potential","Can't create VHartree File!");
        }
    }	

#ifndef __MPI
    int count=0;
    for(int k=0; k<pw.ncz; k++)
    {
        ofs << "\n" << k << " iz";
        double value = 0.0;
        double ave = 0.0;
        for(int j=0; j<pw.ncy; j++)
        {
            for(int i=0; i<pw.ncx; i++)
            {
                if(count%8==0) ofs << "\n";
                value = v(is, i*pw.ncy*pw.ncz + j*pw.ncz + k);
                ofs << " " << value;
                ave += value;
                ++count;
            }
        }
        ofs << "\n" << ave/pw.ncx/pw.ncy << " average";
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);
    // only do in the first pool.
    if(MY_POOL==0)
    {
        // num_z: how many planes on processor 'ip'
        int *num_z = new int[NPROC_IN_POOL];
        ZEROS(num_z, NPROC_IN_POOL);
        for (int iz=0;iz<pw.ncz;iz++)
        {
            int ip = iz % NPROC_IN_POOL;
            num_z[ip]++;
        }

        // start_z: start position of z in
        // processor ip.
        int *start_z = new int[NPROC_IN_POOL];
        ZEROS(start_z, NPROC_IN_POOL);
        for (int ip=1;ip<NPROC_IN_POOL;ip++)
        {
            start_z[ip] = start_z[ip-1]+num_z[ip-1];
        }

        // which_ip: found iz belongs to which ip.
        int *which_ip = new int[pw.ncz];
        ZEROS(which_ip, pw.ncz);
        for(int iz=0; iz<pw.ncz; iz++)
        {
            for(int ip=0; ip<NPROC_IN_POOL; ip++)
            {
                if(iz>=start_z[NPROC_IN_POOL-1])
                {
                    which_ip[iz] = NPROC_IN_POOL-1;
                    break;
                }
                else if(iz>=start_z[ip] && iz<start_z[ip+1])
                {
                    which_ip[iz] = ip;
                    break;
                }
            }
            //ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
        }
        int count=0;
        int nxy = pw.ncx * pw.ncy;
        double* zpiece = new double[nxy];
        // save the rho one z by one z.
        for(int iz=0; iz<pw.ncz; iz++)
        {
            //ofs_running << "\n" << iz << " iz"; //LiuXh modify 20200624
            // tag must be different for different iz.
            ZEROS(zpiece, nxy);
            int tag = iz;
            MPI_Status ierror;

            // case 1: the first part of rho in processor 0.
            if(which_ip[iz] == 0 && RANK_IN_POOL ==0)
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v(is, ir*pw.nczp+iz-start_z[RANK_IN_POOL] );
                    //ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
                }
            }
            // case 2: > first part rho: send the rho to
            // processor 0.
            else if(which_ip[iz] == RANK_IN_POOL )
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v(is, ir*pw.nczp+iz-start_z[RANK_IN_POOL]);
                }
                MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
            }

            // case 2: > first part rho: processor 0 receive the rho
            // from other processors
            else if(RANK_IN_POOL==0)
            {
                MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
                //ofs_running << "\n Receieve First number = " << zpiece[0];
            }

            // write data
            if(MY_RANK==0)
            {
                //ofs << "\niz=" << iz;
                double ave = 0.0;
                for(int ir=0; ir<nxy; ir++)
                {
                    if(count%8==0) ofs << "\n";
                    ofs << " " << zpiece[ir];
                    ave += zpiece[ir];
                    ++count;
                }
                ofs << "\n" << ave/nxy << " average"; 
            }
        }
        delete[] zpiece;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(MY_RANK==0) ofs.close();
    timer::tick("potential","write_potential");
    return;
}


void Potential::write_elecstat_pot(const string &fn, const string &fn_ave)
{
    TITLE("Potential","write_elecstat_pot");
    timer::tick("Potential","write_elecstat_pot");

    double *v_elecstat;
    v_elecstat = new double[pw.nrxx];
    ZEROS(v_elecstat, pw.nrxx);

    complex<double> *Porter = UFFT.porter;
    ZEROS( Porter, pw.nrxx );
    
    int nspin0 = 1;
    if(NSPIN==2) nspin0 = NSPIN;
    for(int is=0; is<nspin0; is++)
    {
        for(int ir=0; ir<pw.nrxx; ir++)
        {
            Porter[ir] += complex<double>( CHR.rho[is][ir], 0.0 );
        }
    }

    //=============================
    //  bring rho (aux) to G space
    //=============================
    pw.FFT_chg.FFT3D(Porter, -1);

    //=======================================================
    // calculate hartree potential in G-space (NB: V(G=0)=0 )
    //=======================================================
    complex<double> *vh_g  = new complex<double>[pw.ngmc];
    ZEROS(vh_g, pw.ngmc);

    for(int ig = pw.gstart; ig<pw.ngmc; ig++)
    {
        const int j = pw.ig2fftc[ig];
        if(pw.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            const double fac = e2 * FOUR_PI / (ucell.tpiba2 * pw.gg [ig]);
            vh_g[ig] = fac * Porter[j];
        }
    }

    ZEROS(Porter, pw.nrxx);

    for (int ig = 0;ig < pw.ngmc;ig++)
    {
        Porter[pw.ig2fftc[ig]] = vh_g[ig];
    }

    //==========================================
    //transform hartree potential to real space
    //==========================================
    pw.FFT_chg.FFT3D(Porter, 1);
    //==========================================
    //Add hartree potential and local pseudopot
    //==========================================
    for (int ir = 0;ir < pw.nrxx;ir++)
    {
        v_elecstat[ir] = Porter[ir].real() + this->vltot[ir];
    }

    //-------------------------------------------
    // output the electrostatic potential into a file.
    //-------------------------------------------
    ofstream ofs;
    ofstream ofs_ave;

    if(MY_RANK==0)
    {
        ofs.open( fn.c_str() );
        ofs_ave.open( fn_ave.c_str() );

        ofs << ucell.latName << endl;//1
        ofs << " " << ucell.lat0 * 0.529177 << endl;
        ofs << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
        ofs << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
        ofs << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

        ofs_ave << ucell.latName << endl;//1
        ofs_ave << " " << ucell.lat0 * 0.529177 << endl;
        ofs_ave << " " << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
        ofs_ave << " " << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
        ofs_ave << " " << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

        for(int it=0; it<ucell.ntype; it++)
        {
            ofs << " " << ucell.atoms[it].label;
            ofs_ave << " " << ucell.atoms[it].label;
        }
        ofs << endl;
        ofs_ave << endl;
        for(int it=0; it<ucell.ntype; it++)
        {
            ofs << " " << ucell.atoms[it].na;
            ofs_ave << " " << ucell.atoms[it].na;
        }
        ofs << endl;
        ofs << "Direct" << endl;

        ofs_ave << endl;
        ofs_ave << "Direct" << endl;

        for(int it=0; it<ucell.ntype; it++)
        {
            for(int ia=0; ia<ucell.atoms[it].na; ia++)
            {
                ofs << " " << ucell.atoms[it].taud[ia].x
                    << " " << ucell.atoms[it].taud[ia].y
                    << " " << ucell.atoms[it].taud[ia].z << endl;

                ofs_ave << " " << ucell.atoms[it].taud[ia].x
                    << " " << ucell.atoms[it].taud[ia].y
                    << " " << ucell.atoms[it].taud[ia].z << endl;
            }
        }

        ofs << pw.ncx << " " << pw.ncy << " " << pw.ncz;
        ofs_ave << pw.ncx << " " << pw.ncy << " " << pw.ncz;

        int precision = 9;
        ofs << setprecision(precision);
        ofs << scientific; 
        ofs_ave << setprecision(precision);
        ofs_ave << scientific; 
        if(!ofs)
        {
            WARNING("potential::write_potential","Can't create VHartree File!");
        }
    }	

#ifndef __MPI
    int count=0;
    for(int k=0; k<pw.ncz; k++)
    {
        ofs << "\n" << k << " iz";
        double value = 0.0;
        double ave = 0.0;
        for(int j=0; j<pw.ncy; j++)
        {
            for(int i=0; i<pw.ncx; i++)
            {
                //if(count%8==0) ofs << "\n";
                if(count%5==0) ofs << "\n";
                value = v_elecstat[i*pw.ncy*pw.ncz + j*pw.ncz + k];
                ofs << " " << value;
                ave += value;
                ++count;
            }
        }
        //ofs << "\n" << ave/pw.ncx/pw.ncy << " average";
        if(k==0) ofs_ave << "iz" << "\taverage";
        ofs_ave << "\n" << k << "\t" << ave/pw.ncx/pw.ncy;
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);
    // only do in the first pool.
    if(MY_POOL==0)
    {
        // num_z: how many planes on processor 'ip'
        int *num_z = new int[NPROC_IN_POOL];
        ZEROS(num_z, NPROC_IN_POOL);
        //for (int iz=0;iz<pw.ncz;iz++)
        //{
        //    int ip = iz % NPROC_IN_POOL;
        //    num_z[ip]++;
        //}
        for (int iz=0;iz<pw.nbz;iz++)
        {
            int ip = iz % NPROC_IN_POOL;
            num_z[ip] += pw.bz;
        }

        // start_z: start position of z in
        // processor ip.
        int *start_z = new int[NPROC_IN_POOL];
        ZEROS(start_z, NPROC_IN_POOL);
        for (int ip=1;ip<NPROC_IN_POOL;ip++)
        {
            start_z[ip] = start_z[ip-1]+num_z[ip-1];
        }

        // which_ip: found iz belongs to which ip.
        int *which_ip = new int[pw.ncz];
        ZEROS(which_ip, pw.ncz);
        for(int iz=0; iz<pw.ncz; iz++)
        {
            for(int ip=0; ip<NPROC_IN_POOL; ip++)
            {
                if(iz>=start_z[NPROC_IN_POOL-1])
                {
                    which_ip[iz] = NPROC_IN_POOL-1;
                    break;
                }
                else if(iz>=start_z[ip] && iz<start_z[ip+1])
                {
                    which_ip[iz] = ip;
                    break;
                }
            }
            //ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
        }
        int count=0;
        int nxy = pw.ncx * pw.ncy;
        double* zpiece = new double[nxy];
        // save the rho one z by one z.
        for(int iz=0; iz<pw.ncz; iz++)
        {
            //ofs_running << "\n" << iz << " iz";
            // tag must be different for different iz.
            ZEROS(zpiece, nxy);
            int tag = iz;
            MPI_Status ierror;

            // case 1: the first part of rho in processor 0.
            if(which_ip[iz] == 0 && RANK_IN_POOL ==0)
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v_elecstat[ir*pw.nczp+iz-start_z[RANK_IN_POOL] ];
                    //ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
                }
            }
            // case 2: > first part rho: send the rho to
            // processor 0.
            else if(which_ip[iz] == RANK_IN_POOL )
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v_elecstat[ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
                }
                MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
            }

            // case 2: > first part rho: processor 0 receive the rho
            // from other processors
            else if(RANK_IN_POOL==0)
            {
                MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
                //ofs_running << "\n Receieve First number = " << zpiece[0];
            }

            // write data
            if(MY_RANK==0)
            {
                //ofs << "\niz=" << iz;
                double ave = 0.0;
                /*
                for(int ir=0; ir<nxy; ir++)
                {
                    //if(count%8==0) ofs << "\n";
                    if(count%5==0) ofs << "\n";
                    //ofs << " " << zpiece[ir];
                    ofs << setw(17) << zpiece[ir];
                    ave += zpiece[ir];
                    ++count;
                }
                //ofs << "\n" << ave/nxy << " average"; 
                */
                for(int iy=0; iy<pw.ncy; iy++)
                {
                    for(int ix=0; ix<pw.ncx; ix++)
                    {
                        //if(count%8==0) ofs << "\n";
                        if(count%5==0) ofs << "\n";
                        //ofs << " " << zpiece[ir];
                        ofs << setw(17) << zpiece[ix*pw.ncy+iy];
                        ave += zpiece[ix*pw.ncy+iy];
                        ++count;
                    }
                }
                if(iz==0) ofs_ave << "\niz" << "\taverage";
                ofs_ave << "\n" << iz << "\t" << ave/pw.ncx/pw.ncy;
            }
        }
        delete[] num_z;
        delete[] start_z;
        delete[] which_ip;
        delete[] zpiece;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(MY_RANK==0) ofs.close();

    delete[] v_elecstat;
    delete[] vh_g;

    timer::tick("Potential","write_potential");
    return;
}
