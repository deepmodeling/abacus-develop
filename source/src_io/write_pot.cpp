#include "../src_pw/potential.h"
#include "../src_pw/global.h"
#include "../module_surchem/efield.h"
#include "../module_base/timer.h"

// translate from write_rho in charge.cpp.
void Potential::write_potential(
	const int &is, 
	const int &iter, 
	const std::string &fn, 
	const ModuleBase::matrix &v, 
	const int &precision, 
	const int &hartree)const
{
    ModuleBase::TITLE("potential","write_potential");

    if(out_pot==0) 
    {
        return;
    }
    else if(out_pot<0)
    {
        if(hartree==0) return;
    }
    else if(iter % out_pot != 0)
    {
        return;
    }
    ModuleBase::timer::tick("potential","write_potential");

    std::ofstream ofs;

    if(GlobalV::MY_RANK==0)
    {
        ofs.open( fn.c_str() );

        ofs << GlobalC::ucell.latName << std::endl;//1
        ofs << " " << GlobalC::ucell.lat0 * 0.529177 << std::endl;
        ofs << " " << GlobalC::ucell.latvec.e11 << " " << GlobalC::ucell.latvec.e12 << " " << GlobalC::ucell.latvec.e13 << std::endl;
        ofs << " " << GlobalC::ucell.latvec.e21 << " " << GlobalC::ucell.latvec.e22 << " " << GlobalC::ucell.latvec.e23 << std::endl;
        ofs << " " << GlobalC::ucell.latvec.e31 << " " << GlobalC::ucell.latvec.e32 << " " << GlobalC::ucell.latvec.e33 << std::endl;

        for(int it=0; it<GlobalC::ucell.ntype; it++)
        {
            ofs << " " << GlobalC::ucell.atoms[it].label;
        }
        ofs << std::endl;
        for(int it=0; it<GlobalC::ucell.ntype; it++)
        {
            ofs << " " << GlobalC::ucell.atoms[it].na;
        }
        ofs << std::endl;
        ofs << "Direct" << std::endl;

        for(int it=0; it<GlobalC::ucell.ntype; it++)
        {
            for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
            {
                ofs << " " << GlobalC::ucell.atoms[it].taud[ia].x
                    << " " << GlobalC::ucell.atoms[it].taud[ia].y
                    << " " << GlobalC::ucell.atoms[it].taud[ia].z << std::endl;
            }
        }
        ofs << GlobalC::rhopw->nx << " " << GlobalC::rhopw->ny << " " << GlobalC::rhopw->nz;
        ofs << std::setprecision(precision);
        ofs << scientific; 
        if(!ofs)
        {
            ModuleBase::WARNING("potential::write_potential","Can't create VHartree File!");
        }
    }	

#ifndef __MPI
    int count=0;
    for(int k=0; k<GlobalC::rhopw->nz; k++)
    {
        ofs << "\n" << k << " iz";
        double value = 0.0;
        double ave = 0.0;
        for(int j=0; j<GlobalC::rhopw->ny; j++)
        {
            for(int i=0; i<GlobalC::rhopw->nx; i++)
            {
                if(count%8==0) ofs << "\n";
                value = v(is, i*GlobalC::rhopw->ny*GlobalC::rhopw->nz + j*GlobalC::rhopw->nz + k);
                ofs << " " << value;
                ave += value;
                ++count;
            }
        }
        ofs << "\n" << ave/GlobalC::rhopw->nx/GlobalC::rhopw->ny << " average";
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);
    // only do in the first pool.
    if(GlobalV::MY_POOL==0)
    {
        // num_z: how many planes on processor 'ip'
        int *num_z = new int[GlobalV::NPROC_IN_POOL];
        ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
        for (int iz=0;iz<GlobalC::rhopw->nz;iz++)
        {
            int ip = iz % GlobalV::NPROC_IN_POOL;
            num_z[ip]++;
        }

        // start_z: start position of z in
        // processor ip.
        int *start_z = new int[GlobalV::NPROC_IN_POOL];
        ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
        for (int ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
        {
            start_z[ip] = start_z[ip-1]+num_z[ip-1];
        }

        // which_ip: found iz belongs to which ip.
        int *which_ip = new int[GlobalC::rhopw->nz];
        ModuleBase::GlobalFunc::ZEROS(which_ip, GlobalC::rhopw->nz);
        for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
        {
            for(int ip=0; ip<GlobalV::NPROC_IN_POOL; ip++)
            {
                if(iz>=start_z[GlobalV::NPROC_IN_POOL-1])
                {
                    which_ip[iz] = GlobalV::NPROC_IN_POOL-1;
                    break;
                }
                else if(iz>=start_z[ip] && iz<start_z[ip+1])
                {
                    which_ip[iz] = ip;
                    break;
                }
            }
            //GlobalV::ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
        }
        int count=0;
        int nxy = GlobalC::rhopw->nx * GlobalC::rhopw->ny;
        double* zpiece = new double[nxy];
        // save the rho one z by one z.
        for(int iz=0; iz<GlobalC::rhopw->nz; iz++)
        {
            //GlobalV::ofs_running << "\n" << iz << " iz"; //LiuXh modify 20200624
            // tag must be different for different iz.
            ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
            int tag = iz;
            MPI_Status ierror;

            // case 1: the first part of rho in processor 0.
            if(which_ip[iz] == 0 && GlobalV::RANK_IN_POOL ==0)
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v(is, ir*GlobalC::rhopw->nplane+iz-GlobalC::rhopw->startz_current );
                    //GlobalV::ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*GlobalC::rhopw->nplane+iz=" << ir*GlobalC::rhopw->nplane+iz;
                }
            }
            // case 2: > first part rho: send the rho to
            // processor 0.
            else if(which_ip[iz] == GlobalV::RANK_IN_POOL )
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v(is, ir*GlobalC::rhopw->nplane+iz-GlobalC::rhopw->startz_current);
                }
                MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
            }

            // case 2: > first part rho: processor 0 receive the rho
            // from other processors
            else if(GlobalV::RANK_IN_POOL==0)
            {
                MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
                //GlobalV::ofs_running << "\n Receieve First number = " << zpiece[0];
            }

            // write data
            if(GlobalV::MY_RANK==0)
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
    if(GlobalV::MY_RANK==0) ofs.close();
    ModuleBase::timer::tick("potential","write_potential");
    return;
}


void Potential::write_elecstat_pot(const std::string &fn, const std::string &fn_ave, ModulePW::PW_Basis* rho_basis)
{
    ModuleBase::TITLE("Potential","write_elecstat_pot");
    ModuleBase::timer::tick("Potential","write_elecstat_pot");

    double *v_elecstat = new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(v_elecstat, rho_basis->nrxx);

    std::complex<double> *vh_r = new std::complex<double>[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS( vh_r, rho_basis->nrxx );
    std::complex<double> *vh_g  = new std::complex<double>[rho_basis->npw];

    
    int nspin0 = 1;
    if(GlobalV::NSPIN==2) nspin0 = GlobalV::NSPIN;
    for(int is=0; is<nspin0; is++)
    {
        for(int ir=0; ir<rho_basis->nrxx; ir++)
        {
            vh_r[ir] += std::complex<double>( GlobalC::CHR.rho[is][ir], 0.0 );
        }
    }

    //=============================
    //  bring rho (aux) to G space
    //=============================
    rho_basis->real2recip(vh_r, vh_g);

    //=======================================================
    // calculate hartree potential in G-space (NB: V(G=0)=0 )
    //=======================================================
    

    for(int ig = 0; ig < rho_basis->npw; ++ig)
    {
        if(rho_basis->ig_gge0==ig)    continue;
        const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (GlobalC::ucell.tpiba2 * rho_basis->gg[ig]);
        vh_g[ig] *= fac;
    }

    //==========================================
    //transform hartree potential to real space
    //==========================================
    rho_basis->recip2real(vh_g, vh_r);

    //==========================================
    // Dipole correction
    //==========================================
    ModuleBase::matrix v_efield;
    if (GlobalV::EFIELD_FLAG && GlobalV::DIP_COR_FLAG)
    {
        v_efield.create(GlobalV::NSPIN, rho_basis->nrxx);
        v_efield = Efield::add_efield(GlobalC::ucell, GlobalC::rhopw, GlobalV::NSPIN, GlobalC::CHR.rho, GlobalC::solvent_model);
    }

    //==========================================
    //Add hartree potential and local pseudopot
    //==========================================
    for (int ir = 0;ir < rho_basis->nrxx;ir++)
    {
        v_elecstat[ir] = vh_r[ir].real() + this->vltot[ir];

        if (GlobalV::EFIELD_FLAG && GlobalV::DIP_COR_FLAG)
        {
            v_elecstat[ir] += v_efield(0, ir);
        }
    }

    //-------------------------------------------
    // output the electrostatic potential into a file.
    //-------------------------------------------
    std::ofstream ofs;
    std::ofstream ofs_ave;

    if(GlobalV::MY_RANK==0)
    {
        ofs.open( fn.c_str() );
        ofs_ave.open( fn_ave.c_str() );

        ofs << GlobalC::ucell.latName << std::endl;//1
        ofs << " " << GlobalC::ucell.lat0 * 0.529177 << std::endl;
        ofs << " " << GlobalC::ucell.latvec.e11 << " " << GlobalC::ucell.latvec.e12 << " " << GlobalC::ucell.latvec.e13 << std::endl;
        ofs << " " << GlobalC::ucell.latvec.e21 << " " << GlobalC::ucell.latvec.e22 << " " << GlobalC::ucell.latvec.e23 << std::endl;
        ofs << " " << GlobalC::ucell.latvec.e31 << " " << GlobalC::ucell.latvec.e32 << " " << GlobalC::ucell.latvec.e33 << std::endl;

        ofs_ave << GlobalC::ucell.latName << std::endl;//1
        ofs_ave << " " << GlobalC::ucell.lat0 * 0.529177 << std::endl;
        ofs_ave << " " << GlobalC::ucell.latvec.e11 << " " << GlobalC::ucell.latvec.e12 << " " << GlobalC::ucell.latvec.e13 << std::endl;
        ofs_ave << " " << GlobalC::ucell.latvec.e21 << " " << GlobalC::ucell.latvec.e22 << " " << GlobalC::ucell.latvec.e23 << std::endl;
        ofs_ave << " " << GlobalC::ucell.latvec.e31 << " " << GlobalC::ucell.latvec.e32 << " " << GlobalC::ucell.latvec.e33 << std::endl;

        for(int it=0; it<GlobalC::ucell.ntype; it++)
        {
            ofs << " " << GlobalC::ucell.atoms[it].label;
            ofs_ave << " " << GlobalC::ucell.atoms[it].label;
        }
        ofs << std::endl;
        ofs_ave << std::endl;
        for(int it=0; it<GlobalC::ucell.ntype; it++)
        {
            ofs << " " << GlobalC::ucell.atoms[it].na;
            ofs_ave << " " << GlobalC::ucell.atoms[it].na;
        }
        ofs << std::endl;
        ofs << "Direct" << std::endl;

        ofs_ave << std::endl;
        ofs_ave << "Direct" << std::endl;

        for(int it=0; it<GlobalC::ucell.ntype; it++)
        {
            for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
            {
                ofs << " " << GlobalC::ucell.atoms[it].taud[ia].x
                    << " " << GlobalC::ucell.atoms[it].taud[ia].y
                    << " " << GlobalC::ucell.atoms[it].taud[ia].z << std::endl;

                ofs_ave << " " << GlobalC::ucell.atoms[it].taud[ia].x
                    << " " << GlobalC::ucell.atoms[it].taud[ia].y
                    << " " << GlobalC::ucell.atoms[it].taud[ia].z << std::endl;
            }
        }

        ofs << rho_basis->nx << " " << rho_basis->ny << " " << rho_basis->nz;
        ofs_ave << rho_basis->nx << " " << rho_basis->ny << " " << rho_basis->nz;

        int precision = 9;
        ofs << std::setprecision(precision);
        ofs << scientific; 
        ofs_ave << std::setprecision(precision);
        ofs_ave << scientific; 
        if(!ofs)
        {
            ModuleBase::WARNING("potential::write_potential","Can't create VHartree File!");
        }
    }	

#ifndef __MPI
    int count=0;
    for(int k=0; k<rho_basis->nz; k++)
    {
        ofs << "\n" << k << " iz";
        double value = 0.0;
        double ave = 0.0;
        for(int j=0; j<rho_basis->ny; j++)
        {
            for(int i=0; i<rho_basis->nx; i++)
            {
                //if(count%8==0) ofs << "\n";
                if(count%5==0) ofs << "\n";
                value = v_elecstat[i*rho_basis->ny*rho_basis->nz + j*rho_basis->nz + k];
                ofs << " " << value;
                ave += value;
                ++count;
            }
        }
        //ofs << "\n" << ave/GlobalC::rhopw->nx/GlobalC::rhopw->ny << " average";
        if(k==0) ofs_ave << "iz" << "\taverage";
        ofs_ave << "\n" << k << "\t" << ave/rho_basis->nx/rho_basis->ny;
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);
    // only do in the first pool.
    if(GlobalV::MY_POOL==0)
    {
        // num_z: how many planes on processor 'ip'
        int *num_z = new int[GlobalV::NPROC_IN_POOL];
        ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
        //for (int iz=0;iz<GlobalC::rhopw->nz;iz++)
        //{
        //    int ip = iz % GlobalV::NPROC_IN_POOL;
        //    num_z[ip]++;
        //}
        for (int iz=0;iz<GlobalC::bigpw->nbz;iz++)
        {
            int ip = iz % GlobalV::NPROC_IN_POOL;
            num_z[ip] += GlobalC::bigpw->bz;
        }

        // start_z: start position of z in
        // processor ip.
        int *start_z = new int[GlobalV::NPROC_IN_POOL];
        ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
        for (int ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
        {
            start_z[ip] = start_z[ip-1]+num_z[ip-1];
        }

        // which_ip: found iz belongs to which ip.
        int *which_ip = new int[rho_basis->nz];
        ModuleBase::GlobalFunc::ZEROS(which_ip, rho_basis->nz);
        for(int iz=0; iz<rho_basis->nz; iz++)
        {
            for(int ip=0; ip<GlobalV::NPROC_IN_POOL; ip++)
            {
                if(iz>=start_z[GlobalV::NPROC_IN_POOL-1])
                {
                    which_ip[iz] = GlobalV::NPROC_IN_POOL-1;
                    break;
                }
                else if(iz>=start_z[ip] && iz<start_z[ip+1])
                {
                    which_ip[iz] = ip;
                    break;
                }
            }
            //GlobalV::ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
        }
        int count=0;
        int nxy = rho_basis->nxy;
        double* zpiece = new double[nxy];
        // save the rho one z by one z.
        for(int iz=0; iz<rho_basis->nz; iz++)
        {
            //GlobalV::ofs_running << "\n" << iz << " iz";
            // tag must be different for different iz.
            ModuleBase::GlobalFunc::ZEROS(zpiece, nxy);
            int tag = iz;
            MPI_Status ierror;

            // case 1: the first part of rho in processor 0.
            if(which_ip[iz] == 0 && GlobalV::RANK_IN_POOL ==0)
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v_elecstat[ir*GlobalC::rhopw->nplane+iz-GlobalC::rhopw->startz_current ];
                    //GlobalV::ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*GlobalC::rhopw->nplane+iz=" << ir*GlobalC::rhopw->nplane+iz;
                }
            }
            // case 2: > first part rho: send the rho to
            // processor 0.
            else if(which_ip[iz] == GlobalV::RANK_IN_POOL )
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v_elecstat[ir*GlobalC::rhopw->nplane+iz-GlobalC::rhopw->startz_current];
                }
                MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
            }

            // case 2: > first part rho: processor 0 receive the rho
            // from other processors
            else if(GlobalV::RANK_IN_POOL==0)
            {
                MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
                //GlobalV::ofs_running << "\n Receieve First number = " << zpiece[0];
            }

            // write data
            if(GlobalV::MY_RANK==0)
            {
                //ofs << "\niz=" << iz;
                double ave = 0.0;
                /*
                for(int ir=0; ir<nxy; ir++)
                {
                    //if(count%8==0) ofs << "\n";
                    if(count%5==0) ofs << "\n";
                    //ofs << " " << zpiece[ir];
                    ofs << std::setw(17) << zpiece[ir];
                    ave += zpiece[ir];
                    ++count;
                }
                //ofs << "\n" << ave/nxy << " average"; 
                */
                for(int iy=0; iy<rho_basis->ny; iy++)
                {
                    for(int ix=0; ix<rho_basis->nx; ix++)
                    {
                        //if(count%8==0) ofs << "\n";
                        if(count%5==0) ofs << "\n";
                        //ofs << " " << zpiece[ir];
                        ofs << std::setw(17) << zpiece[ix*rho_basis->ny+iy];
                        ave += zpiece[ix*rho_basis->ny+iy];
                        ++count;
                    }
                }
                if(iz==0) ofs_ave << "\niz" << "\taverage";
                ofs_ave << "\n" << iz << "\t" << ave/rho_basis->nx/rho_basis->ny;
            }
        }
        delete[] num_z;
        delete[] start_z;
        delete[] which_ip;
        delete[] zpiece;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(GlobalV::MY_RANK==0) ofs.close();

    delete[] v_elecstat;
    delete[] vh_g;
    delete[] vh_r;

    ModuleBase::timer::tick("Potential","write_potential");
    return;
}
