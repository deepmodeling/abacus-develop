#include "potential_new.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/potentials/efield.h"
#include "module_base/timer.h"
#include "module_base/element_name.h"
namespace elecstate
{
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
	if (GlobalV::out_pot!=1) 
	{
		return;
	}
    ModuleBase::timer::tick("Potential","write_potential");

    std::ofstream ofs;

	if(GlobalV::MY_RANK==0)
	{
		ofs.open(fn.c_str());
		
		if (!ofs)
		{
            ModuleBase::WARNING("Potential::write_potential","Can't create Potential File!");
		}	

		/// output header for cube file
		ofs << "Cubefile created from ABACUS SCF calculation. The inner loop is z index, followed by y index, x index in turn." << std::endl;
		ofs << GlobalV::NSPIN << " (nspin) " << std::endl;

        // We assum that SPIN1_POT.cube does not need Fermi Energy.
		// if(GlobalV::NSPIN==1 || GlobalV::NSPIN == 4)
		// {
		// 	ofs << GlobalC::en.ef << " (fermi energy, in Ry)" << std::endl;
		// }
		// else if(GlobalV::NSPIN==2)
		// {
		// 	if (GlobalV::TWO_EFERMI)
		// 	{
		// 		if(is==0)		ofs << GlobalC::en.ef_up << " (fermi energy for spin=1, in Ry)" << std::endl; 
		// 		else if(is==1)	ofs << GlobalC::en.ef_dw << " (fermi energy for spin=2, in Ry)" << std::endl;
		// 	}
		// 	else
		// 	{
		// 		ofs << GlobalC::en.ef << " (fermi energy, in Ry)" << std::endl;
		// 	}
		// }
		// else
		// {
		// 	ModuleBase::WARNING_QUIT("write_potential","check nspin!");
		// }

		ofs << GlobalC::ucell.nat << " 0.0 0.0 0.0 " << std::endl;
		double fac=GlobalC::ucell.lat0;
		ofs << this->rho_basis_->nx 
			<< " " << fac*GlobalC::ucell.latvec.e11/double(this->rho_basis_->nx) 
			<< " " << fac*GlobalC::ucell.latvec.e12/double(this->rho_basis_->nx) 
			<< " " << fac*GlobalC::ucell.latvec.e13/double(this->rho_basis_->nx) << std::endl;
		ofs << this->rho_basis_->ny 
			<< " " << fac*GlobalC::ucell.latvec.e21/double(this->rho_basis_->ny) 
			<< " " << fac*GlobalC::ucell.latvec.e22/double(this->rho_basis_->ny) 
			<< " " << fac*GlobalC::ucell.latvec.e23/double(this->rho_basis_->ny) << std::endl;
		ofs << this->rho_basis_->nz 
			<< " " << fac*GlobalC::ucell.latvec.e31/double(this->rho_basis_->nz) 
			<< " " << fac*GlobalC::ucell.latvec.e32/double(this->rho_basis_->nz) 
			<< " " << fac*GlobalC::ucell.latvec.e33/double(this->rho_basis_->nz) << std::endl;

		std::string element = "";
		for(int it=0; it<GlobalC::ucell.ntype; it++)
		{
			// erase the number in label, such as Fe1.
			element = GlobalC::ucell.atoms[it].label;
			std::string::iterator temp = element.begin();
			while (temp != element.end())
			{
				if ((*temp >= '1') && (*temp <= '9'))
				{
					temp = element.erase(temp);
				}
				else
				{
					temp++;
				}
			}

			for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				//convert from label to atomic number
				int z = 0;
				for(int j=0; j!=ModuleBase::element_name.size(); j++)
				{
					if (element == ModuleBase::element_name[j])
					{
						z=j+1;
						break;
					}
				}
				ofs << " " << z << " " << GlobalC::ucell.atoms[it].ncpp.zv
						 << " " << fac*GlobalC::ucell.atoms[it].tau[ia].x
						 << " " << fac*GlobalC::ucell.atoms[it].tau[ia].y
						 << " " << fac*GlobalC::ucell.atoms[it].tau[ia].z << std::endl;
			}
		}
		ofs << std::setprecision(precision);
		ofs << scientific;
	}

#ifndef __MPI
    int count=0;
	for(int i=0; i<this->rho_basis_->nx; i++)
	{
		for(int j=0; j<this->rho_basis_->ny; j++)
		{
			for(int k=0; k<this->rho_basis_->nz; k++)
			{
                ofs << " " << v(is, k*this->rho_basis_->nx*this->rho_basis_->ny+i*this->rho_basis_->ny+j);
				if(k%6==5 && k!=this->rho_basis_->nz-1) ofs << "\n";
            }
            ofs << "\n";
        }
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);
    // only do in the first pool.
    if(GlobalV::MY_POOL==0)
    {
        int nxyz = this->rho_basis_->nx * this->rho_basis_->ny * this->rho_basis_->nz;
		double* pot_cube = new double[nxyz];
		ModuleBase::GlobalFunc::ZEROS(pot_cube, nxyz);

        // num_z: how many planes on processor 'ip'
        int *num_z = new int[GlobalV::NPROC_IN_POOL];
        ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
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
        int *which_ip = new int[this->rho_basis_->nz];
        ModuleBase::GlobalFunc::ZEROS(which_ip, this->rho_basis_->nz);
        for(int iz=0; iz<this->rho_basis_->nz; iz++)
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
        int nxy = this->rho_basis_->nx * this->rho_basis_->ny;
        double* zpiece = new double[nxy];
        // save the rho one z by one z.
        for(int iz=0; iz<this->rho_basis_->nz; iz++)
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
                    zpiece[ir] = v(is, ir*this->rho_basis_->nplane+iz-this->rho_basis_->startz_current );
                }
            }
            // case 2: > first part rho: send the rho to
            // processor 0.
            else if(which_ip[iz] == GlobalV::RANK_IN_POOL )
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v(is, ir*this->rho_basis_->nplane+iz-this->rho_basis_->startz_current);
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

            if(GlobalV::MY_RANK==0)
            {
                double ave = 0.0;
                for(int ir=0; ir<nxy; ir++)
                {
                    pot_cube[ir+iz*nxy]=zpiece[ir];
                }
            }
        }
        delete[] zpiece;
		delete[] which_ip;
		delete[] num_z;
		delete[] start_z;
        if(GlobalV::MY_RANK==0)
		{
			for(int ix=0; ix<this->rho_basis_->nx; ix++)
			{
				for(int iy=0; iy<this->rho_basis_->ny; iy++)
				{
					for (int iz=0; iz<this->rho_basis_->nz; iz++)
					{
						ofs << " " << pot_cube[iz*this->rho_basis_->nx*this->rho_basis_->ny+ix*this->rho_basis_->ny+iy];
						if(iz%6==5 && iz!=this->rho_basis_->nz-1) ofs << "\n";
					}
					ofs << "\n";
				}
			}
		}
        delete[] pot_cube;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if(GlobalV::MY_RANK==0) ofs.close();
    ModuleBase::timer::tick("Potential","write_potential");
    return;
}


void Potential::write_elecstat_pot(const std::string &fn, const std::string &fn_ave, ModulePW::PW_Basis* rho_basis, const Charge* const chr)
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
            vh_r[ir] += std::complex<double>( chr->rho[is][ir], 0.0 );
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
        const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (this->ucell_->tpiba2 * rho_basis->gg[ig]);
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
        v_efield = elecstate::Efield::add_efield(*(this->ucell_), const_cast<ModulePW::PW_Basis *>(this->rho_basis_), GlobalV::NSPIN, chr->rho, GlobalC::solvent_model);
    }

    //==========================================
    //Add hartree potential and local pseudopot
    //==========================================
    for (int ir = 0;ir < rho_basis->nrxx;ir++)
    {
        v_elecstat[ir] = vh_r[ir].real() + this->v_effective_fixed[ir];

        if (GlobalV::EFIELD_FLAG && GlobalV::DIP_COR_FLAG)
        {
            v_elecstat[ir] += v_efield(0, ir);
        }
        if (GlobalV::imp_sol)
        {
            v_elecstat[ir] += GlobalC::solvent_model.delta_phi[ir];
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

        ofs << this->ucell_->latName << std::endl;//1
        ofs << " " << this->ucell_->lat0 * 0.529177 << std::endl;
        ofs << " " << this->ucell_->latvec.e11 << " " << this->ucell_->latvec.e12 << " " << this->ucell_->latvec.e13 << std::endl;
        ofs << " " << this->ucell_->latvec.e21 << " " << this->ucell_->latvec.e22 << " " << this->ucell_->latvec.e23 << std::endl;
        ofs << " " << this->ucell_->latvec.e31 << " " << this->ucell_->latvec.e32 << " " << this->ucell_->latvec.e33 << std::endl;

        ofs_ave << this->ucell_->latName << std::endl;//1
        ofs_ave << " " << this->ucell_->lat0 * 0.529177 << std::endl;
        ofs_ave << " " << this->ucell_->latvec.e11 << " " << this->ucell_->latvec.e12 << " " << this->ucell_->latvec.e13 << std::endl;
        ofs_ave << " " << this->ucell_->latvec.e21 << " " << this->ucell_->latvec.e22 << " " << this->ucell_->latvec.e23 << std::endl;
        ofs_ave << " " << this->ucell_->latvec.e31 << " " << this->ucell_->latvec.e32 << " " << this->ucell_->latvec.e33 << std::endl;

        for(int it=0; it<this->ucell_->ntype; it++)
        {
            ofs << " " << this->ucell_->atoms[it].label;
            ofs_ave << " " << this->ucell_->atoms[it].label;
        }
        ofs << std::endl;
        ofs_ave << std::endl;
        for(int it=0; it<this->ucell_->ntype; it++)
        {
            ofs << " " << this->ucell_->atoms[it].na;
            ofs_ave << " " << this->ucell_->atoms[it].na;
        }
        ofs << std::endl;
        ofs << "Direct" << std::endl;

        ofs_ave << std::endl;
        ofs_ave << "Direct" << std::endl;

        for(int it=0; it<this->ucell_->ntype; it++)
        {
            for(int ia=0; ia<this->ucell_->atoms[it].na; ia++)
            {
                ofs << " " << this->ucell_->atoms[it].taud[ia].x
                    << " " << this->ucell_->atoms[it].taud[ia].y
                    << " " << this->ucell_->atoms[it].taud[ia].z << std::endl;

                ofs_ave << " " << this->ucell_->atoms[it].taud[ia].x
                    << " " << this->ucell_->atoms[it].taud[ia].y
                    << " " << this->ucell_->atoms[it].taud[ia].z << std::endl;
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
        //ofs << "\n" << ave/this->rho_basis_->nx/this->rho_basis_->ny << " average";
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
        //for (int iz=0;iz<this->rho_basis_->nz;iz++)
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
                    zpiece[ir] = v_elecstat[ir*this->rho_basis_->nplane+iz-this->rho_basis_->startz_current ];
                    //GlobalV::ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*this->rho_basis_->nplane+iz=" << ir*this->rho_basis_->nplane+iz;
                }
            }
            // case 2: > first part rho: send the rho to
            // processor 0.
            else if(which_ip[iz] == GlobalV::RANK_IN_POOL )
            {
                for(int ir=0; ir<nxy; ir++)
                {
                    zpiece[ir] = v_elecstat[ir*this->rho_basis_->nplane+iz-this->rho_basis_->startz_current];
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

    ModuleBase::timer::tick("Potential","write_elecstat_pot");
    return;
}

}//namespace elecstate