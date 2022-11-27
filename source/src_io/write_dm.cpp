#include "../src_lcao/local_orbital_charge.h"
#include "../src_pw/global.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/blas_connector.h"
#include "../module_base/timer.h"

//-------------------------------------------------
// NOTE for Local_Orbital_Charge::write_dm
// I will give an example here, suppose we have a 4*4 
// density matrix (symmetry) which is
// 1.1 2.3 3.6 4.2
// 2.3 5.2 7.1 8.9
// 3.6 7.1 1.3 3.2  
// 4.2 8.9 3.2 2.4
// we use two processors, each one has 3 orbitals
// processor 1 has orbital index 1, 2, 4:
// ('no' means no value on this processor)
// 1.1 2.3 no  4.2
// 2.3 5.2 no  8.9
// no  no  no  no   
// 4.2 8.9 no  2.4
// processor 2 has orbital index 1, 3, 4;
// 1.1 no  3.6 4.2
// no  no  no  no 
// 3.6 no  1.3 3.2  
// 4.2 no  3.2 2.4
// now we want to reduce them and print out,
// we plan to reduce them one row by one row,
// then for the first row, we need to set the
// temparary array to 4 (GlobalV::NLOCAL in code),
// then we reduce first row, it will become
// 2.2 2.3 3.6 8.4,
// we notice that the first element and fourth
// element is doubled, that's because the density
// may have overlap, so we need to first count
// for each element, how many times it is duplicated
// on other processors, this is why there is
// a 'count' integer array in the code.
// UPDATED BY MOHAN 2014-05-18
void Local_Orbital_Charge::write_dm(
	const int &is, 
	const int &iter, 
	const std::string &fn, 
	const int &precision)
{
    ModuleBase::TITLE("Local_Orbital_Charge","write_dm");

	if (out_dm==0)
	{
		return;
	}
	else if(iter % out_dm != 0)
	{
		return; 
	}
	ModuleBase::timer::tick("Local_Orbital_Charge","write_dm");

	time_t start, end;
	std::ofstream ofs;

	if(GlobalV::MY_RANK==0)
	{
		start = time(NULL);

		ofs.open(fn.c_str());
		if (!ofs)
		{
			ModuleBase::WARNING("Charge::write_rho","Can't create Charge File!");
		}

		//GlobalV::ofs_running << "\n Output charge file." << std::endl;

		ofs << GlobalC::ucell.latName << std::endl;//1
		ofs << " " << GlobalC::ucell.lat0 * ModuleBase::BOHR_TO_A << std::endl;
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
			Atom* atom = &GlobalC::ucell.atoms[it];
			ofs << std::setprecision(15);
			for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				ofs << " " << atom->taud[ia].x
					<< " " << atom->taud[ia].y
					<< " " << atom->taud[ia].z << std::endl;
			}
		}

		ofs << "\n " << GlobalV::NSPIN;
		if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
		{
			ofs << "\n " << GlobalC::en.ef << " (fermi energy)";
		}
		else if(GlobalV::NSPIN==2)
		{
			if(is==0)ofs << "\n " << GlobalC::en.ef_up << " (fermi energy for spin=1)";
			else if(is==1)ofs << "\n " << GlobalC::en.ef_dw << " (fermi energy for spin=2)";
		}
		else
		{
			ModuleBase::WARNING_QUIT("write_rho","check nspin!");
		}


		ofs << "\n  " << GlobalV::NLOCAL << " " << GlobalV::NLOCAL << std::endl;

		ofs << std::setprecision(precision);
		ofs << scientific;

	}

    //ofs << "\n " << GlobalV::GAMMA_ONLY_LOCAL << " (GAMMA ONLY LOCAL)" << std::endl;
#ifndef __MPI
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        for(int i=0; i<GlobalV::NLOCAL; ++i)
        {
            for(int j=0; j<GlobalV::NLOCAL; ++j)
            {
                if(j%8==0) ofs << "\n";
                ofs << " " << this->DM[is][i][j];
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("write_dm","not ready yet");
        ofs << " " << GlobalC::GridT.nnrg << " (nnrg)" << std::endl;
        for(int i=0; i<GlobalC::GridT.nnrg; ++i)
        {
            if(i%8==0) ofs << "\n";
            ofs << " " << this->DM_R[is][i];
        }
    }
#else
    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        //xiaohui modify 2014-06-18
        
        double* tmp = new double[GlobalV::NLOCAL];
        int* count = new int[GlobalV::NLOCAL];
        for (int i=0; i<GlobalV::NLOCAL; ++i)
        {
            // when reduce, there may be 'redundance', we need to count them.
            ModuleBase::GlobalFunc::ZEROS(count, GlobalV::NLOCAL);
            const int mu = GlobalC::GridT.trace_lo[i];
            if (mu >= 0)
            {
                for (int j=0; j<GlobalV::NLOCAL; ++j)
                {
                    const int nu = GlobalC::GridT.trace_lo[j];
                    if (nu >= 0)
                    {
                        count[j]=1; 
                    }
                }
            }
            Parallel_Reduce::reduce_int_all( count, GlobalV::NLOCAL );

            // reduce the density matrix for 'i' line.
            ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);
            if (mu >= 0)
            {
                for (int j=0; j<GlobalV::NLOCAL; j++)
                {
                    const int nu = GlobalC::GridT.trace_lo[j];
                    if (nu >=0)
                    {
                        tmp[j] = DM[is][mu][nu];
                        //GlobalV::ofs_running << " dmi=" << i << " j=" << j << " " << DM[is][mu][nu] << std::endl;
                    }
                }
            }
            Parallel_Reduce::reduce_double_all( tmp, GlobalV::NLOCAL );

            if(GlobalV::MY_RANK==0)
            {
                for (int j=0; j<GlobalV::NLOCAL; j++)
                {
                    if(j%8==0) ofs << "\n";
                    if(count[j]>0)
                    {
                        ofs << " " << tmp[j]/(double)count[j];
                    }
                    else
                    {
                        ofs << " 0"; 
                    }
                }
            }
        }
        delete[] tmp;
        delete[] count;
        
        //xiaohui add 2014-06-18
        //for(int i=0; i<GlobalV::NLOCAL; ++i)
        //{
        //  for(int j=0; j<GlobalV::NLOCAL; ++j)
        //  {
        //      if(j%8==0) ofs << "\n";
        //      ofs << " " << this->DM[is][i][j];
        //  }
        //}

    }
    else
    {
        ofs << " " << GlobalC::GridT.nnrg << " (nnrg)" << std::endl;
        ModuleBase::WARNING_QUIT("local_orbital_charge","not ready to output DM_R");
    }
#endif
	if(GlobalV::MY_RANK==0)
	{
		end = time(NULL);
		ModuleBase::GlobalFunc::OUT_TIME("write_rho",start,end);
		ofs.close();
	}
	ModuleBase::timer::tick("Local_Orbital_Charge","write_dm");

    return;
}

void Local_Orbital_Charge::write_dm1(const int &is, const int &istep)
{
    ModuleBase::timer::tick("Local_Orbital_Charge","write_dm");
    ModuleBase::TITLE("Local_Orbital_Charge","write_dm");

    get_dm_sparse(is);
    write_dm_sparse(is, istep);
    destroy_dm_sparse();

    ModuleBase::timer::tick("Local_Orbital_Charge","write_dm");
}

inline void output_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, double>> &XR, const double &sparse_threshold, const Parallel_Orbitals &pv)
{
    double *line = nullptr;
    std::vector<int> indptr;
    indptr.reserve(GlobalV::NLOCAL + 1);
    indptr.push_back(0);

    std::stringstream tem1;
    tem1 << GlobalV::global_out_dir << "temp_sparse_indices.dat";
    std::ofstream ofs_tem1;
    std::ifstream ifs_tem1;

    if (GlobalV::DRANK == 0)
    {
        ofs_tem1.open(tem1.str().c_str());
    }

    line = new double[GlobalV::NLOCAL];
    for(int row = 0; row < GlobalV::NLOCAL; ++row)
    {
        // line = new double[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(line, GlobalV::NLOCAL);

        if(pv.trace_loc_row[row] >= 0)
        {
            auto iter = XR.find(row);
            if (iter != XR.end())
            {
                for (auto &value : iter->second)
                {
                    line[value.first] = value.second;
                }
            }
        }

        Parallel_Reduce::reduce_double_all(line, GlobalV::NLOCAL);

        if(GlobalV::DRANK == 0)
        {
            int nonzeros_count = 0;
            for (int col = 0; col < GlobalV::NLOCAL; ++col)
            {
                if (std::abs(line[col]) > sparse_threshold)
                {
                    ofs << " " << fixed << scientific << std::setprecision(8) << line[col];
                    ofs_tem1 << " " << col;

                    nonzeros_count++;
                }

            }
            nonzeros_count += indptr.back();
            indptr.push_back(nonzeros_count);
        }

        // delete[] line;
        // line = nullptr;

    }

    delete[] line;
    line = nullptr;

    if (GlobalV::DRANK == 0)
    {
        ofs << std::endl;
        ofs_tem1 << std::endl;
        ofs_tem1.close();
        ifs_tem1.open(tem1.str().c_str());
        ofs << ifs_tem1.rdbuf();
        ifs_tem1.close();
        for (auto &i : indptr)
        {
            ofs << " " << i;
        }
        ofs << std::endl;
        
        std::remove(tem1.str().c_str());

    }

}

void Local_Orbital_Charge::get_dm_sparse(const int &is)
{
    ModuleBase::timer::tick("Local_Orbital_Charge","get_dm_sparse");
    ModuleBase::TITLE("Local_Orbital_Charge","get_dm_sparse");

    double sparse_threshold = 1e-10;

    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    double temp_value_double;

    for(int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) adj = true;
                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if( distance1 < rcut1 && distance2 < rcut2 )
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if(adj)
                {
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);

                    Abfs::Vector3_Order<int> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);

                    for(int ii=0; ii<atom1->nw; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = ParaV->trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = ParaV->trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            temp_value_double = this->DM_R[is][index];
                            if (std::abs(temp_value_double) > sparse_threshold)
                            {
                                this->DMR_sparse[dR][iw1_all][iw2_all] = temp_value_double;
                            }

                            ++index;
                        }
                    }
                }
            }
        }
    }

    ModuleBase::timer::tick("Local_Orbital_Charge","get_dm_sparse");
}

void Local_Orbital_Charge::write_dm_sparse(const int &is, const int &istep)
{
    ModuleBase::timer::tick("Local_Orbital_Charge","write_dm_sparse");
    ModuleBase::TITLE("Local_Orbital_Charge","write_dm_sparse");

    double sparse_threshold = 1e-10;

//get list of R
    int R_minX = int(GlobalC::GridD.getD_minX());
    int R_minY = int(GlobalC::GridD.getD_minY());
    int R_minZ = int(GlobalC::GridD.getD_minZ());

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    for(int ix = 0; ix < R_x; ix++)
    {
        for(int iy = 0; iy < R_y; iy++)
        {
            for(int iz = 0; iz < R_z; iz++)
            {
                Abfs::Vector3_Order<int> temp_R(ix+R_minX, iy+R_minY, iz+R_minZ);
                all_R_coor.insert(temp_R);
            }
        }
    }

    int total_R_num = all_R_coor.size();
    int output_R_number = 0;
    int *DMR_nonzero_num = {nullptr};

    DMR_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(DMR_nonzero_num, total_R_num);

    int count = 0;
    for (auto &R_coor : all_R_coor)
    {
        auto iter = DMR_sparse.find(R_coor);
        if (iter != DMR_sparse.end())
        {
            for (auto &row_loop : iter->second)
            {
                DMR_nonzero_num[count] += row_loop.second.size();
            }
        }

        count++;
    }

    Parallel_Reduce::reduce_int_all(DMR_nonzero_num, total_R_num);

    for (int index = 0; index < total_R_num; ++index)
    {
        if (DMR_nonzero_num[index] != 0)
        {
            output_R_number++;
        }
    }

    std::stringstream ssdm;
    if(GlobalV::CALCULATION == "md")
    {
        ssdm << GlobalV::global_matrix_dir << istep << "_" << "data-DMR-sparse_SPIN" << is << ".csr";
    }
    else
    {
        ssdm << GlobalV::global_out_dir << "data-DMR-sparse_SPIN" << is << ".csr";
    }
    std::ofstream g1;

    if(GlobalV::DRANK==0)
    {
        g1.open(ssdm.str().c_str(), ios::app);
        g1 << "STEP: " << istep << std::endl;
        g1 << "Matrix Dimension of DM(R): " << GlobalV::NLOCAL <<std::endl;
        g1 << "Matrix number of DM(R): " << output_R_number << std::endl;
    }

    count = 0;
    for (auto &R_coor : all_R_coor)
    {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (DMR_nonzero_num[count] == 0)
        {
            count++;
            continue;
        }

        if (GlobalV::DRANK == 0)
        {
            g1 << dRx << " " << dRy << " " << dRz << " " << DMR_nonzero_num[count] << std::endl;
        }

        if (DMR_nonzero_num[count] != 0)
        {
            output_single_R(g1, DMR_sparse[R_coor], sparse_threshold, *ParaV);
        }

        count++;

    }

    if(GlobalV::DRANK==0) 
    {
        g1.close();
    }
    
    delete[] DMR_nonzero_num;
    DMR_nonzero_num = nullptr;

    ModuleBase::timer::tick("Local_Orbital_Charge","write_dm_sparse");
}

void Local_Orbital_Charge::destroy_dm_sparse()
{
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_DMR_sparse;
    DMR_sparse.swap(empty_DMR_sparse);
}