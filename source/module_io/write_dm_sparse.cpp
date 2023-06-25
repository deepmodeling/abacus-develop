#include "module_io/write_dm_sparse.h"

#include "module_base/blas_connector.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void ModuleIO::write_dm1(const int &is, const int &istep, double** dm2d, const Parallel_Orbitals* ParaV,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &DMR_sparse)
{
    ModuleBase::timer::tick("ModuleIO","write_dm");
    ModuleBase::TITLE("ModuleIO","write_dm");

    get_dm_sparse(is, dm2d, ParaV, DMR_sparse);
    write_dm_sparse(is, istep, ParaV, DMR_sparse);
    destroy_dm_sparse(DMR_sparse);

    ModuleBase::timer::tick("ModuleIO","write_dm");
}

void ModuleIO::get_dm_sparse(const int &is, double** dm2d, const Parallel_Orbitals* ParaV,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &DMR_sparse)
{
    ModuleBase::timer::tick("ModuleIO","get_dm_sparse");
    ModuleBase::TITLE("ModuleIO","get_dm_sparse");

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

                            temp_value_double = dm2d[is][index];
                            if (std::abs(temp_value_double) > sparse_threshold)
                            {
                                DMR_sparse[dR][iw1_all][iw2_all] = temp_value_double;
                            }

                            ++index;
                        }
                    }
                }
            }
        }
    }

    ModuleBase::timer::tick("ModuleIO","get_dm_sparse");
}

void ModuleIO::write_dm_sparse(const int &is, const int &istep, const Parallel_Orbitals* ParaV,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &DMR_sparse)
{
    ModuleBase::timer::tick("ModuleIO","write_dm_sparse");
    ModuleBase::TITLE("ModuleIO","write_dm_sparse");

    double sparse_threshold = 1e-10;
    int step = istep;

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
    if(GlobalV::CALCULATION == "md" && !GlobalV::out_app_flag)
    {
        ssdm << GlobalV::global_matrix_dir << step << "_" << "data-DMR-sparse_SPIN" << is << ".csr";
    }
    else
    {
        ssdm << GlobalV::global_out_dir << "data-DMR-sparse_SPIN" << is << ".csr";
    }
    std::ofstream g1;

    if(GlobalV::DRANK==0)
    {
        if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
        {
            g1.open(ssdm.str().c_str(), ios::app);
        }
        else
        {
            g1.open(ssdm.str().c_str());
        }
        g1 << "STEP: " << step << std::endl;
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
            ModuleIO::output_single_R(g1, DMR_sparse[R_coor], sparse_threshold,false, *ParaV);
        }

        count++;

    }

    if(GlobalV::DRANK==0) 
    {
        g1.close();
    }
    
    delete[] DMR_nonzero_num;
    DMR_nonzero_num = nullptr;

    ModuleBase::timer::tick("ModuleIO","write_dm_sparse");
}

void ModuleIO::destroy_dm_sparse(std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &DMR_sparse)
{
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_DMR_sparse;
    DMR_sparse.swap(empty_DMR_sparse);
}
