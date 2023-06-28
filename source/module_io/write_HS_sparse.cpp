#include "write_HS_sparse.h"

#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "single_R_io.h"

void ModuleIO::save_HSR_sparse(
    const int &istep,
    LCAO_Matrix &lm,
    const double& sparse_threshold,
    const bool &binary,  
    const std::string &SR_filename, 
    const std::string &HR_filename_up, 
    const std::string &HR_filename_down = ""
)
{
    ModuleBase::TITLE("ModuleIO","save_HSR_sparse");
    ModuleBase::timer::tick("ModuleIO","save_HSR_sparse");

    auto &all_R_coor_ptr = lm.all_R_coor;
    auto &output_R_coor_ptr = lm.output_R_coor;
    auto &HR_sparse_ptr = lm.HR_sparse;
    auto &SR_sparse_ptr = lm.SR_sparse;
    auto &HR_soc_sparse_ptr = lm.HR_soc_sparse;
    auto &SR_soc_sparse_ptr = lm.SR_soc_sparse;

    int total_R_num = all_R_coor_ptr.size();
    int output_R_number = 0;
    int *H_nonzero_num[2] = {nullptr, nullptr};
    int *S_nonzero_num = nullptr;
    int step = istep;

    S_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(S_nonzero_num, total_R_num);

    int spin_loop = 1;
    if (GlobalV::NSPIN == 2)
    {
        spin_loop = 2;
    }

    for (int ispin = 0; ispin < spin_loop; ++ispin)
    {
        H_nonzero_num[ispin] = new int[total_R_num];
        ModuleBase::GlobalFunc::ZEROS(H_nonzero_num[ispin], total_R_num);
    }

    int count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        if (GlobalV::NSPIN != 4)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                auto iter = HR_sparse_ptr[ispin].find(R_coor);
                if (iter != HR_sparse_ptr[ispin].end())
                {
                    for (auto &row_loop : iter->second)
                    {
                        H_nonzero_num[ispin][count] += row_loop.second.size();
                    }
                }
            }

            auto iter = SR_sparse_ptr.find(R_coor);
            if (iter != SR_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    S_nonzero_num[count] += row_loop.second.size();
                }
            }
        }
        else
        {
            auto iter = HR_soc_sparse_ptr.find(R_coor);
            if (iter != HR_soc_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    H_nonzero_num[0][count] += row_loop.second.size();
                }
            }

            iter = SR_soc_sparse_ptr.find(R_coor);
            if (iter != SR_soc_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    S_nonzero_num[count] += row_loop.second.size();
                }
            }
        }

        count++;
    }

    Parallel_Reduce::reduce_int_all(S_nonzero_num, total_R_num);
    for (int ispin = 0; ispin < spin_loop; ++ispin)
    {
        Parallel_Reduce::reduce_int_all(H_nonzero_num[ispin], total_R_num);
    }

    if (GlobalV::NSPIN == 2)
    {
        for (int index = 0; index < total_R_num; ++index)
        {
            if (H_nonzero_num[0][index] != 0 || H_nonzero_num[1][index] != 0 || S_nonzero_num[index] != 0)
            {
                output_R_number++;
            }
        }
    }
    else
    {
        for (int index = 0; index < total_R_num; ++index)
        {
            if (H_nonzero_num[0][index] != 0 || S_nonzero_num[index] != 0)
            {
                output_R_number++;
            }
        }
    }

    std::stringstream ssh[2];
    std::stringstream sss;
    if(GlobalV::CALCULATION == "md" && !GlobalV::out_app_flag)
    {
        ssh[0] << GlobalV::global_matrix_dir << step << "_" << HR_filename_up;
        ssh[1] << GlobalV::global_matrix_dir << step << "_" << HR_filename_down;
        sss << GlobalV::global_matrix_dir << step << "_" << SR_filename;
    }
    else
    {
        ssh[0] << GlobalV::global_out_dir << HR_filename_up;
        ssh[1] << GlobalV::global_out_dir << HR_filename_down;
        sss << GlobalV::global_out_dir << SR_filename;
    }
    std::ofstream g1[2];
    std::ofstream g2;

    if(GlobalV::DRANK==0)
    {
        if (binary)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
                {
                    g1[ispin].open(ssh[ispin].str().c_str(), std::ios::binary | std::ios::app);
                }
                else
                {
                    g1[ispin].open(ssh[ispin].str().c_str(), std::ios::binary);
                }
                g1[ispin].write(reinterpret_cast<char *>(&step), sizeof(int));
                g1[ispin].write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
                g1[ispin].write(reinterpret_cast<char *>(&output_R_number), sizeof(int));
            }

            if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
            {
                g2.open(sss.str().c_str(), std::ios::binary | std::ios::app);
            }
            else
            {
                g2.open(sss.str().c_str(), std::ios::binary);
            }
            g2.write(reinterpret_cast<char *>(&step), sizeof(int));
            g2.write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
            g2.write(reinterpret_cast<char *>(&output_R_number), sizeof(int));
        }
        else
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
                {
                    g1[ispin].open(ssh[ispin].str().c_str(), std::ios::app);
                }
                else
                {
                    g1[ispin].open(ssh[ispin].str().c_str());
                }
                g1[ispin] << "STEP: " << step << std::endl;
                g1[ispin] << "Matrix Dimension of H(R): " << GlobalV::NLOCAL <<std::endl;
                g1[ispin] << "Matrix number of H(R): " << output_R_number << std::endl;
            }

            if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
            {
                g2.open(sss.str().c_str(), std::ios::app);
            }
            else
            {
                g2.open(sss.str().c_str());
            }
            g2 << "STEP: " << step <<std::endl;
            g2 << "Matrix Dimension of S(R): " << GlobalV::NLOCAL <<std::endl;
            g2 << "Matrix number of S(R): " << output_R_number << std::endl;
        }
    }

    output_R_coor_ptr.clear();

    count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (GlobalV::NSPIN == 2)
        {
            if (H_nonzero_num[0][count] == 0 && H_nonzero_num[1][count] == 0 && S_nonzero_num[count] == 0)
            {
                count++;
                continue;
            }
        }
        else
        {
            if (H_nonzero_num[0][count] == 0 && S_nonzero_num[count] == 0)
            {
                count++;
                continue;
            }
        }

        output_R_coor_ptr.insert(R_coor);

        if (GlobalV::DRANK == 0)
        {
            if (binary)
            {
                for (int ispin = 0; ispin < spin_loop; ++ispin)
                {
                    g1[ispin].write(reinterpret_cast<char *>(&dRx), sizeof(int));
                    g1[ispin].write(reinterpret_cast<char *>(&dRy), sizeof(int));
                    g1[ispin].write(reinterpret_cast<char *>(&dRz), sizeof(int));
                    g1[ispin].write(reinterpret_cast<char *>(&H_nonzero_num[ispin][count]), sizeof(int));
                }

                g2.write(reinterpret_cast<char *>(&dRx), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRy), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRz), sizeof(int));
                g2.write(reinterpret_cast<char *>(&S_nonzero_num[count]), sizeof(int));
            }
            else
            {
                for (int ispin = 0; ispin < spin_loop; ++ispin)
                {
                    g1[ispin] << dRx << " " << dRy << " " << dRz << " " << H_nonzero_num[ispin][count] << std::endl;
                }
                g2 << dRx << " " << dRy << " " << dRz << " " << S_nonzero_num[count] << std::endl;
            }
        }

        for (int ispin = 0; ispin < spin_loop; ++ispin)
        {
            if (H_nonzero_num[ispin][count] == 0)
            {
                // if (GlobalV::DRANK == 0)
                // {
                //     if (!binary)
                //     {
                //         g1[ispin] << std::endl;
                //         g1[ispin] << std::endl;
                //         for (int index = 0; index < GlobalV::NLOCAL+1; ++index)
                //         {
                //             g1[ispin] << 0 << " ";
                //         }
                //         g1[ispin] << std::endl;
                //     }
                // }
            }
            else
            {
                if (GlobalV::NSPIN != 4)
                {
                    output_single_R(g1[ispin], HR_sparse_ptr[ispin][R_coor], sparse_threshold, binary, *lm.ParaV);
                }
                else
                {
                    output_soc_single_R(g1[ispin], HR_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
                }
            }
        }

        if (S_nonzero_num[count] == 0)
        {
            // if (!binary)
            // {
            //     if (GlobalV::DRANK == 0)
            //     {
            //         g2 << std::endl;
            //         g2 << std::endl;
            //         for (int index = 0; index < GlobalV::NLOCAL+1; ++index)
            //         {
            //             g2 << 0 << " ";
            //         }
            //         g2 << std::endl;
            //     }
            // }
        }
        else
        {
            if (GlobalV::NSPIN != 4)
            {
                output_single_R(g2, SR_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
            }
            else
            {
                output_soc_single_R(g2, SR_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
            }
        }

        count++;

    }

    if(GlobalV::DRANK==0) 
    {
        for (int ispin = 0; ispin < spin_loop; ++ispin) g1[ispin].close();
        g2.close();
    }
    
    for (int ispin = 0; ispin < spin_loop; ++ispin) 
    {
        delete[] H_nonzero_num[ispin];
        H_nonzero_num[ispin] = nullptr;
    }
    delete[] S_nonzero_num;
    S_nonzero_num = nullptr;

    ModuleBase::timer::tick("ModuleIO","save_HSR_sparse");
    return;
}

void ModuleIO::save_SR_sparse(
    LCAO_Matrix &lm,
    const double& sparse_threshold,
    const bool &binary,  
    const std::string &SR_filename
)
{
    ModuleBase::TITLE("ModuleIO","save_SR_sparse");
    ModuleBase::timer::tick("ModuleIO","save_SR_sparse");

    auto &all_R_coor_ptr = lm.all_R_coor;
    auto &SR_sparse_ptr = lm.SR_sparse;
    auto &SR_soc_sparse_ptr = lm.SR_soc_sparse;

    int total_R_num = all_R_coor_ptr.size();
    int output_R_number = 0;
    int *S_nonzero_num = nullptr;

    S_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(S_nonzero_num, total_R_num);

    int count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        if (GlobalV::NSPIN != 4)
        {
            auto iter = SR_sparse_ptr.find(R_coor);
            if (iter != SR_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    S_nonzero_num[count] += row_loop.second.size();
                }
            }
        }
        else
        {
            auto iter = SR_soc_sparse_ptr.find(R_coor);
            if (iter != SR_soc_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    S_nonzero_num[count] += row_loop.second.size();
                }
            }
        }

        count++;
    }

    Parallel_Reduce::reduce_int_all(S_nonzero_num, total_R_num);

    for (int index = 0; index < total_R_num; ++index)
    {
        if (S_nonzero_num[index] != 0)
        {
            output_R_number++;
        }
    }

    std::stringstream sss;
    sss << SR_filename;
    std::ofstream g2;

    if(GlobalV::DRANK==0)
    {
        if (binary)
        {
            g2.open(sss.str().c_str(), std::ios::binary);
            g2.write(reinterpret_cast<char *>(0), sizeof(int));
            g2.write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
            g2.write(reinterpret_cast<char *>(&output_R_number), sizeof(int));
        }
        else
        {
            g2.open(sss.str().c_str());
            g2 << "STEP: " << 0 << std::endl;
            g2 << "Matrix Dimension of S(R): " << GlobalV::NLOCAL <<std::endl;
            g2 << "Matrix number of S(R): " << output_R_number << std::endl;
        }
    }

    count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (S_nonzero_num[count] == 0)
        {
            count++;
            continue;
        }

        if (GlobalV::DRANK == 0)
        {
            if (binary)
            {
                g2.write(reinterpret_cast<char *>(&dRx), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRy), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRz), sizeof(int));
                g2.write(reinterpret_cast<char *>(&S_nonzero_num[count]), sizeof(int));
            }
            else
            {
                g2 << dRx << " " << dRy << " " << dRz << " " << S_nonzero_num[count] << std::endl;
            }
        }

        if (GlobalV::NSPIN != 4)
        {
            output_single_R(g2, SR_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
        }
        else
        {
            output_soc_single_R(g2, SR_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
        }

        count++;

    }

    if(GlobalV::DRANK==0) 
    {
        g2.close();
    }

    delete[] S_nonzero_num;
    S_nonzero_num = nullptr;

    ModuleBase::timer::tick("ModuleIO","save_SR_sparse");
    return;
}

void ModuleIO::save_TR_sparse(
    const int &istep,
    LCAO_Matrix &lm,
    const double& sparse_threshold,
    const bool &binary,  
    const std::string &TR_filename
)
{
    ModuleBase::TITLE("ModuleIO","save_TR_sparse");
    ModuleBase::timer::tick("ModuleIO","save_TR_sparse");

    auto &all_R_coor_ptr = lm.all_R_coor;
    auto &TR_sparse_ptr = lm.TR_sparse;
    auto &TR_soc_sparse_ptr = lm.TR_soc_sparse;

    int total_R_num = all_R_coor_ptr.size();
    int output_R_number = 0;
    int *T_nonzero_num = nullptr;
    int step = istep;

    T_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(T_nonzero_num, total_R_num);

    int count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        if (GlobalV::NSPIN != 4)
        {
            auto iter = TR_sparse_ptr.find(R_coor);
            if (iter != TR_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    T_nonzero_num[count] += row_loop.second.size();
                }
            }
        }
        else
        {
            auto iter = TR_soc_sparse_ptr.find(R_coor);
            if (iter != TR_soc_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    T_nonzero_num[count] += row_loop.second.size();
                }
            }
        }

        count++;
    }

    Parallel_Reduce::reduce_int_all(T_nonzero_num, total_R_num);

    for (int index = 0; index < total_R_num; ++index)
    {
        if (T_nonzero_num[index] != 0)
        {
            output_R_number++;
        }
    }

    std::stringstream sss;
    sss << TR_filename;
    std::ofstream g2;

    if(GlobalV::DRANK==0)
    {
        if (binary)
        {
            if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
            {
                g2.open(sss.str().c_str(), std::ios::binary | std::ios::app);
            }
            else
            {
                g2.open(sss.str().c_str(), std::ios::binary);
            }
            g2.write(reinterpret_cast<char *>(&step), sizeof(int));
            g2.write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
            g2.write(reinterpret_cast<char *>(&output_R_number), sizeof(int));
        }
        else
        {
            if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
            {
                g2.open(sss.str().c_str(), std::ios::app);
            }
            else
            {
                g2.open(sss.str().c_str());
            }
            g2 << "STEP: " << step << std::endl;
            g2 << "Matrix Dimension of T(R): " << GlobalV::NLOCAL <<std::endl;
            g2 << "Matrix number of T(R): " << output_R_number << std::endl;
        }
    }

    count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (T_nonzero_num[count] == 0)
        {
            count++;
            continue;
        }

        if (GlobalV::DRANK == 0)
        {
            if (binary)
            {
                g2.write(reinterpret_cast<char *>(&dRx), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRy), sizeof(int));
                g2.write(reinterpret_cast<char *>(&dRz), sizeof(int));
                g2.write(reinterpret_cast<char *>(&T_nonzero_num[count]), sizeof(int));
            }
            else
            {
                g2 << dRx << " " << dRy << " " << dRz << " " << T_nonzero_num[count] << std::endl;
            }
        }

        if (GlobalV::NSPIN != 4)
        {
            output_single_R(g2, TR_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
        }
        else
        {
            output_soc_single_R(g2, TR_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
        }

        count++;

    }

    if(GlobalV::DRANK==0) 
    {
        g2.close();
    }

    delete[] T_nonzero_num;
    T_nonzero_num = nullptr;

    ModuleBase::timer::tick("ModuleIO","save_TR_sparse");
    return;
}

void ModuleIO::save_dH_sparse(
    const int &istep,
    LCAO_Matrix &lm,
    const double& sparse_threshold,
    const bool &binary
)
{
    ModuleBase::TITLE("ModuleIO","save_dH_sparse");
    ModuleBase::timer::tick("ModuleIO","save_dH_sparse");

    auto &all_R_coor_ptr = lm.all_R_coor;
    auto &output_R_coor_ptr = lm.output_R_coor;
    auto &dHRx_sparse_ptr = lm.dHRx_sparse;
    auto &dHRx_soc_sparse_ptr = lm.dHRx_soc_sparse;
    auto &dHRy_sparse_ptr = lm.dHRy_sparse;
    auto &dHRy_soc_sparse_ptr = lm.dHRy_soc_sparse;
    auto &dHRz_sparse_ptr = lm.dHRz_sparse;
    auto &dHRz_soc_sparse_ptr = lm.dHRz_soc_sparse;

    int total_R_num = all_R_coor_ptr.size();
    int output_R_number = 0;
    int *dHx_nonzero_num[2] = {nullptr, nullptr};
    int *dHy_nonzero_num[2] = {nullptr, nullptr};
    int *dHz_nonzero_num[2] = {nullptr, nullptr};
    int step = istep;

    int spin_loop = 1;
    if (GlobalV::NSPIN == 2)
    {
        spin_loop = 2;
    }

    for (int ispin = 0; ispin < spin_loop; ++ispin)
    {
        dHx_nonzero_num[ispin] = new int[total_R_num];
        ModuleBase::GlobalFunc::ZEROS(dHx_nonzero_num[ispin], total_R_num);
        dHy_nonzero_num[ispin] = new int[total_R_num];
        ModuleBase::GlobalFunc::ZEROS(dHy_nonzero_num[ispin], total_R_num);
        dHz_nonzero_num[ispin] = new int[total_R_num];
        ModuleBase::GlobalFunc::ZEROS(dHz_nonzero_num[ispin], total_R_num);                
    }

    int count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        if (GlobalV::NSPIN != 4)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                auto iter1 = dHRx_sparse_ptr[ispin].find(R_coor);
                if (iter1 != dHRx_sparse_ptr[ispin].end())
                {
                    for (auto &row_loop : iter1->second)
                    {
                        dHx_nonzero_num[ispin][count] += row_loop.second.size();
                    }
                }
                
                auto iter2 = dHRy_sparse_ptr[ispin].find(R_coor);
                if (iter2 != dHRy_sparse_ptr[ispin].end())
                {
                    for (auto &row_loop : iter2->second)
                    {
                        dHy_nonzero_num[ispin][count] += row_loop.second.size();
                    }
                }
                
                auto iter3 = dHRz_sparse_ptr[ispin].find(R_coor);
                if (iter3 != dHRz_sparse_ptr[ispin].end())
                {
                    for (auto &row_loop : iter3->second)
                    {
                        dHz_nonzero_num[ispin][count] += row_loop.second.size();
                    }
                }
            }
        }
        else
        {
            auto iter = dHRx_soc_sparse_ptr.find(R_coor);
            if (iter != dHRx_soc_sparse_ptr.end())
            {
                for (auto &row_loop : iter->second)
                {
                    dHx_nonzero_num[0][count] += row_loop.second.size();
                }
            }
        }

        count++;
    }

    for (int ispin = 0; ispin < spin_loop; ++ispin)
    {
        Parallel_Reduce::reduce_int_all(dHx_nonzero_num[ispin], total_R_num);
        Parallel_Reduce::reduce_int_all(dHy_nonzero_num[ispin], total_R_num);
        Parallel_Reduce::reduce_int_all(dHz_nonzero_num[ispin], total_R_num); 
    }

    if (GlobalV::NSPIN == 2)
    {
        for (int index = 0; index < total_R_num; ++index)
        {
            if (dHx_nonzero_num[0][index] != 0 || dHx_nonzero_num[1][index] != 0 ||
                dHy_nonzero_num[0][index] != 0 || dHy_nonzero_num[1][index] != 0 ||
                dHz_nonzero_num[0][index] != 0 || dHz_nonzero_num[1][index] != 0)
            {
                output_R_number++;
            }
        }
    }
    else
    {
        for (int index = 0; index < total_R_num; ++index)
        {
            if (dHx_nonzero_num[0][index] != 0 || dHy_nonzero_num[0][index] != 0 || dHz_nonzero_num[0][index] != 0)
            {
                output_R_number++;
            }
        }
    }

    std::stringstream sshx[2];
    std::stringstream sshy[2];
    std::stringstream sshz[2];
    if(GlobalV::CALCULATION == "md" && !GlobalV::out_app_flag)
    {
        sshx[0] << GlobalV::global_matrix_dir << step << "_" << "data-dHRx-sparse_SPIN0.csr";
        sshx[1] << GlobalV::global_matrix_dir << step << "_" << "data-dHRx-sparse_SPIN1.csr";
        sshy[0] << GlobalV::global_matrix_dir << step << "_" << "data-dHRy-sparse_SPIN0.csr";
        sshy[1] << GlobalV::global_matrix_dir << step << "_" << "data-dHRy-sparse_SPIN1.csr";
        sshz[0] << GlobalV::global_matrix_dir << step << "_" << "data-dHRz-sparse_SPIN0.csr";
        sshz[1] << GlobalV::global_matrix_dir << step << "_" << "data-dHRz-sparse_SPIN1.csr";                
    }
    else
    {
        sshx[0] << GlobalV::global_out_dir << "data-dHRx-sparse_SPIN0.csr";
        sshx[1] << GlobalV::global_out_dir << "data-dHRx-sparse_SPIN1.csr";
        sshy[0] << GlobalV::global_out_dir << "data-dHRy-sparse_SPIN0.csr";
        sshy[1] << GlobalV::global_out_dir << "data-dHRy-sparse_SPIN1.csr";
        sshz[0] << GlobalV::global_out_dir << "data-dHRz-sparse_SPIN0.csr";
        sshz[1] << GlobalV::global_out_dir << "data-dHRz-sparse_SPIN1.csr";                
    }
    std::ofstream g1x[2];
    std::ofstream g1y[2];
    std::ofstream g1z[2];

    if(GlobalV::DRANK==0)
    {
        if (binary)
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
                {
                    g1x[ispin].open(sshx[ispin].str().c_str(), std::ios::binary | std::ios::app);
                    g1y[ispin].open(sshy[ispin].str().c_str(), std::ios::binary | std::ios::app);
                    g1z[ispin].open(sshz[ispin].str().c_str(), std::ios::binary | std::ios::app);
                }
                else
                {
                    g1x[ispin].open(sshx[ispin].str().c_str(), std::ios::binary);
                    g1y[ispin].open(sshy[ispin].str().c_str(), std::ios::binary);
                    g1z[ispin].open(sshz[ispin].str().c_str(), std::ios::binary);
                }

                g1x[ispin].write(reinterpret_cast<char *>(&step), sizeof(int));
                g1x[ispin].write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
                g1x[ispin].write(reinterpret_cast<char *>(&output_R_number), sizeof(int));

                g1y[ispin].write(reinterpret_cast<char *>(&step), sizeof(int));
                g1y[ispin].write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
                g1y[ispin].write(reinterpret_cast<char *>(&output_R_number), sizeof(int));

                g1z[ispin].write(reinterpret_cast<char *>(&step), sizeof(int));
                g1z[ispin].write(reinterpret_cast<char *>(&GlobalV::NLOCAL), sizeof(int));
                g1z[ispin].write(reinterpret_cast<char *>(&output_R_number), sizeof(int));                                
            }
        }
        else
        {
            for (int ispin = 0; ispin < spin_loop; ++ispin)
            {
                if(GlobalV::CALCULATION == "md" && GlobalV::out_app_flag && step)
                {
                    g1x[ispin].open(sshx[ispin].str().c_str(), std::ios::app);
                    g1y[ispin].open(sshy[ispin].str().c_str(), std::ios::app);
                    g1z[ispin].open(sshz[ispin].str().c_str(), std::ios::app);
                }
                else
                {
                    g1x[ispin].open(sshx[ispin].str().c_str());
                    g1y[ispin].open(sshy[ispin].str().c_str());
                    g1z[ispin].open(sshz[ispin].str().c_str());
                }

                g1x[ispin] << "STEP: " << step << std::endl;
                g1x[ispin] << "Matrix Dimension of dHx(R): " << GlobalV::NLOCAL <<std::endl;
                g1x[ispin] << "Matrix number of dHx(R): " << output_R_number << std::endl;

                g1y[ispin] << "STEP: " << step << std::endl;
                g1y[ispin] << "Matrix Dimension of dHy(R): " << GlobalV::NLOCAL <<std::endl;
                g1y[ispin] << "Matrix number of dHy(R): " << output_R_number << std::endl;

                g1z[ispin] << "STEP: " << step << std::endl;
                g1z[ispin] << "Matrix Dimension of dHz(R): " << GlobalV::NLOCAL <<std::endl;
                g1z[ispin] << "Matrix number of dHz(R): " << output_R_number << std::endl;                                
            }
        }
    }

    output_R_coor_ptr.clear();

    count = 0;
    for (auto &R_coor : all_R_coor_ptr)
    {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (GlobalV::NSPIN == 2)
        {
            if (dHx_nonzero_num[0][count] == 0 && dHx_nonzero_num[1][count] == 0 &&
                dHy_nonzero_num[0][count] == 0 && dHy_nonzero_num[1][count] == 0 &&
                dHz_nonzero_num[0][count] == 0 && dHz_nonzero_num[1][count] == 0)
            {
                count++;
                continue;
            }
        }
        else
        {
            if (dHx_nonzero_num[0][count] == 0 && dHy_nonzero_num[0][count] == 0 && dHz_nonzero_num[0][count] == 0)
            {
                count++;
                continue;
            }
        }

        output_R_coor_ptr.insert(R_coor);

        if (GlobalV::DRANK == 0)
        {
            if (binary)
            {
                for (int ispin = 0; ispin < spin_loop; ++ispin)
                {
                    g1x[ispin].write(reinterpret_cast<char *>(&dRx), sizeof(int));
                    g1x[ispin].write(reinterpret_cast<char *>(&dRy), sizeof(int));
                    g1x[ispin].write(reinterpret_cast<char *>(&dRz), sizeof(int));
                    g1x[ispin].write(reinterpret_cast<char *>(&dHx_nonzero_num[ispin][count]), sizeof(int));

                    g1y[ispin].write(reinterpret_cast<char *>(&dRx), sizeof(int));
                    g1y[ispin].write(reinterpret_cast<char *>(&dRy), sizeof(int));
                    g1y[ispin].write(reinterpret_cast<char *>(&dRz), sizeof(int));
                    g1y[ispin].write(reinterpret_cast<char *>(&dHy_nonzero_num[ispin][count]), sizeof(int));

                    g1z[ispin].write(reinterpret_cast<char *>(&dRx), sizeof(int));
                    g1z[ispin].write(reinterpret_cast<char *>(&dRy), sizeof(int));
                    g1z[ispin].write(reinterpret_cast<char *>(&dRz), sizeof(int));
                    g1z[ispin].write(reinterpret_cast<char *>(&dHz_nonzero_num[ispin][count]), sizeof(int));                                        
                }
            }
            else
            {
                for (int ispin = 0; ispin < spin_loop; ++ispin)
                {
                    g1x[ispin] << dRx << " " << dRy << " " << dRz << " " << dHx_nonzero_num[ispin][count] << std::endl;
                    g1y[ispin] << dRx << " " << dRy << " " << dRz << " " << dHy_nonzero_num[ispin][count] << std::endl;
                    g1z[ispin] << dRx << " " << dRy << " " << dRz << " " << dHz_nonzero_num[ispin][count] << std::endl;
                }
            }
        }

        for (int ispin = 0; ispin < spin_loop; ++ispin)
        {
            if (dHx_nonzero_num[ispin][count] > 0)
            {
                if (GlobalV::NSPIN != 4)
                {
                    output_single_R(g1x[ispin], dHRx_sparse_ptr[ispin][R_coor], sparse_threshold, binary, *lm.ParaV);
                }
                else
                {
                    output_soc_single_R(g1x[ispin], dHRx_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
                }
            }
            if (dHy_nonzero_num[ispin][count] > 0)
            {
                if (GlobalV::NSPIN != 4)
                {
                    output_single_R(g1y[ispin], dHRy_sparse_ptr[ispin][R_coor], sparse_threshold, binary, *lm.ParaV);
                }
                else
                {
                    output_soc_single_R(g1y[ispin], dHRy_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
                }
            }
            if (dHz_nonzero_num[ispin][count] > 0)
            {
                if (GlobalV::NSPIN != 4)
                {
                    output_single_R(g1z[ispin], dHRz_sparse_ptr[ispin][R_coor], sparse_threshold, binary, *lm.ParaV);
                }
                else
                {
                    output_soc_single_R(g1z[ispin], dHRz_soc_sparse_ptr[R_coor], sparse_threshold, binary, *lm.ParaV);
                }
            }                        
        }

          count++;

    }

    if(GlobalV::DRANK==0) 
    {
        for (int ispin = 0; ispin < spin_loop; ++ispin) g1x[ispin].close();
        for (int ispin = 0; ispin < spin_loop; ++ispin) g1y[ispin].close();
        for (int ispin = 0; ispin < spin_loop; ++ispin) g1z[ispin].close();
    }
    
    for (int ispin = 0; ispin < spin_loop; ++ispin) 
    {
        delete[] dHx_nonzero_num[ispin];
        dHx_nonzero_num[ispin] = nullptr;
        delete[] dHy_nonzero_num[ispin];
        dHy_nonzero_num[ispin] = nullptr;
        delete[] dHz_nonzero_num[ispin];
        dHz_nonzero_num[ispin] = nullptr;                
    }

    ModuleBase::timer::tick("ModuleIO","save_dH_sparse");
    return;
}
