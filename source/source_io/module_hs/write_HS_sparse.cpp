#include "write_HS_sparse.h"

#include "source_base/global_function.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_rt/td_info.h"
#include "single_R_io.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace
{
template <typename Tdata>
std::vector<long long> count_nonzeros_by_R(
    const ModuleIO::SparseRMatrix<Tdata>& smat,
    const std::set<ModuleIO::RCoordinate>& all_R_coor,
    const double threshold,
    const bool reduce)
{
    std::vector<long long> nonzero_num(all_R_coor.size(), 0);
    int count = 0;
    for (const auto& R_coor: all_R_coor)
    {
        const auto iter = smat.find(R_coor);
        if (iter != smat.end())
        {
            for (const auto& row_loop: iter->second)
            {
                for (const auto& col_value: row_loop.second)
                {
                    if (std::abs(col_value.second) > threshold)
                    {
                        ++nonzero_num[count];
                    }
                }
            }
        }
        ++count;
    }

    if (reduce)
    {
        Parallel_Reduce::reduce_all(nonzero_num.data(), static_cast<int>(nonzero_num.size()));
    }
    return nonzero_num;
}

int count_output_R(const std::vector<long long>& nonzero_num)
{
    int output_R_number = 0;
    for (const long long count: nonzero_num)
    {
        if (count != 0)
        {
            ++output_R_number;
        }
    }
    return output_R_number;
}

void open_sparse_file(std::ofstream& ofs, const ModuleIO::SparseWriteOptions& options)
{
    std::ios_base::openmode mode = std::ios::out;
    if (options.binary)
    {
        mode |= std::ios::binary;
    }
    if (PARAM.inp.calculation == "md" && PARAM.inp.out_app_flag && options.istep)
    {
        mode |= std::ios::app;
    }
    ofs.open(options.filename.c_str(), mode);
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::open_sparse_file",
                                 "Cannot open sparse matrix file: " + options.filename);
    }
}

void write_sparse_header(std::ofstream& ofs,
                         const ModuleIO::SparseWriteOptions& options,
                         const int nlocal,
                         const int output_R_number)
{
    const int step = std::max(options.istep, 0);
    if (options.binary)
    {
        ofs.write(reinterpret_cast<const char*>(&step), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&nlocal), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&output_R_number), sizeof(int));
    }
    else
    {
        ofs << "STEP: " << step << std::endl;
        ofs << "Matrix Dimension of " + options.label + "(R): " << nlocal
            << std::endl;
        ofs << "Matrix number of " + options.label + "(R): "
            << output_R_number << std::endl;
    }
}

void write_R_record(std::ofstream& ofs,
                    const ModuleIO::RCoordinate& R_coor,
                    const long long nonzero_count,
                    const bool binary)
{
    int dRx = R_coor.x;
    int dRy = R_coor.y;
    int dRz = R_coor.z;
    if (binary)
    {
        const int count = static_cast<int>(nonzero_count);
        ofs.write(reinterpret_cast<char*>(&dRx), sizeof(int));
        ofs.write(reinterpret_cast<char*>(&dRy), sizeof(int));
        ofs.write(reinterpret_cast<char*>(&dRz), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&count), sizeof(int));
    }
    else
    {
        ofs << dRx << " " << dRy << " " << dRz << " " << nonzero_count
            << std::endl;
    }
}

void check_output_file_open(const std::ofstream& ofs,
                            const std::string& filename,
                            const std::string& context)
{
    if (!ofs.is_open())
    {
        ModuleBase::WARNING_QUIT(context, "Cannot open sparse matrix file: " + filename);
    }
}
} // namespace

void ModuleIO::save_dH_sparse(const int& istep,
                              const Parallel_Orbitals& pv,
                              LCAO_HS_Arrays& HS_Arrays,
                              const double& sparse_thr,
                              const bool& binary,
                              const std::string& fileflag,
                              const int precision) {
    ModuleBase::TITLE("ModuleIO", "save_dH_sparse");
    ModuleBase::timer::start("ModuleIO", "save_dH_sparse");
    SparseWriteOptions single_R_options;
    single_R_options.threshold = sparse_thr;
    single_R_options.binary = binary;
    single_R_options.precision = precision;
    single_R_options.reduce = true;
    single_R_options.temp_dir = PARAM.globalv.global_out_dir;

    auto& all_R_coor_ptr = HS_Arrays.all_R_coor;
    auto& output_R_coor_ptr = HS_Arrays.output_R_coor;
    auto& dHRx_sparse_ptr = HS_Arrays.dHRx_sparse;
    auto& dHRx_soc_sparse_ptr = HS_Arrays.dHRx_soc_sparse;
    auto& dHRy_sparse_ptr = HS_Arrays.dHRy_sparse;
    auto& dHRy_soc_sparse_ptr = HS_Arrays.dHRy_soc_sparse;
    auto& dHRz_sparse_ptr = HS_Arrays.dHRz_sparse;
    auto& dHRz_soc_sparse_ptr = HS_Arrays.dHRz_soc_sparse;

    const int total_R_num = static_cast<int>(all_R_coor_ptr.size());
    int output_R_number = 0;
    std::vector<long long> dHx_nonzero_num[2];
    std::vector<long long> dHy_nonzero_num[2];
    std::vector<long long> dHz_nonzero_num[2];
    int step = istep;

    int spin_loop = 1;
    if (PARAM.inp.nspin == 2) {
        spin_loop = 2;
    }

    if (PARAM.inp.nspin != 4)
    {
        for (int ispin = 0; ispin < spin_loop; ++ispin)
        {
            dHx_nonzero_num[ispin] = count_nonzeros_by_R(dHRx_sparse_ptr[ispin], all_R_coor_ptr, sparse_thr, true);
            dHy_nonzero_num[ispin] = count_nonzeros_by_R(dHRy_sparse_ptr[ispin], all_R_coor_ptr, sparse_thr, true);
            dHz_nonzero_num[ispin] = count_nonzeros_by_R(dHRz_sparse_ptr[ispin], all_R_coor_ptr, sparse_thr, true);
        }
    }
    else
    {
        dHx_nonzero_num[0] = count_nonzeros_by_R(dHRx_soc_sparse_ptr, all_R_coor_ptr, sparse_thr, true);
        dHy_nonzero_num[0] = count_nonzeros_by_R(dHRy_soc_sparse_ptr, all_R_coor_ptr, sparse_thr, true);
        dHz_nonzero_num[0] = count_nonzeros_by_R(dHRz_soc_sparse_ptr, all_R_coor_ptr, sparse_thr, true);
    }

    const auto has_output_R = [&](const int index) {
        for (int ispin = 0; ispin < spin_loop; ++ispin)
        {
            if (dHx_nonzero_num[ispin][index] != 0
                || dHy_nonzero_num[ispin][index] != 0
                || dHz_nonzero_num[ispin][index] != 0)
            {
                return true;
            }
        }
        return false;
    };

    for (int index = 0; index < total_R_num; ++index)
    {
        if (has_output_R(index))
        {
            output_R_number++;
        }
    }

    std::stringstream sshx[2];
    std::stringstream sshy[2];
    std::stringstream sshz[2];

	if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag) 
	{
		sshx[0] << PARAM.globalv.global_matrix_dir
			<< "d"<<fileflag<<"rxs1g" << step << "_nao.csr";
		sshx[1] << PARAM.globalv.global_matrix_dir
			<< "d"<<fileflag<<"rxs2g" << step << "_nao.csr";
		sshy[0] << PARAM.globalv.global_matrix_dir
			<< "d"<<fileflag<<"rys1g" << step << "_nao.csr";
		sshy[1] << PARAM.globalv.global_matrix_dir
			<< "d"<<fileflag<<"rys2g" << step << "_nao.csr";
		sshz[0] << PARAM.globalv.global_matrix_dir
			<< "d"<<fileflag<<"rzs1g" << step << "_nao.csr";
		sshz[1] << PARAM.globalv.global_matrix_dir
			<< "d"<<fileflag<<"rzs2g" << step << "_nao.csr";
	} 
	else 
	{
		sshx[0] << PARAM.globalv.global_out_dir << "d"<<fileflag<<"rxs1_nao.csr";
        sshx[1] << PARAM.globalv.global_out_dir << "d"<<fileflag<<"rxs2_nao.csr";
        sshy[0] << PARAM.globalv.global_out_dir << "d"<<fileflag<<"rys1_nao.csr";
        sshy[1] << PARAM.globalv.global_out_dir << "d"<<fileflag<<"rys2_nao.csr";
        sshz[0] << PARAM.globalv.global_out_dir << "d"<<fileflag<<"rzs1_nao.csr";
        sshz[1] << PARAM.globalv.global_out_dir << "d"<<fileflag<<"rzs2_nao.csr";
    }
    std::ofstream g1x[2];
    std::ofstream g1y[2];
    std::ofstream g1z[2];

	if (GlobalV::DRANK == 0) 
	{
		if (binary) // binary format 
		{
			int nlocal = PARAM.globalv.nlocal;
			for (int ispin = 0; ispin < spin_loop; ++ispin) 
			{
				if (PARAM.inp.calculation == "md" && PARAM.inp.out_app_flag
						&& step) 
				{
					g1x[ispin].open(sshx[ispin].str().c_str(),
                                    std::ios::binary | std::ios::app);
                    g1y[ispin].open(sshy[ispin].str().c_str(),
                                    std::ios::binary | std::ios::app);
                    g1z[ispin].open(sshz[ispin].str().c_str(),
                                    std::ios::binary | std::ios::app);
				} 
				else 
				{
                    g1x[ispin].open(sshx[ispin].str().c_str(),std::ios::binary);
                    g1y[ispin].open(sshy[ispin].str().c_str(),std::ios::binary);
                    g1z[ispin].open(sshz[ispin].str().c_str(),std::ios::binary);
                }
                check_output_file_open(g1x[ispin], sshx[ispin].str(), "ModuleIO::save_dH_sparse");
                check_output_file_open(g1y[ispin], sshy[ispin].str(), "ModuleIO::save_dH_sparse");
                check_output_file_open(g1z[ispin], sshz[ispin].str(), "ModuleIO::save_dH_sparse");

                g1x[ispin].write(reinterpret_cast<char*>(&step), sizeof(int));
                g1x[ispin].write(reinterpret_cast<char*>(&nlocal),
                                 sizeof(int));
                g1x[ispin].write(reinterpret_cast<char*>(&output_R_number),
                                 sizeof(int));

                g1y[ispin].write(reinterpret_cast<char*>(&step), sizeof(int));
                g1y[ispin].write(reinterpret_cast<char*>(&nlocal),
                                 sizeof(int));
                g1y[ispin].write(reinterpret_cast<char*>(&output_R_number),
                                 sizeof(int));

                g1z[ispin].write(reinterpret_cast<char*>(&step), sizeof(int));
                g1z[ispin].write(reinterpret_cast<char*>(&nlocal),
                                 sizeof(int));
                g1z[ispin].write(reinterpret_cast<char*>(&output_R_number),
                                 sizeof(int));
            }
		} 
		else 
		{
			for (int ispin = 0; ispin < spin_loop; ++ispin) 
			{
				if (PARAM.inp.calculation == "md" && PARAM.inp.out_app_flag && step) 
				{
					g1x[ispin].open(sshx[ispin].str().c_str(), std::ios::app);
                    g1y[ispin].open(sshy[ispin].str().c_str(), std::ios::app);
                    g1z[ispin].open(sshz[ispin].str().c_str(), std::ios::app);
				} 
				else 
				{
					GlobalV::ofs_running << " dH/dRx data are in file: " << sshx[ispin].str() << std::endl;
					GlobalV::ofs_running << " dH/dRy data are in file: " << sshy[ispin].str() << std::endl;
					GlobalV::ofs_running << " dH/dRz data are in file: " << sshz[ispin].str() << std::endl;
                    g1x[ispin].open(sshx[ispin].str().c_str());
                    g1y[ispin].open(sshy[ispin].str().c_str());
                    g1z[ispin].open(sshz[ispin].str().c_str());
                }
                check_output_file_open(g1x[ispin], sshx[ispin].str(), "ModuleIO::save_dH_sparse");
                check_output_file_open(g1y[ispin], sshy[ispin].str(), "ModuleIO::save_dH_sparse");
                check_output_file_open(g1z[ispin], sshz[ispin].str(), "ModuleIO::save_dH_sparse");

                g1x[ispin] << "STEP: " << step << std::endl;
                g1x[ispin] << "Matrix Dimension of dHx(R): " << PARAM.globalv.nlocal
                           << std::endl;
                g1x[ispin] << "Matrix number of dHx(R): " << output_R_number
                           << std::endl;

                g1y[ispin] << "STEP: " << step << std::endl;
                g1y[ispin] << "Matrix Dimension of dHy(R): " << PARAM.globalv.nlocal
                           << std::endl;
                g1y[ispin] << "Matrix number of dHy(R): " << output_R_number
                           << std::endl;

                g1z[ispin] << "STEP: " << step << std::endl;
                g1z[ispin] << "Matrix Dimension of dHz(R): " << PARAM.globalv.nlocal
                           << std::endl;
                g1z[ispin] << "Matrix number of dHz(R): " << output_R_number
                           << std::endl;
            }
        }
    }

    output_R_coor_ptr.clear();

    int count = 0;
    for (auto& R_coor: all_R_coor_ptr) {
        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        if (!has_output_R(count))
        {
            count++;
            continue;
        }

        output_R_coor_ptr.insert(R_coor);

        if (GlobalV::DRANK == 0) {
            if (binary) {
                for (int ispin = 0; ispin < spin_loop; ++ispin) {
                    const int dHx_count = static_cast<int>(dHx_nonzero_num[ispin][count]);
                    const int dHy_count = static_cast<int>(dHy_nonzero_num[ispin][count]);
                    const int dHz_count = static_cast<int>(dHz_nonzero_num[ispin][count]);
                    g1x[ispin].write(reinterpret_cast<char*>(&dRx),
                                     sizeof(int));
                    g1x[ispin].write(reinterpret_cast<char*>(&dRy),
                                     sizeof(int));
                    g1x[ispin].write(reinterpret_cast<char*>(&dRz),
                                     sizeof(int));
                    g1x[ispin].write(reinterpret_cast<const char*>(&dHx_count),
                                     sizeof(int));

                    g1y[ispin].write(reinterpret_cast<char*>(&dRx),
                                     sizeof(int));
                    g1y[ispin].write(reinterpret_cast<char*>(&dRy),
                                     sizeof(int));
                    g1y[ispin].write(reinterpret_cast<char*>(&dRz),
                                     sizeof(int));
                    g1y[ispin].write(reinterpret_cast<const char*>(&dHy_count),
                                     sizeof(int));

                    g1z[ispin].write(reinterpret_cast<char*>(&dRx),
                                     sizeof(int));
                    g1z[ispin].write(reinterpret_cast<char*>(&dRy),
                                     sizeof(int));
                    g1z[ispin].write(reinterpret_cast<char*>(&dRz),
                                     sizeof(int));
                    g1z[ispin].write(reinterpret_cast<const char*>(&dHz_count),
                                     sizeof(int));
                }
            } else {
                for (int ispin = 0; ispin < spin_loop; ++ispin) {
                    g1x[ispin] << dRx << " " << dRy << " " << dRz << " "
                               << dHx_nonzero_num[ispin][count] << std::endl;
                    g1y[ispin] << dRx << " " << dRy << " " << dRz << " "
                               << dHy_nonzero_num[ispin][count] << std::endl;
                    g1z[ispin] << dRx << " " << dRy << " " << dRz << " "
                               << dHz_nonzero_num[ispin][count] << std::endl;
                }
            }
        }

        for (int ispin = 0; ispin < spin_loop; ++ispin) {
            if (dHx_nonzero_num[ispin][count] > 0) {
                if (PARAM.inp.nspin != 4) {
                    output_single_R(g1x[ispin],
                                    dHRx_sparse_ptr[ispin][R_coor],
                                    pv,
                                    single_R_options);
                } else {
                    output_single_R(g1x[ispin],
                                    dHRx_soc_sparse_ptr[R_coor],
                                    pv,
                                    single_R_options);
                }
            }
            if (dHy_nonzero_num[ispin][count] > 0) {
                if (PARAM.inp.nspin != 4) {
                    output_single_R(g1y[ispin],
                                    dHRy_sparse_ptr[ispin][R_coor],
                                    pv,
                                    single_R_options);
                } else {
                    output_single_R(g1y[ispin],
                                    dHRy_soc_sparse_ptr[R_coor],
                                    pv,
                                    single_R_options);
                }
            }
            if (dHz_nonzero_num[ispin][count] > 0) {
                if (PARAM.inp.nspin != 4) {
                    output_single_R(g1z[ispin],
                                    dHRz_sparse_ptr[ispin][R_coor],
                                    pv,
                                    single_R_options);
                } else {
                    output_single_R(g1z[ispin],
                                    dHRz_soc_sparse_ptr[R_coor],
                                    pv,
                                    single_R_options);
                }
            }
        }

        count++;
    }

    if (GlobalV::DRANK == 0) {
        for (int ispin = 0; ispin < spin_loop; ++ispin) {
            g1x[ispin].close();
        }
        for (int ispin = 0; ispin < spin_loop; ++ispin) {
            g1y[ispin].close();
        }
        for (int ispin = 0; ispin < spin_loop; ++ispin) {
            g1z[ispin].close();
        }
    }

    ModuleBase::timer::end("ModuleIO", "save_dH_sparse");
    return;
}

template <typename Tdata>
void ModuleIO::save_sparse(
    const SparseRMatrix<Tdata>& smat,
    const std::set<RCoordinate>& all_R_coor,
    const Parallel_Orbitals& pv,
    const SparseWriteOptions& options) {
    ModuleBase::TITLE("ModuleIO", "save_sparse");
    ModuleBase::timer::start("ModuleIO", "save_sparse");
    const int nlocal = pv.get_global_row_size();
    if (nlocal <= 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::save_sparse",
                                 "Parallel_Orbitals global row size must be positive.");
    }

    const std::vector<long long> nonzero_num
        = count_nonzeros_by_R(smat, all_R_coor, options.threshold, options.reduce);
    const int output_R_number = count_output_R(nonzero_num);
    std::ofstream ofs;
    if (!options.reduce || GlobalV::DRANK == 0)
    {
        open_sparse_file(ofs, options);
        write_sparse_header(ofs, options, nlocal, output_R_number);
    }

    int count = 0;
    for (const auto& R_coor: all_R_coor)
    {
        if (nonzero_num[count] == 0)
        {
            count++;
            continue;
        }

        if (!options.reduce || GlobalV::DRANK == 0)
        {
            write_R_record(ofs, R_coor, nonzero_num[count], options.binary);
        }

        if (smat.count(R_coor))
        {
            output_single_R(ofs, smat.at(R_coor), pv, options);
        }
        else
        {
            SparseRBlock<Tdata> empty_map;
            output_single_R(ofs, empty_map, pv, options);
        }
        ++count;
    }
    if (!options.reduce || GlobalV::DRANK == 0)
    {
        ofs.close();
    }

    ModuleBase::timer::end("ModuleIO", "save_sparse");
}

template void ModuleIO::save_sparse<double>(
    const SparseRMatrix<double>&,
    const std::set<RCoordinate>&,
    const Parallel_Orbitals&,
    const SparseWriteOptions&);

template void ModuleIO::save_sparse<std::complex<double>>(
    const SparseRMatrix<std::complex<double>>&,
    const std::set<RCoordinate>&,
    const Parallel_Orbitals&,
    const SparseWriteOptions&);
