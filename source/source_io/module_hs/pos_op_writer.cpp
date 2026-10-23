#include "pos_op_writer.h"

#include "source_base/global_variable.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_io/module_hs/lat_r_csr.h"
#include "source_io/module_hs/pos_op_csr.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"

#include <cstdio>
#include <sstream>

PosOpWriter::PosOpWriter(const PosOpBasis& basis, const PosOpCalc& calc, const Parallel_Orbitals& pv)
    : basis_(basis), calc_(calc), pv_(pv)
{
}

void PosOpWriter::out_lat_r(const UnitCell& ucell,
                             const Grid_Driver& gd,
                             int istep,
                             int precision,
                             const std::string& global_out_dir,
                             const std::string& global_matrix_dir,
                             const std::string& calculation,
                             bool out_app_flag,
                             int nlocal,
                             int npol,
                             double sparse_threshold,
                             bool binary,
                             std::ofstream& ofs_running)
{
    ModuleBase::TITLE("PosOpWriter", "out_lat_r");
    ModuleBase::timer::start("PosOpWriter", "out_lat_r");

    int step = istep;
    // set R coor range
    int R_minX = int(-gd.getGlayerX_minus());
    int R_minY = int(-gd.getGlayerY_minus());
    int R_minZ = int(-gd.getGlayerZ_minus());

    int R_x = gd.getGlayerX() + gd.getGlayerX_minus();
    int R_y = gd.getGlayerY() + gd.getGlayerY_minus();
    int R_z = gd.getGlayerZ() + gd.getGlayerZ_minus();

    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    for (int ix = 0; ix < R_x; ix++)
    {
        for (int iy = 0; iy < R_y; iy++)
        {
            for (int iz = 0; iz < R_z; iz++)
            {
                Abfs::Vector3_Order<int> temp_R(ix + R_minX, iy + R_minY, iz + R_minZ);
                all_R_coor.insert(temp_R);
            }
        }
    }

    // calculate rR matrix
    ModuleBase::Vector3<double> origin_point(0.0, 0.0, 0.0);
    double factor = sqrt(ModuleBase::FOUR_PI / 3.0);
    int output_R_number = 0;
    ModuleIO::SparseWriteOptions lat_r_options;
    lat_r_options.threshold = sparse_threshold;
    lat_r_options.binary = binary;
    lat_r_options.precision = precision;
    lat_r_options.reduce = true;

    std::stringstream tem1;
    tem1 << global_out_dir << "tmp-rr_nao.txt";
    std::ofstream ofs_tem1;

    if (GlobalV::DRANK == 0)
    {
        if (binary)
        {
            ofs_tem1.open(tem1.str().c_str(), std::ios::binary);
        }
        else
        {
            ofs_tem1.open(tem1.str().c_str());
        }
        if (!ofs_tem1.is_open())
        {
            ModuleBase::WARNING_QUIT("PosOpWriter::out_lat_r", "Cannot open temporary sparse matrix file: " + tem1.str());
        }
    }

    for (auto& R_coor: all_R_coor)
    {
        std::map<size_t, std::map<size_t, double>> psi_r_psi_sparse[3];

        int dRx = R_coor.x;
        int dRy = R_coor.y;
        int dRz = R_coor.z;

        ModuleBase::Vector3<double> R_car = ModuleBase::Vector3<double>(dRx, dRy, dRz) * ucell.latvec;

        int ir, ic;
        for (int iw1 = 0; iw1 < nlocal; iw1++)
        {
            ir = pv_.global2local_row(iw1);
            if (ir >= 0)
            {
                for (int iw2 = 0; iw2 < nlocal; iw2++)
                {
                    ic = pv_.global2local_col(iw2);
                    if (ic >= 0)
                    {
                        int orb_index_row = iw1 / npol;
                        int orb_index_col = iw2 / npol;

                        // The off-diagonal term in SOC calculaiton is zero, and the two diagonal terms are the same
                        int new_index
                            = iw1 - npol * orb_index_row + (iw2 - npol * orb_index_col) * npol;

                        if (new_index == 0 || new_index == 3)
                        {
                            int it1 = basis_.get_iw2it(orb_index_row);
                            int ia1 = basis_.get_iw2ia(orb_index_row);
                            int iN1 = basis_.get_iw2iN(orb_index_row);
                            int iL1 = basis_.get_iw2iL(orb_index_row);
                            int im1 = basis_.get_iw2im(orb_index_row);

                            int it2 = basis_.get_iw2it(orb_index_col);
                            int ia2 = basis_.get_iw2ia(orb_index_col);
                            int iN2 = basis_.get_iw2iN(orb_index_col);
                            int iL2 = basis_.get_iw2iL(orb_index_col);
                            int im2 = basis_.get_iw2im(orb_index_col);

                            // PosOpCalc::pos_matrix expects absolute Cartesian positions
                            // of both centers and computes the inter-center distance itself.
                            // The second center is atom ia2 translated by the lattice vector R.
                            ModuleBase::Vector3<double> tau1_car_ = ucell.atoms[it1].tau[ia1] * ucell.lat0;
                            ModuleBase::Vector3<double> tau2_car_ = (ucell.atoms[it2].tau[ia2] + R_car) * ucell.lat0;

                            ModuleBase::Vector3<double> temp_prp = calc_.pos_matrix(tau1_car_,
                                                                                     it1, iL1, im1, iN1,
                                                                                     tau2_car_,
                                                                                     it2, iL2, im2, iN2);

                            if (std::abs(temp_prp.x) > sparse_threshold)
                            {
                                psi_r_psi_sparse[0][iw1][iw2] = temp_prp.x;
                            }

                            if (std::abs(temp_prp.y) > sparse_threshold)
                            {
                                psi_r_psi_sparse[1][iw1][iw2] = temp_prp.y;
                            }

                            if (std::abs(temp_prp.z) > sparse_threshold)
                            {
                                psi_r_psi_sparse[2][iw1][iw2] = temp_prp.z;
                            }
                        }
                    }
                }
            }
        }

        int rR_nonzero_num[3] = {0, 0, 0};
        for (int direction = 0; direction < 3; ++direction)
        {
            for (auto& row_loop: psi_r_psi_sparse[direction])
            {
                rR_nonzero_num[direction] += row_loop.second.size();
            }
        }

        Parallel_Reduce::reduce_all(rR_nonzero_num, 3);

        if (ModuleIO::detail::lat_r_nonempty(rR_nonzero_num))
        {
            output_R_number++;

            if (GlobalV::DRANK == 0)
            {
                if (binary)
                {
                    ofs_tem1.write(reinterpret_cast<char*>(&dRx), sizeof(int));
                    ofs_tem1.write(reinterpret_cast<char*>(&dRy), sizeof(int));
                    ofs_tem1.write(reinterpret_cast<char*>(&dRz), sizeof(int));
                }
                else
                {
                    ofs_tem1 << dRx << " " << dRy << " " << dRz << std::endl;
                }
            }

            for (int direction = 0; direction < 3; ++direction)
            {
                if (GlobalV::DRANK == 0)
                {
                    if (binary)
                    {
                        ofs_tem1.write(reinterpret_cast<char*>(&rR_nonzero_num[direction]), sizeof(int));
                    }
                    else
                    {
                        ofs_tem1 << rR_nonzero_num[direction] << std::endl;
                    }
                }

                if (rR_nonzero_num[direction])
                {
                    ModuleIO::save_lat_r(ofs_tem1, psi_r_psi_sparse[direction], pv_, lat_r_options);
                }
                else
                {
                    // do nothing
                }
            }
        }
    }

    if (GlobalV::DRANK == 0)
    {
        std::stringstream ssr;
        const bool md_no_append = (calculation == "md") && !out_app_flag;
        if (step >= 0)
        {
            ssr << (md_no_append ? global_matrix_dir : global_out_dir) << "rrg" << (step + 1) << "_nao.txt";
        }
        else
        {
            ssr << global_out_dir << "rr_nao.txt";
        }

        ofs_tem1.close();
        const bool open_in_append = (calculation == "md") && out_app_flag && step >= 0;
        const int header_step = std::max(step, 0);
        ModuleIO::detail::assemble_csr(ssr.str(),
                                       tem1.str(),
                                       header_step,
                                       nlocal,
                                       output_R_number,
                                       binary,
                                       open_in_append,
                                       "PosOpWriter::out_lat_r");
        ofs_running << " Write r(R) matrix in NAO basis to file: " << ssr.str() << std::endl;

        std::remove(tem1.str().c_str());
    }

    ModuleBase::timer::end("PosOpWriter", "out_lat_r");
    return;
}
