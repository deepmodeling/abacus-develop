#include "pos_op_writer.h"

#include "source_base/global_variable.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_io/module_hs/lat_r_csr.h"
#include "source_io/module_hs/pos_op_csr.h"
#include "source_io/module_output/write_notice.h"
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
                             int out_app_flag,
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

    // One temporary payload per Cartesian component; the final files are
    // assembled after the R loop, once the total block count is known.
    const char axis_label[3] = {'x', 'y', 'z'};
    std::stringstream tem[3];
    std::ofstream ofs_tem[3];

    if (GlobalV::DRANK == 0)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            tem[direction] << global_out_dir << "tmp-rr" << axis_label[direction] << "_nao.txt";
            std::ios_base::openmode mode = std::ios::out;
            if (binary)
            {
                mode |= std::ios::binary;
            }
            const std::string temp_name = tem[direction].str();
            ofs_tem[direction].open(temp_name.c_str(), mode);
            if (!ofs_tem[direction].is_open())
            {
                const std::string message = "Cannot open temporary sparse matrix file: " + temp_name;
                ModuleBase::WARNING_QUIT("PosOpWriter::out_lat_r", message);
            }
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
            ir = pv_.para_2d().irloc(iw1);
            if (ir >= 0)
            {
                for (int iw2 = 0; iw2 < nlocal; iw2++)
                {
                    ic = pv_.para_2d().icloc(iw2);
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

            for (int direction = 0; direction < 3; ++direction)
            {
                const int nnz = rR_nonzero_num[direction];
                if (GlobalV::DRANK == 0)
                {
                    if (binary)
                    {
                        ofs_tem[direction].write(reinterpret_cast<char*>(&dRx), sizeof(int));
                        ofs_tem[direction].write(reinterpret_cast<char*>(&dRy), sizeof(int));
                        ofs_tem[direction].write(reinterpret_cast<char*>(&dRz), sizeof(int));
                        ofs_tem[direction].write(reinterpret_cast<char*>(&rR_nonzero_num[direction]),
                                                 sizeof(int));
                    }
                    else
                    {
                        ofs_tem[direction] << dRx << " " << dRy << " " << dRz << " " << nnz << std::endl;
                    }
                }

                if (binary)
                {
                    if (nnz > 0)
                    {
                        ModuleIO::save_lat_r(ofs_tem[direction],
                                             psi_r_psi_sparse[direction],
                                             pv_,
                                             lat_r_options);
                    }
                }
                else
                {
                    // Text output always emits the three CSR sections so that
                    // an empty component block remains a valid zero matrix record.
                    ModuleIO::save_lat_r_csr(ofs_tem[direction],
                                             psi_r_psi_sparse[direction],
                                             pv_,
                                             lat_r_options);
                }
            }
        }
    }

    // Append from the second ionic step onward in append mode; the first
    // step (step == 0) truncates so a rerun does not retain prior frames.
    // Single-point SCF passes step = -1 and also starts a fresh file.
    const bool open_in_append = (out_app_flag == 1) && (step > 0);
    // ionic step printed one-based per governance rule 22
    const int header_step = std::max(step, 0) + 1;

    if (GlobalV::DRANK == 0)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            // Modes 0 and 1 keep the fixed filename (mode 0 overwrites, mode 1
            // appends all ionic steps); only mode 2 embeds the one-based step.
            std::stringstream ssr;
            if (step >= 0 && out_app_flag == 2)
            {
                ssr << global_matrix_dir << "rr" << axis_label[direction] << "g"
                    << (step + 1) << "_nao.csr";
            }
            else
            {
                ssr << global_matrix_dir << "rr" << axis_label[direction] << "_nao.csr";
            }

            const std::string output_name = ssr.str();
            const std::string payload_name = tem[direction].str();
            const std::string axis(1, axis_label[direction]);

            ofs_tem[direction].close();
            ModuleIO::detail::assemble_csr(output_name,
                                           payload_name,
                                           header_step,
                                           nlocal,
                                           output_R_number,
                                           binary,
                                           open_in_append,
                                           axis,
                                           &ucell,
                                           "PosOpWriter::out_lat_r");
            ModuleIO::write_matrix_notice(ofs_running, "r(R)", output_name);

            std::remove(payload_name.c_str());
        }
    }

    ModuleBase::timer::end("PosOpWriter", "out_lat_r");
    return;
}
