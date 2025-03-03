#include "check_atomic_stru.h"
#include "module_base/timer.h"
#include "module_base/element_covalent_radius.h"
#include "module_base/timer.h"
#include "omp.h"
void Check_Atomic_Stru::check_atomic_stru2(UnitCell& ucell, const double& factor)
{
    // First we calculate all bond length in the structure,
    // and compare with the covalent_bond_length,
    // if there has bond length is shorter than covalent_bond_length * factor,
    // we think this structure is unreasonable.
    ModuleBase::timer::tick("Check_Atomic_Stru", "Check_Atomic_Stru");
    assert(ucell.ntype > 0);
    const int ntype = ucell.ntype;
    const double lat0 = ucell.lat0;
    bool all_pass = true;
    bool no_warning = true;
    const double warning_coef = 0.6;
    const double min_factor_coef = std::min(warning_coef, factor);

    std::stringstream errorlog;
    errorlog.setf(std::ios_base::fixed, std::ios_base::floatfield);

    std::vector<double> symbol_covalent_radiuss(ntype);
    for (int it = 0; it < ntype; it++)
    {
        std::string symbol1 = "";
        for (char ch: ucell.atoms[it].label)
        {
            if (std::isalpha(ch))
            {
                symbol1.push_back(ch);
            }
        }

        if (ModuleBase::CovalentRadius.find(symbol1) != ModuleBase::CovalentRadius.end())
        {
            symbol_covalent_radiuss[it] = ModuleBase::CovalentRadius.at(symbol1);
        }
        else
        {
            std::stringstream mess;
            mess << "Notice: symbol '" << symbol1 << "' is not an element symbol!!!! ";
            mess << "set the covalent radius to be 0." << std::endl;
            GlobalV::ofs_running << mess.str();
            std::cout << mess.str();
        }
    }
    double* latvec = new double[9];
    latvec[0] = ucell.a1.x;
    latvec[1] = ucell.a2.x;
    latvec[2] = ucell.a3.x;
    latvec[3] = ucell.a1.y;
    latvec[4] = ucell.a2.y;
    latvec[5] = ucell.a3.y;
    latvec[6] = ucell.a1.z;
    latvec[7] = ucell.a2.z;
    latvec[8] = ucell.a3.z;
    double* A = new double[27 * 3];
    std::vector<std::string> cell(27);
    std::vector<std::string> label(ntype);
    for (int i = 0; i < 27; i++)
    {
        int a = (i / 9) % 3 - 1; // 计算a的值
        int b = (i / 3) % 3 - 1; // 计算b的值
        int c = i % 3 - 1;       // 计算c的值
        A[3 * i] = a * latvec[0] + b * latvec[1] + c * latvec[2];
        A[3 * i + 1] = a * latvec[3] + b * latvec[4] + c * latvec[5];
        A[3 * i + 2] = a * latvec[6] + b * latvec[7] + c * latvec[8];
        std::ostringstream tmp_oss;
        tmp_oss << " (cell:" << std::setw(2) << a << " " << std::setw(2) << b << " " << std::setw(2) << c
                << "), distance= ";
        cell[i] = tmp_oss.str();
    }
    for (int it = 0; it < ntype; it++)
    {
        std::ostringstream tmp_oss;
        tmp_oss << std::setw(3) << ucell.atoms[it].label;
        label[it] = tmp_oss.str();
    }
    #pragma omp parallel
    {
        #pragma omp single
    for (int it1 = 0; it1 < ntype; it1++)
    {
        #pragma omp task firstprivate(it1)
        for (int ia1 = 0; ia1 < ucell.atoms[it1].na; ia1++)
        {
            const double symbol1_covalent_radius = symbol_covalent_radiuss[it1];
            double x1 = ucell.atoms[it1].taud[ia1].x;
            double y1 = ucell.atoms[it1].taud[ia1].y;
            double z1 = ucell.atoms[it1].taud[ia1].z;
            
            for (int it2 = it1; it2 < ntype; it2++)
            {
                double symbol2_covalent_radius = symbol_covalent_radiuss[it2];
                const double bohr_to_a= ModuleBase::BOHR_TO_A;
                double covalent_length = (symbol1_covalent_radius + symbol2_covalent_radius) / bohr_to_a;
                const double min_error = covalent_length * min_factor_coef / ucell.lat0;
                const double min_error_2 = min_error * min_error;
                const double factor_error = covalent_length * factor / ucell.lat0;
                for (int ia2 = ia1; ia2 < ucell.atoms[it2].na; ia2++)
                {
                    double delta_x = ucell.atoms[it2].taud[ia2].x - x1;
                    double delta_y = ucell.atoms[it2].taud[ia2].y - y1;
                    double delta_z = ucell.atoms[it2].taud[ia2].z - z1;
                    const double delta_x_lat = delta_x * latvec[0] + delta_y * latvec[1] + delta_z * latvec[2];
                    const double delta_y_lat = delta_x * latvec[3] + delta_y * latvec[4] + delta_z * latvec[5];
                    const double delta_z_lat = delta_x * latvec[6] + delta_y * latvec[7] + delta_z * latvec[8];

                    double* current_ptr = &A[0];
                    const bool is_same_atom = (it1 == it2) && (ia1 == ia2);
                        
                        for (int i = 0; i < 27; i++)
                        {
                            if (is_same_atom && i==13)
                                continue;
                            const double part1 = delta_x_lat + *(current_ptr++);
                            const double part2 = delta_y_lat + *(current_ptr++);
                            const double part3 = delta_z_lat + *(current_ptr++);
                            const double bond_length = part1 * part1 + part2 * part2 + part3 * part3;
                            if (bond_length < min_error_2)
                            {
                                const double sqrt_bon = sqrt(bond_length) / lat0;
                                no_warning = false;
                                all_pass   = bond_length < factor_error ? false :true;
                                // errorlog << std::setw(3) << ia1 + 1 << "-th " << &label[it1] << ", " << std::setw(3)
                                //         << ia2 + 1 << "-th " << &label[it2] << &cell[i]
                                //         << std::setprecision(3) << sqrt_bon
                                //         << " Bohr (" << sqrt_bon * bohr_to_a << " Angstrom)\n";
                            }
                        } // a

                } // ia2
            } // it2
        } // ia1
    } // it1
    }
    delete[] latvec;
    ModuleBase::timer::tick("Check_Atomic_Stru", "Check_Atomic_Stru");
    if (!all_pass || !no_warning)
    {
        std::stringstream mess;
        mess << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "%%%%%% WARNING  WARNING  WARNING  WARNING  WARNING  %%%%%%" << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "!!! WARNING: Some atoms are too close!!!" << std::endl;
        mess << "!!! Please check the nearest-neighbor list in log file." << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "%%%%%% WARNING  WARNING  WARNING  WARNING  WARNING  %%%%%%" << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

        GlobalV::ofs_running << mess.str() << mess.str() << mess.str() << errorlog.str();
        std::cout << mess.str() << mess.str() << mess.str() << std::endl;
        std::ofstream outFile("example2.txt");
        if (outFile.is_open())
        {

            // 将stringstream中的内容写入到文件
            outFile << errorlog.str();

            // 关闭文件
            outFile.close();
            std::cout << "数据成功写入到文件。\n";
        }
        else
        {
            std::cerr << "无法打开文件进行写入。\n";
        }
        if (!all_pass)
        {
            mess.clear();
            mess.str("");
            mess << "If this structure is what you want, you can set "
                    "'min_dist_coef'"
                 << std::endl;
            mess << "as a smaller value (the current value is " << factor << ") in INPUT file." << std::endl;
            GlobalV::ofs_running << mess.str();
            std::cout << mess.str();
            // ModuleBase::WARNING_QUIT("Input", "The structure is unreasonable!");
        }
    }
}


void matrix_multiply(const std::vector<std::vector<double>>& A,
                     const std::vector<std::vector<double>>& B,
                     std::vector<std::vector<double>>& C,
                     int N) {
    // 使用OpenMP并行化外层循环，即对C中的每一行进行并行处理
    ModuleBase::timer::tick("Check_Atomic_Stru", "matrix_multiply");
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int k = 0; k < N; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    ModuleBase::timer::tick("Check_Atomic_Stru", "matrix_multiply");
}

void test()
{
    ModuleBase::timer::tick("Check_Atomic_Stru", "test");
    std::stringstream errorlog;
    const double delta_x_lat=10;
    const double delta_y_lat=10;
    const double delta_z_lat=10;
    const double min_error_2=1.0;
    const double bohr_to_a=1.0;
    const double lat0 =1.0;
    double *A =new double [27*3];
    bool no_warning=0;
    bool all_pass=0;
    std::vector<std::stringstream> local_results(omp_get_num_threads());
    #pragma omp parallel
    {
        int tid=omp_get_thread_num();
        std::vector<std::string> local_result;
    #pragma omp parallel for
    for (int j=0;j<1024*1024;j++)
    for (int i = 0; i < 27; i++)
    {
        // if (is_same_atom && i==13)
        //     continue;
        const double part1 = delta_x_lat + A[3*i];
        const double part2 = delta_y_lat + A[3*i+1];
        const double part3 = delta_z_lat + A[3*i+2];
        const double bond_length = part1 * part1 + part2 * part2 + part3 * part3;
        // if (bond_length < min_error_2)
        // {
            const double sqrt_bon = sqrt(bond_length) / lat0;
            no_warning = false;
            all_pass   = bond_length < min_error_2 ? false :true;
            // errorlog << std::setw(3) << ia1 + 1 << "-th " << &label[it1] << ", " << std::setw(3)
            //         << ia2 + 1 << "-th " << &label[it2] << &cell[i]
                // local_results[tid] << std::setprecision(3) << sqrt_bon
                //     << " Bohr (" << sqrt_bon * bohr_to_a << " Angstrom)\n";
                local_result.push_back( " Bohr (");

        // }
    }
    }
    ModuleBase::timer::tick("Check_Atomic_Stru", "test");
}


void process_data() {
    ModuleBase::timer::tick("Check_Atomic_Stru", "process_data");
    std::vector<std::string> local_results;
    const size_t total_iterations = 1024 * 1024 *1024;
    const double lat0 = 1.0;
    const double bohr_to_a = 0.529177; // 示例值
    const double min_error_2 = 1.0; // 示例值

    std::vector<std::string> results;
    results.reserve(total_iterations); // 预分配空间以提高性能

    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0) {
            local_results.clear(); // 初始化线程局部存储
        }

        #pragma omp barrier // 确保所有线程都已初始化它们的局部存储

        #pragma omp for
        for (size_t i = 0; i < total_iterations; ++i) {
            const double delta_x_lat = sqrt(i + i);
            const double delta_y_lat = sqrt(i * i);
            const double delta_z_lat = sqrt(i * i + 10);
            const double sqrt_bon = sqrt(delta_x_lat + delta_y_lat + delta_z_lat) / lat0;

            if (sqrt_bon < min_error_2) {
                std::stringstream ss;
                ss << std::setprecision(3) << sqrt_bon
                   << " Bohr (" << sqrt_bon * bohr_to_a << " Angstrom)\n";
                local_results.push_back(ss.str());
            }
        }

        // 合并线程局部存储到全局存储
        #pragma omp critical
        {
            results.insert(results.end(), local_results.begin(), local_results.end());
        }
    }

    // 统一输出结果
    for (const auto& result : results) {
        std::cerr << result;
    }
    ModuleBase::timer::tick("Check_Atomic_Stru", "process_data");
}

// 注意：这里假设 local_results 已经在某个地方定义了
std::vector<std::string> local_results;

void Check_Atomic_Stru::check_atomic_stru1(UnitCell& ucell, const double& factor)
{
    const double lat0=1.0;
    ModuleBase::timer::tick("Check_Atomic_Stru", "Check_Atomic_Stru1");
    const int N = 1000; // 矩阵大小

    // 初始化矩阵A和B
    std::vector<std::vector<double>> A(N, std::vector<double>(N));
    std::vector<std::vector<double>> B(N, std::vector<double>(N));
    std::vector<std::vector<double>> C(N, std::vector<double>(N));
    // 执行矩阵乘法
    matrix_multiply(A, B, C, N);
    // test();
    process_data();
    ModuleBase::timer::tick("Check_Atomic_Stru", "Check_Atomic_Stru1");
}