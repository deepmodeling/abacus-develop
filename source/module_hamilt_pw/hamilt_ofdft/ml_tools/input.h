#ifndef INPUT_H
#define INPUT_H

#include <torch/torch.h>

class Input
{
    // ---------- read in the settings from nnINPUT --------
  public:
    Input(){};
    ~Input()
    {
        delete[] this->train_dir;
        delete[] this->train_cell;
        delete[] this->train_a;
        delete[] this->validation_dir;
        delete[] this->validation_cell;
        delete[] this->validation_a;

        delete[] this->ml_gammanl;
        delete[] this->ml_pnl;
        delete[] this->ml_qnl;
        delete[] this->ml_xi;
        delete[] this->ml_tanhxi;
        delete[] this->ml_tanhxi_nl;
        delete[] this->ml_tanh_pnl;
        delete[] this->ml_tanh_qnl;
        delete[] this->ml_tanhp_nl;
        delete[] this->ml_tanhq_nl;
        delete[] this->chi_xi;
        delete[] this->chi_pnl;
        delete[] this->chi_qnl;
        delete[] this->kernel_type;
        delete[] this->kernel_scaling;
        delete[] this->yukawa_alpha;
        delete[] this->kernel_file;

    };

    void readInput();

    template <class T> static void read_value(std::ifstream &ifs, T &var)
    {
        ifs >> var;
        ifs.ignore(150, '\n');
        return;
    }

    template <class T> static void read_values(std::ifstream &ifs, const int length, T *var)
    {
        for (int i = 0; i < length; ++i)
        {
            ifs >> var[i];
        }
        ifs.ignore(150, '\n');
        return;
    }

    // training
    int fftdim = 0;
    int nbatch = 0;
    int ntrain = 1;
    int nvalidation = 0;
    std::string *train_dir = nullptr;
    std::string *train_cell = nullptr;
    double *train_a = nullptr;
    std::string *validation_dir = nullptr;
    std::string *validation_cell = nullptr;
    double *validation_a = nullptr;
    std::string loss = "both";
    int nepoch = 1000;
    double lr_start = 0.01; // learning rate 2023-02-24
    double lr_end = 1e-4;
    int lr_fre = 5000;
    double exponent = 5.; // exponent of weight rho^{exponent/3.}

    // output
    int dump_fre = 1;
    int print_fre = 1;

    // descriptors
    // semi-local descriptors
    bool ml_gamma = false;
    bool ml_p = false;
    bool ml_q = false;
    bool ml_tanhp = false;
    bool ml_tanhq = false;
    double chi_p = 1.;
    double chi_q = 1.;
    // non-local descriptors
    bool* ml_gammanl = nullptr;
    bool* ml_pnl = nullptr;
    bool* ml_qnl = nullptr;
    bool* ml_xi = nullptr;
    bool* ml_tanhxi = nullptr;
    bool* ml_tanhxi_nl = nullptr;
    bool* ml_tanh_pnl = nullptr;
    bool* ml_tanh_qnl = nullptr;
    bool* ml_tanhp_nl = nullptr;
    bool* ml_tanhq_nl = nullptr;
    double* chi_xi = nullptr;
    double* chi_pnl = nullptr;
    double* chi_qnl = nullptr;

    int feg_limit = 0; // Free Electron Gas
    int change_step = 0; // when feg_limit=3, change the output of net after change_step

    // coefficients in loss function
    double coef_e = 1.;
    double coef_p = 1.;
    double coef_feg_e = 1.;
    double coef_feg_p = 1.;

    // size of nn
    int nnode = 10;
    int nlayer = 3;

    // kernel
    int nkernel = 1;
    int* kernel_type = nullptr;
    double* kernel_scaling = nullptr;
    double* yukawa_alpha = nullptr;
    std::string* kernel_file = nullptr;

    // GPU
    std::string device_type = "gpu";
    bool check_pot = false;

    static void print(std::string message)
    {
        std::cout << message << std::endl;
    }
};
#endif