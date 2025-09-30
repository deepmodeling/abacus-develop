#ifndef TRAIN_KEDF_H
#define TRAIN_KEDF_H

#include "./data.h"
#include "./grid.h"
#include "./input.h"
#include "./kernel.h"
#include "./nn_of.h"
#include "./pauli_potential.h"

#include <torch/torch.h>

class Train_KEDF
{
  public:
    Train_KEDF(){};
    ~Train_KEDF();

    std::shared_ptr<NN_OFImpl> nn;
    Input input;
    Grid grid_train;
    Grid grid_vali;
    Kernel *kernel_train = nullptr;
    Kernel *kernel_vali = nullptr;
    PauliPotential potential;
    //----------- training set -----------
    Data data_train;
    double *train_volume = nullptr;
    //---------validation set ------------
    Data data_vali;
    double *vali_volume = nullptr;
    // ------------------------------------

    torch::Device device = torch::Device(torch::kCUDA);
    int ninput = 0;
    std::vector<std::string> descriptor_type = {};
    std::vector<int> kernel_index = {};

    // -------- free electron gas ---------
    torch::Tensor feg_inpt;
    torch::Tensor feg_predict;
    torch::Tensor feg_dFdgamma;

    // ----------- constants ---------------
    double feg3_correct = 0.541324854612918; // ln(e - 1)
    const double cTF
        = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)
          * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    const double pqcoef = 1.0 / (4.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0)); // coefficient of p and q

    void train();
    void potTest();
    void setUpFFT();
    void set_device();
    void init_input_index();
    void init();

    torch::Tensor lossFunction(torch::Tensor enhancement, torch::Tensor target, torch::Tensor coef = torch::ones(1));
    torch::Tensor lossFunction_new(torch::Tensor enhancement,
                                   torch::Tensor target,
                                   torch::Tensor weight,
                                   torch::Tensor coef = torch::ones(1));
};

// class OF_data : public torch::data::Dataset<OF_data>
// {
// private:
//     torch::Tensor input;
//     torch::Tensor target;

// public:
//     explicit OF_data(torch::Tensor &input, torch::Tensor &target)
//     {
//         this->input = input.clone();
//         this->target = target.clone();
//     }

//     torch::data::Example<> get(size_t index) override
//     {
//         return {this->input[index], this->target[index]};
//     }

//     torch::optional<size_t> size() const override
//     {
//         return this->input.size(0);
//     }
// };

#endif