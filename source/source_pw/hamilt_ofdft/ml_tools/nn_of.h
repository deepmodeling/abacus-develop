#ifndef NN_OF_H
#define NN_OF_H

#include <torch/torch.h>

struct NN_OFImpl:torch::nn::Module{
    // three hidden layers and one output layer
    NN_OFImpl(
        int nrxx, 
        int nrxx_vali, 
        int ninpt, 
        int nnode,
        int nlayer,
        torch::Device device
        );
    ~NN_OFImpl()
    {
        // delete[] this->fcs;
    };


    template <class T>
    void set_data(
        T *data,
        const std::vector<std::string> &descriptor_type,
        const std::vector<int> &kernel_index,
        torch::Tensor &nn_input
    )
    {
        if (data->nx_tot <= 0) return;
        for (int i = 0; i < descriptor_type.size(); ++i)
        {
            nn_input.index({"...", i}) = data->get_data(descriptor_type[i], kernel_index[i]);
        }
    }

    torch::Tensor forward(torch::Tensor inpt);

    // torch::nn::Linear fc1{nullptr}, fc2{nullptr}, fc3{nullptr}, fc4{nullptr}, fc5{nullptr};
    // torch::nn::Linear fcs[5] = {fc1, fc2, fc3, fc4, fc5};

    torch::nn::Linear fc1{nullptr}, fc2{nullptr}, fc3{nullptr}, fc4{nullptr};

    torch::Tensor inputs;
    torch::Tensor input_vali;
    torch::Tensor F; // enhancement factor, output of NN

    int nrxx = 10;
    int nrxx_vali = 0;
    int ninpt = 6;
    int nnode = 10;
    int nlayer = 3;
    int nfc = 4;
};
TORCH_MODULE(NN_OF);

#endif