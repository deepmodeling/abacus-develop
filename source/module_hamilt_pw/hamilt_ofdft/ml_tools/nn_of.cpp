#include "nn_of.h"

NN_OFImpl::NN_OFImpl(int nrxx, int nrxx_vali, int ninpt, int nnode, int nlayer, torch::Device device)
{
    this->nrxx = nrxx;
    this->nrxx_vali = nrxx_vali;
    this->ninpt = ninpt;
    this->nnode = nnode;
    std::cout << "nnode = " << this->nnode << std::endl;
    this->nlayer = nlayer;
    std::cout << "nlayer = " << this->nlayer << std::endl;
    this->nfc = nlayer + 1;

    this->inputs = torch::zeros({this->nrxx, this->ninpt}).to(device);
    this->F = torch::zeros({this->nrxx, 1}).to(device);
    if (nrxx_vali > 0) this->input_vali = torch::zeros({nrxx_vali, this->ninpt}).to(device);

    fc1 = register_module("fc1", torch::nn::Linear(ninpt, nnode));
    fc2 = register_module("fc2", torch::nn::Linear(nnode, nnode));
    fc3 = register_module("fc3", torch::nn::Linear(nnode, nnode));
    fc4 = register_module("fc4", torch::nn::Linear(nnode, 1));

    this->to(device);
}


torch::Tensor NN_OFImpl::forward(torch::Tensor inpt)
{
    inpt = torch::tanh(fc1->forward(inpt)); // covert data into (-1,1)
    inpt = torch::tanh(fc2->forward(inpt));
    inpt = torch::tanh(fc3->forward(inpt));
    inpt = fc4->forward(inpt); // for feg = 3

    return inpt;
}