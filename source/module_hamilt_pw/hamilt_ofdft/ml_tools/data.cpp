#include "./data.h"
#include "npy.hpp"

Data::~Data()
{
    delete[] load_gammanl;
    delete[] load_pnl;
    delete[] load_qnl;
    delete[] load_xi;
    delete[] load_tanhxi;
    delete[] load_tanhxi_nl;
    delete[] load_tanh_pnl;
    delete[] load_tanh_qnl;
    delete[] load_tanhp_nl;
    delete[] load_tanhq_nl;
}

void Data::load_data(Input &input, const int ndata, std::string *dir, const torch::Device device)
{
    if (ndata <= 0) { return;
}
    this->init_label(input);
    this->init_data(input.nkernel, ndata, input.fftdim, device);
    this->load_data_(input, ndata, input.fftdim, dir);
    std::cout << "enhancement mean: " << this->enhancement_mean << std::endl;
    std::cout << "exponent: " << input.exponent << std::endl;
    std::cout << "tau mean: " << this->tau_mean << std::endl;
    std::cout << "pauli potential mean: " << this->pauli_mean << std::endl;
    std::cout << "Load data done" << std::endl;
}

torch::Tensor Data::get_data(std::string parameter, const int ikernel){
    if (parameter == "gamma"){
        return this->gamma.reshape({this->nx_tot});
    }
    if (parameter == "p"){
        return this->p.reshape({this->nx_tot});
    }
    if (parameter == "q"){
        return this->q.reshape({this->nx_tot});
    }
    if (parameter == "tanhp"){
        return this->tanhp.reshape({this->nx_tot});
    }
    if (parameter == "tanhq"){
        return this->tanhq.reshape({this->nx_tot});
    }
    if (parameter == "gammanl"){
        return this->gammanl[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "pnl"){
        return this->pnl[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "qnl"){
        return this->qnl[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "xi"){
        return this->xi[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "tanhxi"){
        return this->tanhxi[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "tanhxi_nl"){
        return this->tanhxi_nl[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "tanh_pnl"){
        return this->tanh_pnl[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "tanh_qnl"){
        return this->tanh_qnl[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "tanhp_nl"){
        return this->tanhp_nl[ikernel].reshape({this->nx_tot});
    }
    if (parameter == "tanhq_nl"){
        return this->tanhq_nl[ikernel].reshape({this->nx_tot});
    }
    return torch::zeros({});
}

void Data::init_label(Input &input)
{
    // Input::print("init_label");
    this->load_gammanl = new bool[input.nkernel];
    this->load_pnl = new bool[input.nkernel];
    this->load_qnl = new bool[input.nkernel];
    this->load_xi = new bool[input.nkernel];
    this->load_tanhxi = new bool[input.nkernel];
    this->load_tanhxi_nl = new bool[input.nkernel];
    this->load_tanh_pnl = new bool[input.nkernel];
    this->load_tanh_qnl = new bool[input.nkernel];
    this->load_tanhp_nl = new bool[input.nkernel];
    this->load_tanhq_nl = new bool[input.nkernel];

    bool load_gammanl_tot = false;
    bool load_pnl_tot = false;
    bool load_qnl_tot = false;
    // bool load_xi_tot = false;
    // bool load_tanhxi_tot = false;
    // bool load_tanhxi_nl_tot = false;
    bool load_tanh_pnl_tot = false;
    bool load_tanh_qnl_tot = false;
    bool load_tanhp_nl_tot = false;
    bool load_tanhq_nl_tot = false;

    for (int ik = 0; ik < input.nkernel; ++ik)
    {
        this->load_gammanl[ik] = input.ml_gammanl[ik];
        this->load_pnl[ik] = input.ml_pnl[ik];
        this->load_qnl[ik] = input.ml_qnl[ik];
        this->load_tanhxi_nl[ik] = input.ml_tanhxi_nl[ik];
        this->load_tanhxi[ik] = input.ml_tanhxi[ik] || input.ml_tanhxi_nl[ik];
        this->load_xi[ik] = input.ml_xi[ik] || this->load_tanhxi[ik];
        this->load_tanh_pnl[ik] = input.ml_tanh_pnl[ik];
        this->load_tanh_qnl[ik] = input.ml_tanh_qnl[ik];
        this->load_tanhp_nl[ik] = input.ml_tanhp_nl[ik];
        this->load_tanhq_nl[ik] = input.ml_tanhq_nl[ik];
        // this->load_pnl[ik] = input.ml_pnl[ik] || input.ml_tanh_pnl[ik];

        load_gammanl_tot = load_gammanl_tot || this->load_gammanl[ik];
        load_pnl_tot = load_pnl_tot || this->load_pnl[ik];
        load_qnl_tot = load_qnl_tot || this->load_qnl[ik];
        // load_xi_tot = load_xi_tot || this->load_xi[ik];
        // load_tanhxi_tot = load_tanhxi_tot || this->load_tanhxi[ik];
        // load_tanhxi_nl_tot = load_tanhxi_nl_tot || this->load_tanhxi_nl[ik];
        load_tanh_pnl_tot = load_tanh_pnl_tot || this->load_tanh_pnl[ik];
        load_tanh_qnl_tot = load_tanh_qnl_tot || this->load_tanh_qnl[ik];
        load_tanhp_nl_tot = load_tanhp_nl_tot || this->load_tanhp_nl[ik];
        load_tanhq_nl_tot = load_tanhq_nl_tot || this->load_tanhq_nl[ik];
        
        // std::cout << "load_gammanl    " << this->load_gammanl[ik] << std::endl;
        // std::cout << "load_pnl    " << this->load_pnl[ik] << std::endl;
        // std::cout << "load_qnl    " << this->load_qnl[ik] << std::endl;
        // std::cout << "load_tanhxi_nl    " << this->load_tanhxi_nl[ik] << std::endl;
        // std::cout << "load_tanhxi    " << this->load_tanhxi[ik] << std::endl;
        // std::cout << "load_xi    " << this->load_xi[ik] << std::endl;
        // std::cout << "load_tanh_pnl    " << this->load_tanh_pnl[ik] << std::endl;
        // std::cout << "load_tanh_qnl    " << this->load_tanh_qnl[ik] << std::endl;
        // std::cout << "load_tanhp_nl    " << this->load_tanhp_nl[ik] << std::endl;
        // std::cout << "load_tanhq_nl    " << this->load_tanhq_nl[ik] << std::endl;
    }
    this->load_gamma = input.ml_gamma || load_gammanl_tot;
    this->load_tanhp = input.ml_tanhp || load_tanhp_nl_tot || load_tanh_pnl_tot;
    this->load_tanhq = input.ml_tanhq || load_tanhq_nl_tot || load_tanh_qnl_tot;
    this->load_p = input.ml_p || this->load_tanhp || load_pnl_tot;
    this->load_q = input.ml_q || this->load_tanhq || load_qnl_tot;
    // Input::print("init_label done");
}

void Data::init_data(const int nkernel, const int ndata, const int fftdim, const torch::Device device)
{
    // Input::print("init_data");
    this->nx = std::pow(fftdim, 3);
    this->nx_tot = this->nx * ndata;

    this->rho = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
    if (this->load_p){
        this->nablaRho = torch::zeros({ndata, 3, fftdim, fftdim, fftdim}).to(device);
    }

    this->enhancement = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
    this->enhancement_mean = torch::zeros(ndata).to(device);
    this->tau_mean = torch::zeros(ndata).to(device);
    this->pauli = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
    this->pauli_mean = torch::zeros(ndata).to(device);

    if (this->load_gamma){
        this->gamma = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
    }
    if (this->load_p){
        this->p = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
    }
    if (this->load_q){
        this->q = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
    }
    if (this->load_tanhp){
        this->tanhp = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
    }
    if (this->load_tanhq){
        this->tanhq = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
    }

    for (int ik = 0; ik < nkernel; ++ik)
    {
        this->gammanl.push_back(torch::zeros({}).to(device));
        this->pnl.push_back(torch::zeros({}).to(device));
        this->qnl.push_back(torch::zeros({}).to(device));
        this->xi.push_back(torch::zeros({}).to(device));
        this->tanhxi.push_back(torch::zeros({}).to(device));
        this->tanhxi_nl.push_back(torch::zeros({}).to(device));
        this->tanh_pnl.push_back(torch::zeros({}).to(device));
        this->tanh_qnl.push_back(torch::zeros({}).to(device));
        this->tanhp_nl.push_back(torch::zeros({}).to(device));
        this->tanhq_nl.push_back(torch::zeros({}).to(device));

        if (this->load_gammanl[ik]){
            this->gammanl[ik]   = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_pnl[ik]){
            this->pnl[ik]       = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_qnl[ik]){
            this->qnl[ik]       = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_xi[ik]){
            this->xi[ik]        = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_tanhxi[ik]){
            this->tanhxi[ik]    = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_tanhxi_nl[ik{
            this->tanhxi_nl[ik] = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_tanh_pnl[ik]){
            this->tanh_pnl[ik]  = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_tanh_qnl[ik]){
            this->tanh_qnl[ik]  = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_tanhp_nl[ik]){
            this->tanhp_nl[ik]  = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
        if (this->load_tanhq_nl[ik]){
            this->tanhq_nl[ik]  = torch::zeros({ndata, fftdim, fftdim, fftdim}).to(device);
        }
    }

    // Input::print("init_data done");
}

void Data::load_data_(
    Input &input,
    const int ndata,
    const int fftdim,
    std::string *dir 
)
{
    // Input::print("load_data");
    if (ndata <= 0){ 
        return;
    }

    std::vector<long unsigned int> cshape = {(long unsigned) nx};
    std::vector<double> container(nx);
    bool fortran_order = false;

    for (int idata = 0; idata < ndata; ++idata)
    {
        this->loadTensor(dir[idata] + "/rho.npy", cshape, fortran_order, container, idata, fftdim, rho);
        if (this->load_gamma){
            this->loadTensor(dir[idata] + "/gamma.npy", cshape, fortran_order, container, idata, fftdim, gamma);
        }
        if (this->load_p){
            this->loadTensor(dir[idata] + "/p.npy", cshape, fortran_order, container, idata, fftdim, p);
            npy::LoadArrayFromNumpy(dir[idata] + "/nablaRhox.npy", cshape, fortran_order, container);
            nablaRho[idata][0] = torch::tensor(container).reshape({fftdim, fftdim, fftdim});
            npy::LoadArrayFromNumpy(dir[idata] + "/nablaRhoy.npy", cshape, fortran_order, container);
            nablaRho[idata][1] = torch::tensor(container).reshape({fftdim, fftdim, fftdim});
            npy::LoadArrayFromNumpy(dir[idata] + "/nablaRhoz.npy", cshape, fortran_order, container);
            nablaRho[idata][2] = torch::tensor(container).reshape({fftdim, fftdim, fftdim});
        }
        if (this->load_q){
            this->loadTensor(dir[idata] + "/q.npy", cshape, fortran_order, container, idata, fftdim, q);
        }
        if (this->load_tanhp){
            this->loadTensor(dir[idata] + "/tanhp.npy", cshape, fortran_order, container, idata, fftdim, tanhp);
        }
        if (this->load_tanhq){
            this->loadTensor(dir[idata] + "/tanhq.npy", cshape, fortran_order, container, idata, fftdim, tanhq);
        }

        for (int ik = 0; ik < input.nkernel; ++ik)
        {
            int ktype = input.kernel_type[ik];
            double kscaling = input.kernel_scaling[ik];

            if (this->load_gammanl[ik]){
                this->loadTensor(dir[idata] + this->file_name("gammanl", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, gammanl[ik]);
            }
            if (this->load_pnl[ik]){
                this->loadTensor(dir[idata] + this->file_name("pnl", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, pnl[ik]);
            }
            if (this->load_qnl[ik]){
                this->loadTensor(dir[idata] + this->file_name("qnl", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, qnl[ik]);
            }
            if (this->load_xi[ik]){
                this->loadTensor(dir[idata] + this->file_name("xi", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, xi[ik]);
            }
            if (this->load_tanhxi[ik]){
                this->loadTensor(dir[idata] + this->file_name("tanhxi", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, tanhxi[ik]);
            }
            if (this->load_tanhxi_nl[ik]){
                this->loadTensor(dir[idata] + this->file_name("tanhxi_nl", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, tanhxi_nl[ik]);
            }
            if (this->load_tanh_pnl[ik]){
                this->loadTensor(dir[idata] + this->file_name("tanh_pnl", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, tanh_pnl[ik]);
            }
            if (this->load_tanh_qnl[ik]){
                this->loadTensor(dir[idata] + this->file_name("tanh_qnl", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, tanh_qnl[ik]);
            }
            if (this->load_tanhp_nl[ik]){
                this->loadTensor(dir[idata] + this->file_name("tanhp_nl", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, tanhp_nl[ik]);
            }
            if (this->load_tanhq_nl[ik]){
                this->loadTensor(dir[idata] + this->file_name("tanhq_nl", ktype, kscaling), cshape, fortran_order, container, idata, fftdim, tanhq_nl[ik]);
            }
        }

        this->loadTensor(dir[idata] + "/enhancement.npy", cshape, fortran_order, container, idata, fftdim, enhancement);
        enhancement_mean[idata] = torch::mean(enhancement[idata]);
        tau_mean[idata] = torch::mean(torch::pow(rho[idata], input.exponent/3.) * enhancement[idata]);

        if (input.loss == "potential" || input.loss == "both" || input.loss == "both_new")
        {
            this->loadTensor(dir[idata] + "/pauli.npy", cshape, fortran_order, container, idata, fftdim, pauli);
            pauli_mean[idata] = torch::mean(pauli[idata]);
        }
    }
    enhancement.resize_({this->nx_tot, 1});
    pauli.resize_({nx_tot, 1});

    this->tau_tf = this->cTF * torch::pow(this->rho, 5./3.);
    // Input::print("load_data done");
}

void Data::loadTensor(
    std::string file,
    std::vector<long unsigned int> cshape,
    bool fortran_order, 
    std::vector<double> &container,
    const int index,
    const int fftdim,
    torch::Tensor &data
)
{
    npy::LoadArrayFromNumpy(file, cshape, fortran_order, container);
    data[index] = torch::tensor(container).reshape({fftdim, fftdim, fftdim});
}

void Data::dumpTensor(const torch::Tensor &data, std::string filename, int nx)
{
    std::vector<double> v(nx);
    for (int ir = 0; ir < nx; ++ir){
        v[ir] = data[ir].item<double>();
    }
    // std::vector<double> v(data.data_ptr<float>(), data.data_ptr<float>() + data.numel()); // this works, but only supports float tensor
    const long unsigned cshape[] = {(long unsigned) nx}; // shape
    npy::SaveArrayAsNumpy(filename, false, 1, cshape, v);
    std::cout << "Dumping " << filename << " done" << std::endl;
}

std::string Data::file_name(std::string parameter, const int kernel_type, const double kernel_scaling)
{
    std::stringstream ss;
    ss << "/" << parameter << "_" << kernel_type << "_" << kernel_scaling << ".npy";
    return ss.str();
}
