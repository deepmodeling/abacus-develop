#include "./train_kedf.h"
#include <sstream>
#include <math.h>
#include <chrono>

Train_KEDF::~Train_KEDF()
{
    delete[] this->train_volume;
    delete[] this->vali_volume;
    delete[] this->kernel_train;
    delete[] this->kernel_vali;
}

void Train_KEDF::setUpFFT()
{
    this->train_volume = new double[this->input.ntrain];
    this->grid_train.initGrid(
        this->input.fftdim,
        this->input.ntrain,
        this->input.train_cell,
        this->input.train_a,
        this->device,
        this->train_volume
    );
    this->kernel_train = new Kernel[this->input.nkernel];
    for (int ik = 0; ik < this->input.nkernel; ++ik)
    {
        this->kernel_train[ik].set_para(this->input.kernel_type[ik], this->input.kernel_scaling[ik], this->input.yukawa_alpha[ik], this->input.kernel_file[ik]);
        this->kernel_train[ik].fill_kernel(
            this->input.fftdim,
            this->input.ntrain,
            this->data_train.rho,
            this->train_volume,
            this->input.train_cell,
            this->device,
            this->grid_train.fft_gg
        );
    }
    if (this->input.nvalidation > 0){
        this->vali_volume = new double[this->input.nvalidation];
        this->grid_vali.initGrid(
            this->input.fftdim,
            this->input.nvalidation,
            this->input.validation_cell,
            this->input.validation_a,
            this->device,
            this->vali_volume
        );
        this->kernel_vali = new Kernel[this->input.nkernel];
        for (int ik = 0; ik < this->input.nkernel; ++ik)
        {
            this->kernel_vali[ik].set_para(this->input.kernel_type[ik], this->input.kernel_scaling[ik], this->input.yukawa_alpha[ik], this->input.kernel_file[ik]);
            this->kernel_vali[ik].fill_kernel(
                this->input.fftdim,
                this->input.nvalidation,
                this->data_vali.rho,
                this->vali_volume,
                this->input.validation_cell,
                this->device,
                this->grid_vali.fft_gg
            );
        }
    }
    
    // this->dumpTensor(this->fft_kernel_train[0].reshape({this->data_train.nx}), "kernel_fcc.npy", this->data_train.nx);
    // this->dumpTensor(this->fft_kernel_vali[0].reshape({this->data_train.nx}), "kernel_bcc.npy", this->data_train.nx);
}

void Train_KEDF::set_device()
{
    if (this->input.device_type == "cpu")
    {
        std::cout << "---------- Running on CPU ----------" << std::endl;
        this->device = torch::Device(torch::kCPU);
    }
    else if (this->input.device_type == "gpu")
    {
        if (torch::cuda::cudnn_is_available())
        {
            std::cout << "---------- Running on GPU ----------" << std::endl;
            this->device = torch::Device(torch::kCUDA);
        }
        else
        {
            std::cout << "---- Warning: GPU is unaviable -----" << std::endl;
            std::cout << "---------- Running on CPU ----------" << std::endl;
            this->device = torch::Device(torch::kCPU);
        }
    }
}

void Train_KEDF::init_input_index()
{
    this->ninput = 0;

    // --------- semi-local descriptors ---------
    if (this->input.ml_gamma){
        this->descriptor_type.push_back("gamma");
        this->kernel_index.push_back(-1);
        ninput++;
    } 
    if (this->input.ml_p){
        this->descriptor_type.push_back("p");
        this->kernel_index.push_back(-1);
        ninput++;
    }
    if (this->input.ml_q){
        this->descriptor_type.push_back("q");
        this->kernel_index.push_back(-1);
        ninput++;
    }
    // --------- non-local descriptors ---------
    for (int ik = 0; ik < this->input.nkernel; ++ik)
    {
        if (this->input.ml_gammanl[ik]){
            this->descriptor_type.push_back("gammanl");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
        if (this->input.ml_pnl[ik]){
            this->descriptor_type.push_back("pnl");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
        if (this->input.ml_qnl[ik]){
            this->descriptor_type.push_back("qnl");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
        if (this->input.ml_xi[ik]){
            this->descriptor_type.push_back("xi");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
        if (this->input.ml_tanhxi[ik]){
            this->descriptor_type.push_back("tanhxi");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
        if (this->input.ml_tanhxi_nl[ik]){
            this->descriptor_type.push_back("tanhxi_nl");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
    }
    // --------- semi-local descriptors ---------
    if (this->input.ml_tanhp){
        this->descriptor_type.push_back("tanhp");
        this->kernel_index.push_back(-1);
        ninput++;
    }
    if (this->input.ml_tanhq){
        this->descriptor_type.push_back("tanhq");
        this->kernel_index.push_back(-1);
        ninput++;
    }
    // --------- non-local descriptors ---------
    for (int ik = 0; ik < this->input.nkernel; ++ik)
    {
        if (this->input.ml_tanh_pnl[ik]){
            this->descriptor_type.push_back("tanh_pnl");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
        if (this->input.ml_tanh_qnl[ik]){
            this->descriptor_type.push_back("tanh_qnl");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
        if (this->input.ml_tanhp_nl[ik]){
            this->descriptor_type.push_back("tanhp_nl");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
        if (this->input.ml_tanhq_nl[ik]){
            this->descriptor_type.push_back("tanhq_nl");
            this->kernel_index.push_back(ik);
            this->ninput++;
        }
    }

    std::cout << "ninput = " << ninput << std::endl;

    if (this->input.feg_limit != 0)
    {
        this->feg_inpt = torch::zeros(ninput).to(device);
        for (int i = 0; i < this->ninput; ++i)
        {
            if (this->descriptor_type[i] == "gamma") this->feg_inpt[i] = 1.;
        }

        this->feg_predict = torch::zeros(1).to(device);
        this->feg_dFdgamma = torch::zeros(1).to(device);
    }
    std::cout << "feg_limit = " << this->input.feg_limit << std::endl;
}

void Train_KEDF::init()
{
    this->set_device();
    this->init_input_index();
    this->data_train.load_data(this->input, this->input.ntrain, this->input.train_dir, this->device);
    this->data_vali.load_data(this->input, this->input.nvalidation, this->input.validation_dir, this->device);
    // Input::print("LOAD DATA done");

    this->potential.init(this->input, this->ninput, this->descriptor_type, this->kernel_index);
    // Input::print("init potential done");
    if (this->input.loss == "potential" || this->input.loss == "both" || this->input.loss == "both_new") this->setUpFFT();
    // Input::print("init fft done");
    
    this->nn = std::make_shared<NN_OFImpl>(this->data_train.nx_tot, this->data_vali.nx_tot, this->ninput, this->input.nnode, this->input.nlayer, this->device);
    // Input::print("init_nn done");
    this->nn->set_data(&(this->data_train), this->descriptor_type, this->kernel_index, this->nn->inputs);
    this->nn->set_data(&(this->data_vali), this->descriptor_type, this->kernel_index, this->nn->input_vali);
}

torch::Tensor Train_KEDF::lossFunction(torch::Tensor enhancement, torch::Tensor target, torch::Tensor coef)
{
    return torch::sum(torch::pow(enhancement - target, 2))/this->data_train.nx/coef/coef;
}

torch::Tensor Train_KEDF::lossFunction_new(torch::Tensor enhancement, torch::Tensor target, torch::Tensor weight, torch::Tensor coef)
{
    return torch::sum(torch::pow(weight * (enhancement - target), 2.))/this->data_train.nx/coef/coef;
}


void Train_KEDF::train()
{
    // time
    double tot = 0.;
    double totF = 0.;
    double totFback = 0.;
    double totLoss = 0.;
    double totLback = 0.;
    double totP = 0.;
    double totStep = 0.;
    std::chrono::_V2::system_clock::time_point start, startF, startFB, startL, startLB, startP, startStep, end, endF, endFB, endL, endLB, endP, endStep;

    start = std::chrono::high_resolution_clock::now();

    std::cout << "========== Train_KEDF begin ==========" << std::endl;
    // torch::Tensor target = (this->input.loss=="energy") ? this->data_train.enhancement : this->data_train.pauli;
    if (this->input.loss == "potential" || this->input.loss == "both" || this->input.loss == "both_new")
    {
        this->data_train.pauli.resize_({this->input.ntrain, this->input.fftdim, this->input.fftdim, this->input.fftdim});
    }

    torch::optim::SGD optimizer(this->nn->parameters(), this->input.lr_start);
    double update_coef = this->input.nepoch/std::log(this->input.lr_end/this->input.lr_start); // used to reduce the learning rate

    std::cout << "Epoch\tLoss\tValidation\tLoss_pot\tLoss_E\tLoss_FEG_pot\tLoss_FEG_E\n";
    double lossTrain = 0.;
    double lossPot = 0.;
    double lossE = 0.;
    double lossFEG_pot = 0.;
    double lossFEG_E = 0.;
    double lossVali = 0.;
    double maxLoss = 100.;

    // bool increase_coef_feg_e = false;
    torch::Tensor weight = torch::pow(this->data_train.rho, this->input.exponent/3.);
    for (size_t epoch = 1; epoch <= this->input.nepoch; ++epoch)
    {
        for (int batch_index = 0; batch_index < this->input.ntrain; ++batch_index)
        {
            startF = std::chrono::high_resolution_clock::now();

            optimizer.zero_grad();
            if (this->input.loss == "energy")
            {
                torch::Tensor inpt = torch::slice(this->nn->inputs, 0, batch_index*this->data_train.nx, (batch_index + 1)*this->data_train.nx);
                inpt.requires_grad_(true);
                torch::Tensor prediction = this->nn->forward(inpt);
                startL = std::chrono::high_resolution_clock::now();
                torch::Tensor loss = this->lossFunction(prediction, torch::slice(this->data_train.enhancement, 0, batch_index*this->data_train.nx, (batch_index + 1)*this->data_train.nx), this->data_train.enhancement_mean[batch_index]) * this->input.coef_e;
                lossTrain = loss.item<double>();
                endL = std::chrono::high_resolution_clock::now();
                totLoss += double(std::chrono::duration_cast<std::chrono::microseconds>(endL - startL).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

                startLB = std::chrono::high_resolution_clock::now();
                loss.backward();
                endLB = std::chrono::high_resolution_clock::now();
                totLback += double(std::chrono::duration_cast<std::chrono::microseconds>(endLB - startLB).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
            }
            else if (this->input.loss == "potential" || this->input.loss == "both" || this->input.loss == "both_new")
            {
                torch::Tensor inpt = torch::slice(this->nn->inputs, 0, batch_index*this->data_train.nx, (batch_index + 1)*this->data_train.nx);
                inpt.requires_grad_(true);
                torch::Tensor prediction = this->nn->forward(inpt);

                if (this->input.feg_limit != 3)
                {
                    prediction = torch::softplus(prediction);
                }
                if (this->input.feg_limit != 0)
                {
                    // if (this->feg_inpt.grad().numel()) this->feg_inpt.grad().zero_();
                    this->feg_predict = this->nn->forward(this->feg_inpt);
                    // if (this->input.ml_gamma) this->feg_dFdgamma = torch::autograd::grad({this->feg_predict}, {this->feg_inpt},
                    //                                                                 {torch::ones_like(this->feg_predict)}, true, true)[0][this->nn_input_index["gamma"]];
                    if (this->input.feg_limit == 1) prediction = prediction - torch::softplus(this->feg_predict) + 1.;
                    if (this->input.feg_limit == 3 && epoch < this->input.change_step) prediction = torch::softplus(prediction);
                    if (this->input.feg_limit == 3 && epoch >= this->input.change_step)  prediction = torch::softplus(prediction - this->feg_predict + this->feg3_correct);
                }
                endF = std::chrono::high_resolution_clock::now();
                totF += double(std::chrono::duration_cast<std::chrono::microseconds>(endF - startF).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

                startFB = std::chrono::high_resolution_clock::now();
                torch::Tensor gradient = torch::autograd::grad({prediction}, {inpt},
                                                           {torch::ones_like(prediction)}, true, true)[0];
                endFB = std::chrono::high_resolution_clock::now();
                totFback += double(std::chrono::duration_cast<std::chrono::microseconds>(endFB - startFB).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

                startL = std::chrono::high_resolution_clock::now();
                torch::Tensor pot = this->potential.get_potential(batch_index, this->data_train, prediction.reshape({this->input.fftdim, this->input.fftdim, this->input.fftdim}), gradient, this->kernel_train, this->grid_train);
                torch::Tensor loss = this->lossFunction(pot, this->data_train.pauli[batch_index], this->data_train.pauli_mean[batch_index])
                                     * this->input.coef_p;
                lossPot = loss.item<double>();
                if (this->input.loss == "both")
                {
                    loss = loss + this->input.coef_e * this->lossFunction(prediction, torch::slice(this->data_train.enhancement, 0, batch_index*this->data_train.nx, (batch_index + 1)*this->data_train.nx), this->data_train.enhancement_mean[batch_index]);
                    lossE = loss.item<double>() - lossPot;
                }
                if (this->input.loss == "both_new")
                {
                    loss = loss + this->input.coef_e * this->lossFunction_new(prediction, torch::slice(this->data_train.enhancement, 0, batch_index*this->data_train.nx, (batch_index + 1)*this->data_train.nx), weight[batch_index].reshape({this->data_train.nx, 1}), this->data_train.tau_mean[batch_index]);
                    lossE = loss.item<double>() - lossPot;
                }
                if (this->input.feg_limit != 0)
                {
                    if (this->input.feg_limit == 1 || this->input.feg_limit == 2)
                    {
                        loss = loss + torch::pow(this->feg_predict - 1., 2) * this->input.coef_feg_e;
                        lossFEG_E = loss.item<double>() - (lossPot + lossE + lossFEG_pot);
                        // if (lossFEG_E/lossE < 1e-3 && increase_coef_feg_e == false)
                        // {
                        //     this->input.coef_feg_e *= 2.;
                        //     increase_coef_feg_e = true;
                        //     std::cout << "---------ICREASE COEF FEG E--------" << std::endl;
                        // }
                    }
                    if (this->input.feg_limit == 3)
                    {
                        loss = loss + torch::pow(this->feg_predict - this->feg3_correct, 2) * this->input.coef_feg_e;
                        lossFEG_E = loss.item<double>() - (lossPot + lossE + lossFEG_pot);
                    }
                }

                lossTrain = loss.item<double>();
                endL = std::chrono::high_resolution_clock::now();
                totLoss += double(std::chrono::duration_cast<std::chrono::microseconds>(endL - startL).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
                
                startLB = std::chrono::high_resolution_clock::now();
                loss.backward();
                endLB = std::chrono::high_resolution_clock::now();
                totLback += double(std::chrono::duration_cast<std::chrono::microseconds>(endLB - startLB).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
                // this->dumpTensor(pot.reshape({this->data_train.nx}), "pot_fcc.npy", this->data_train.nx);
                // this->dumpTensor(torch::slice(prediction, 0, batch_index*this->data_train.nx, (batch_index + 1)*this->data_train.nx).reshape({this->data_train.nx}), "F_fcc.npy", this->data_train.nx);
            }
            
            startP = std::chrono::high_resolution_clock::now();
            if (epoch % this->input.print_fre == 0) {
                if (this->input.nvalidation > 0)
                {
                    torch::Tensor valid_pre = this->nn->forward(this->nn->input_vali);
                    if (this->input.feg_limit == 3)
                    {
                        valid_pre = torch::softplus(valid_pre - this->feg_predict + this->feg3_correct);
                    }
                    lossVali = this->lossFunction(valid_pre, this->data_vali.enhancement, this->data_vali.enhancement_mean).item<double>();
                }
                std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(3) << epoch 
                          << std::setw(12) << lossTrain 
                          << std::setw(12) << lossVali 
                          << std::setw(12) << lossPot
                          << std::setw(12) << lossE 
                          << std::setw(12) << lossFEG_pot
                          << std::setw(12) << lossFEG_E
                          << std::endl;
            }
            
            if (lossTrain > maxLoss){
                std::cout << "ERROR: too large loss: " << lossTrain << std::endl;
                exit(0);
            }
            endP = std::chrono::high_resolution_clock::now();
            totP += double(std::chrono::duration_cast<std::chrono::microseconds>(endP - startP).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

            startStep = std::chrono::high_resolution_clock::now();
            optimizer.step();
            endStep = std::chrono::high_resolution_clock::now();
            totStep += double(std::chrono::duration_cast<std::chrono::microseconds>(endStep - startStep).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
        }
        if (epoch % this->input.dump_fre == 0)
        {
            std::stringstream file;
            file << "model/net" << epoch << ".pt";
            torch::save(this->nn, file.str());
        }
        // Reduce the learning_rate
        if (epoch % this->input.lr_fre == 0)
        {
            for (auto &group : optimizer.param_groups())
            {
                if(group.has_options())
                {
                    auto &options = static_cast<torch::optim::SGDOptions &>(group.options());
                    options.lr(this->input.lr_start * std::exp(epoch/update_coef));
                }
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();

    tot = double(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

    std::cout << "=========== Done ============" << std::endl;
    std::cout << std::setprecision(2) << std::setiosflags(std::ios::fixed) << "Item\t\t\tTime (s)\tPercentage" << std::endl;
    std::cout.unsetf(std::ios::scientific);
    std::cout << "Total\t\t\t"          << tot      << "\t\t" << tot/tot        * 100. << " %" << std::endl;
    std::cout << "Forward\t\t\t"        << totF     << "\t\t" << totF/tot       * 100. << " %" << std::endl;
    std::cout << "Enhancement back\t"   << totFback << "\t\t" << totFback/tot   * 100. << " %" << std::endl;
    std::cout << "Loss function\t\t"    << totLoss  << "\t\t" << totLoss/tot    * 100. << " %" << std::endl;
    std::cout << "Loss backward\t\t"    << totLback << "\t\t" << totLback/tot   * 100. << " %" << std::endl;
    std::cout << "Print\t\t\t"          << totP     << "\t\t" << totP/tot       * 100. << " %" << std::endl;
    std::cout << "Step\t\t\t"           << totStep  << "\t\t" << totStep/tot    * 100. << " %" << std::endl;
}

void Train_KEDF::potTest()
{
    this->set_device();
    this->init_input_index();
    this->data_train.load_data(this->input, this->input.ntrain, this->input.train_dir, this->device);
    this->data_vali.load_data(this->input, this->input.nvalidation, this->input.validation_dir, this->device);
    Input::print("LOAD DATA done");

    this->potential.init(this->input, this->ninput, this->descriptor_type, this->kernel_index);
    Input::print("init potential done");
    if (this->input.loss == "potential" || this->input.loss == "both" || this->input.loss == "both_new") this->setUpFFT();
    Input::print("init fft done");
    
    std::chrono::_V2::system_clock::time_point start, end;
    this->nn = std::make_shared<NN_OFImpl>(this->data_train.nx_tot, this->data_vali.nx_tot, this->ninput, this->input.nnode, this->input.nlayer, this->device);
    torch::DeviceType device_type = torch::kCPU;
    torch::load(this->nn, "net.pt", device_type);
    // this->nn->to(this->device);
    Input::print("init_nn done");
    this->nn->set_data(&(this->data_train), this->descriptor_type, this->kernel_index, this->nn->inputs);
    this->nn->set_data(&(this->data_vali), this->descriptor_type, this->kernel_index, this->nn->input_vali);

    this->nn->inputs.requires_grad_(true);

    if (this->input.loss == "potential" || this->input.loss == "both" || this->input.loss == "both_new") this->data_train.pauli.resize_({this->input.ntrain, this->input.fftdim, this->input.fftdim, this->input.fftdim});

    for (int batch_index = 0; batch_index < this->input.ntrain; ++batch_index)
    {
        for (int ii = 0; ii < 1; ++ii)
        {
            torch::Tensor inpts = torch::slice(this->nn->inputs, 0, ii*this->data_train.nx, (ii + 1)*this->data_train.nx);
            inpts.requires_grad_(true);
            torch::Tensor prediction = this->nn->forward(inpts);
            if (this->input.feg_limit != 3)
            {
                prediction = torch::softplus(prediction);
            }
            if (this->input.feg_limit != 0)
            {
                // if (this->input.ml_gamma) if (this->feg_inpt[this->nn_input_index["gamma"]].grad().numel()) this->feg_inpt[this->nn_input_index["gamma"]].grad().zero_();
                if (this->feg_inpt.grad().numel()) this->feg_inpt.grad().zero_();
                this->feg_predict = this->nn->forward(this->feg_inpt);
                // if (this->input.ml_gamma) this->feg_dFdgamma = torch::autograd::grad({this->feg_predict}, {this->feg_inpt[this->nn_input_index["gamma"]]},
                //                                                                 {torch::ones(1)}, true, true)[0];
                // if (this->input.ml_gamma) this->feg_dFdgamma = torch::autograd::grad({this->feg_predict}, {this->feg_inpt},
                //                                                                 {torch::ones_like(this->feg_predict)}, true, true)[0][this->nn_input_index["gamma"]];
                if (this->input.feg_limit == 1) prediction = prediction - torch::softplus(this->feg_predict) + 1.;
                if (this->input.feg_limit == 3) prediction = torch::softplus(prediction - this->feg_predict + this->feg3_correct);
            }

            start = std::chrono::high_resolution_clock::now();
            torch::Tensor gradient = torch::autograd::grad({prediction}, {inpts},
                                                    {torch::ones_like(prediction)}, true, true)[0];
            end = std::chrono::high_resolution_clock::now();
            std::cout << "spend " << double(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << " s" << std::endl;
            
            std::cout << "begin potential" << std::endl;
            torch::Tensor weight = torch::pow(this->data_train.rho, this->input.exponent/3.);

            torch::Tensor pot = this->potential.get_potential(ii, this->data_train, torch::slice(prediction, 0, ii*this->data_train.nx, (ii + 1)*this->data_train.nx).reshape({this->input.fftdim, this->input.fftdim, this->input.fftdim}), gradient, this->kernel_train, this->grid_train);
            std::cout << "after potential" << std::endl;

            torch::Tensor loss = this->lossFunction(pot, this->data_train.pauli[ii], this->data_train.pauli_mean[ii]) * this->input.coef_p;
            if (this->input.loss == "both")
            {
                loss = loss + this->input.coef_e * this->lossFunction(prediction, torch::slice(this->data_train.enhancement, 0, ii*this->data_train.nx, (ii + 1)*this->data_train.nx), this->data_train.enhancement_mean[ii]);
            }
            if (this->input.loss == "both_new")
            {
                loss = loss + this->input.coef_e * this->lossFunction_new(prediction, torch::slice(this->data_train.enhancement, 0, ii*this->data_train.nx, (ii + 1)*this->data_train.nx), weight[ii].reshape({this->data_train.nx, 1}), this->data_train.tau_mean[ii]);
            }
            std::cout << "after loss" << std::endl;
            // loss = loss + this->input.coef_e * this->lossFunction(prediction, torch::slice(this->data_train.enhancement, 0, ii*this->data_train.nx, (ii + 1)*this->data_train.nx));
            double lossTrain = loss.item<double>();
            std::cout << "loss = " << lossTrain << std::endl;
            this->data_train.dumpTensor(pot.reshape({this->data_train.nx}), "potential-nnof.npy", this->data_train.nx);
            this->data_train.dumpTensor(torch::slice(prediction, 0, ii*this->data_train.nx, (ii + 1)*this->data_train.nx).reshape({this->data_train.nx}), "enhancement-nnof.npy", this->data_train.nx);
            // this->dumpTensor(torch::slice(this->nn->F, 0, ii*this->data_train.nx, (ii + 1)*this->data_train.nx).reshape({this->data_train.nx}), "F_fcc.npy", this->data_train.nx);
        }
        // std::cout << std::setiosflags(std::ios::scientific) <<  std::setprecision(12) << this->nn->parameters() << std::endl;
    }
}

