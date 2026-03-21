#include "./input.h"

void Input::readInput()
{
    std::ifstream ifs("nnINPUT", std::ios::in);
    if (!ifs)
    {
        std::cout << " Can't find the nnINPUT file." << std::endl;
        exit(0);
    }

    std::string word;
    int ierr = 0;

    ifs.rdstate();
    while (ifs.good()
    {
        ifs >> word;
        if (ifs.eof()
            break;

        if(word=="fftdim")
        {
            this->read_value(ifs, this->fftdim);
        }
        else if(word=="nbatch")
        {
            this->read_value(ifs, this->nbatch);
        }
        else if(word=="ntrain")
        {
            this->read_value(ifs, this->ntrain);
            this->train_dir = new std::string[this->ntrain];
            this->train_cell = new std::string[this->ntrain];
            this->train_a = new double[this->ntrain];
        }
        else if(word=="nvalidation")
        {
            this->read_value(ifs, this->nvalidation);
            if (this->nvalidation > 0)
            {
                this->validation_dir = new std::string[this->nvalidation];
                this->validation_cell = new std::string[this->nvalidation];
                this->validation_a = new double[this->nvalidation];
            }
        }
        else if(word=="train_dir")
        {
            this->read_values(ifs, this->ntrain, this->train_dir);
        }
        else if(word=="train_cell")
        {
            this->read_values(ifs, this->ntrain, this->train_cell);
        }
        else if(word=="train_a")
        {
            this->read_values(ifs, this->ntrain, this->train_a);
        }
        else if(word=="validation_dir")
        {
            this->read_values(ifs, this->nvalidation, this->validation_dir);
        }
        else if(word=="validation_cell" && this->nvalidation > 0)
        {
            this->read_values(ifs, this->nvalidation, this->validation_cell);
        }
        else if(word=="validation_a" && this->nvalidation > 0)
        {
            this->read_values(ifs, this->nvalidation, this->validation_a);
        }
        else if(word=="loss")
        {
            this->read_value(ifs, this->loss);
        }
        else if(word=="exponent")
        {
            this->read_value(ifs, this->exponent);
        }
        else if(word=="nepoch")
        {
            this->read_value(ifs, this->nepoch);
        }
        else if(word=="lr_start")
        {
            this->read_value(ifs, this->lr_start);
        }
        else if(word=="lr_end")
        {
            this->read_value(ifs, this->lr_end);
        }
        else if(word=="lr_fre")
        {
            this->read_value(ifs, this->lr_fre);
        }
        else if(word=="dump_fre")
        {
            this->read_value(ifs, this->dump_fre);
        }
        else if(word=="print_fre")
        {
            this->read_value(ifs, this->print_fre);
        }
        else if(word=="gamma")
        {
            this->read_value(ifs, this->ml_gamma);
        }
        else if(word=="p")
        {
            this->read_value(ifs, this->ml_p);
        }
        else if(word=="q")
        {
            this->read_value(ifs, this->ml_q);
        }
        else if(word=="gammanl")
        {
            this->read_values(ifs, this->nkernel, this->ml_gammanl);
        }
        else if(word=="pnl")
        {
            this->read_values(ifs, this->nkernel, this->ml_pnl);
        }
        else if(word=="qnl")
        {
            this->read_values(ifs, this->nkernel, this->ml_qnl);
        }
        else if(word=="xi")
        {
            this->read_values(ifs, this->nkernel, this->ml_xi);
        }
        else if(word=="tanhxi")
        {
            this->read_values(ifs, this->nkernel, this->ml_tanhxi);
        }
        else if(word=="tanhxi_nl")
        {
            this->read_values(ifs, this->nkernel, this->ml_tanhxi_nl);
        }
        else if(word=="tanhp")
        {
            this->read_value(ifs, this->ml_tanhp);
        }
        else if(word=="tanhq")
        {
            this->read_value(ifs, this->ml_tanhq);
        }
        else if(word=="tanh_pnl")
        {
            this->read_values(ifs, this->nkernel, this->ml_tanh_pnl);
        }
        else if(word=="tanh_qnl")
        {
            this->read_values(ifs, this->nkernel, this->ml_tanh_qnl);
        }
        else if(word=="tanhp_nl")
        {
            this->read_values(ifs, this->nkernel, this->ml_tanhp_nl);
        }
        else if(word=="tanhq_nl")
        {
            this->read_values(ifs, this->nkernel, this->ml_tanhq_nl);
        }
        else if(word=="chi_xi")
        {
            this->read_values(ifs, this->nkernel, this->chi_xi);
        }
        else if(word=="chi_p")
        {
            this->read_value(ifs, this->chi_p);
        }
        else if(word=="chi_q")
        {
            this->read_value(ifs, this->chi_q);
        }
        else if(word=="chi_pnl")
        {
            this->read_values(ifs, this->nkernel, this->chi_pnl);
        }
        else if(word=="chi_qnl")
        {
            this->read_values(ifs, this->nkernel, this->chi_qnl);
        }
        else if(word=="feg_limit")
        {
            this->read_value(ifs, this->feg_limit);
        }
        else if(word=="change_step")
        {
            this->read_value(ifs, this->change_step);
        }
        else if(word=="coef_e")
        {
            this->read_value(ifs, this->coef_e);
        }
        else if(word=="coef_p")
        {
            this->read_value(ifs, this->coef_p);
        }
        else if(word=="coef_feg_e")
        {
            this->read_value(ifs, this->coef_feg_e);
        }
        else if(word=="coef_feg_p")
        {
            this->read_value(ifs, this->coef_feg_p);
        }
        else if(word=="check_pot")
        {
            this->read_value(ifs, this->check_pot);
        }
        else if(word=="nnode")
        {
            this->read_value(ifs, this->nnode);
        }
        else if(word=="nlayer")
        {
            this->read_value(ifs, this->nlayer);
        }
        else if(word=="nkernel")
        {
            this->read_value(ifs, this->nkernel);
            this->ml_gammanl = new bool[this->nkernel];
            this->ml_pnl = new bool[this->nkernel];
            this->ml_qnl = new bool[this->nkernel];
            this->ml_xi = new bool[this->nkernel];
            this->ml_tanhxi = new bool[this->nkernel];
            this->ml_tanhxi_nl = new bool[this->nkernel];
            this->ml_tanh_pnl = new bool[this->nkernel];
            this->ml_tanh_qnl = new bool[this->nkernel];
            this->ml_tanhp_nl = new bool[this->nkernel];
            this->ml_tanhq_nl = new bool[this->nkernel];
            this->chi_xi = new double[this->nkernel];
            this->chi_pnl = new double[this->nkernel];
            this->chi_qnl = new double[this->nkernel];
            this->kernel_type = new int[this->nkernel];
            this->kernel_scaling = new double[this->nkernel];
            this->yukawa_alpha = new double[this->nkernel];
            this->kernel_file = new std::string[this->nkernel];
            for (int ik = 0; ik < this->nkernel; ++ik)
            {
                this->ml_gammanl[ik] = 0;
                this->ml_pnl[ik] = 0;
                this->ml_qnl[ik] = 0;
                this->ml_xi[ik] = 0;
                this->ml_tanhxi[ik] = 0;
                this->ml_tanhxi_nl[ik] = 0;
                this->ml_tanh_pnl[ik] = 0;
                this->ml_tanh_qnl[ik] = 0;
                this->ml_tanhp_nl[ik] = 0;
                this->ml_tanhq_nl[ik] = 0;
                this->chi_xi[ik] = 1.;
                this->chi_pnl[ik] = 1.;
                this->chi_qnl[ik] = 1.;
                this->kernel_type[ik] = 1;
                this->kernel_scaling[ik] = 1.;
                this->yukawa_alpha[ik] = 1.;
                this->kernel_file[ik] = "none";
            }
        }
        else if(word=="kernel_type")
        {
            this->read_values(ifs, this->nkernel, this->kernel_type);
        }
        else if(word=="yukawa_alpha")
        {
            this->read_values(ifs, this->nkernel, this->yukawa_alpha);
        }
        else if(word=="kernel_scaling")
        {
            this->read_values(ifs, this->nkernel, this->kernel_scaling);
        }
        else if(word=="kernel_file")
        {
            this->read_values(ifs, this->nkernel, this->kernel_file);
        }
        else if(word=="device_type")
        {
            this->read_value(ifs, this->device_type);
        }
        else if (strcmp("energy_type", word) == 0)
        {
            this->read_value(ifs, this->energy_type);
        }
    }

    std::cout << "Read nnINPUT done" << std::endl;
}
