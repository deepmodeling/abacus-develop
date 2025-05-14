#include "./input.h"

void Input::readInput()
{
    std::ifstream ifs("nnINPUT", std::ios::in);
    if (!ifs)
    {
        std::cout << " Can't find the nnINPUT file." << std::endl;
        exit(0);
    }

    char word[80];
    int ierr = 0;

    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word;
        if (ifs.eof())
            break;

        if (strcmp("fftdim", word) == 0)
        {
            this->read_value(ifs, this->fftdim);
        }
        else if (strcmp("nbatch", word) == 0)
        {
            this->read_value(ifs, this->nbatch);
        }
        else if (strcmp("ntrain", word) == 0)
        {
            this->read_value(ifs, this->ntrain);
            this->train_dir = new std::string[this->ntrain];
            this->train_cell = new std::string[this->ntrain];
            this->train_a = new double[this->ntrain];
        }
        else if (strcmp("nvalidation", word) == 0)
        {
            this->read_value(ifs, this->nvalidation);
            if (this->nvalidation > 0)
            {
                this->validation_dir = new std::string[this->nvalidation];
                this->validation_cell = new std::string[this->nvalidation];
                this->validation_a = new double[this->nvalidation];
            }
        }
        else if (strcmp("train_dir", word) == 0)
        {
            this->read_values(ifs, this->ntrain, this->train_dir);
        }
        else if (strcmp("train_cell", word) == 0)
        {
            this->read_values(ifs, this->ntrain, this->train_cell);
        }
        else if (strcmp("train_a", word) == 0)
        {
            this->read_values(ifs, this->ntrain, this->train_a);
        }
        else if (strcmp("validation_dir", word) == 0)
        {
            this->read_values(ifs, this->nvalidation, this->validation_dir);
        }
        else if (strcmp("validation_cell", word) == 0 && this->nvalidation > 0)
        {
            this->read_values(ifs, this->nvalidation, this->validation_cell);
        }
        else if (strcmp("validation_a", word) == 0 && this->nvalidation > 0)
        {
            this->read_values(ifs, this->nvalidation, this->validation_a);
        }
        else if (strcmp("loss", word) == 0)
        {
            this->read_value(ifs, this->loss);
        }
        else if (strcmp("exponent", word) == 0)
        {
            this->read_value(ifs, this->exponent);
        }
        else if (strcmp("nepoch", word) == 0)
        {
            this->read_value(ifs, this->nepoch);
        }
        else if (strcmp("lr_start", word) == 0)
        {
            this->read_value(ifs, this->lr_start);
        }
        else if (strcmp("lr_end", word) == 0)
        {
            this->read_value(ifs, this->lr_end);
        }
        else if (strcmp("lr_fre", word) == 0)
        {
            this->read_value(ifs, this->lr_fre);
        }
        else if (strcmp("dump_fre", word) == 0)
        {
            this->read_value(ifs, this->dump_fre);
        }
        else if (strcmp("print_fre", word) == 0)
        {
            this->read_value(ifs, this->print_fre);
        }
        else if (strcmp("gamma", word) == 0)
        {
            this->read_value(ifs, this->ml_gamma);
        }
        else if (strcmp("p", word) == 0)
        {
            this->read_value(ifs, this->ml_p);
        }
        else if (strcmp("q", word) == 0)
        {
            this->read_value(ifs, this->ml_q);
        }
        else if (strcmp("gammanl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_gammanl);
        }
        else if (strcmp("pnl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_pnl);
        }
        else if (strcmp("qnl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_qnl);
        }
        else if (strcmp("xi", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_xi);
        }
        else if (strcmp("tanhxi", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_tanhxi);
        }
        else if (strcmp("tanhxi_nl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_tanhxi_nl);
        }
        else if (strcmp("tanhp", word) == 0)
        {
            this->read_value(ifs, this->ml_tanhp);
        }
        else if (strcmp("tanhq", word) == 0)
        {
            this->read_value(ifs, this->ml_tanhq);
        }
        else if (strcmp("tanh_pnl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_tanh_pnl);
        }
        else if (strcmp("tanh_qnl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_tanh_qnl);
        }
        else if (strcmp("tanhp_nl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_tanhp_nl);
        }
        else if (strcmp("tanhq_nl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->ml_tanhq_nl);
        }
        else if (strcmp("chi_xi", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->chi_xi);
        }
        else if (strcmp("chi_p", word) == 0)
        {
            this->read_value(ifs, this->chi_p);
        }
        else if (strcmp("chi_q", word) == 0)
        {
            this->read_value(ifs, this->chi_q);
        }
        else if (strcmp("chi_pnl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->chi_pnl);
        }
        else if (strcmp("chi_qnl", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->chi_qnl);
        }
        else if (strcmp("feg_limit", word) == 0)
        {
            this->read_value(ifs, this->feg_limit);
        }
        else if (strcmp("change_step", word) == 0)
        {
            this->read_value(ifs, this->change_step);
        }
        else if (strcmp("coef_e", word) == 0)
        {
            this->read_value(ifs, this->coef_e);
        }
        else if (strcmp("coef_p", word) == 0)
        {
            this->read_value(ifs, this->coef_p);
        }
        else if (strcmp("coef_feg_e", word) == 0)
        {
            this->read_value(ifs, this->coef_feg_e);
        }
        else if (strcmp("coef_feg_p", word) == 0)
        {
            this->read_value(ifs, this->coef_feg_p);
        }
        else if (strcmp("check_pot", word) == 0)
        {
            this->read_value(ifs, this->check_pot);
        }
        else if (strcmp("nnode", word) == 0)
        {
            this->read_value(ifs, this->nnode);
        }
        else if (strcmp("nlayer", word) == 0)
        {
            this->read_value(ifs, this->nlayer);
        }
        else if (strcmp("nkernel", word) == 0)
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
        else if (strcmp("kernel_type", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->kernel_type);
        }
        else if (strcmp("yukawa_alpha", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->yukawa_alpha);
        }
        else if (strcmp("kernel_scaling", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->kernel_scaling);
        }
        else if (strcmp("kernel_file", word) == 0)
        {
            this->read_values(ifs, this->nkernel, this->kernel_file);
        }
        else if (strcmp("device_type", word) == 0)
        {
            this->read_value(ifs, this->device_type);
        }
    }

    std::cout << "Read nnINPUT done" << std::endl;
}
