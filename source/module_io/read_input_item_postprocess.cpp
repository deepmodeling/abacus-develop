#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_postprocess()
{
    // DOS
    {
        Input_Item item("dos_emin_ev");
        item.annotation = "minimal range for dos";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.dos_emin_ev = doublevalue;
            para.sys.dos_setemin = true;
        };
        sync_double(input.dos_emin_ev);
        add_bool_bcast(sys.dos_setemin);
        this->add_item(item);
    }
    {
        Input_Item item("dos_emax_ev");
        item.annotation = "maximal range for dos";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.dos_emax_ev = doublevalue;
            para.sys.dos_setemax = true;
        };
        sync_double(input.dos_emax_ev);
        add_bool_bcast(sys.dos_setemax);
        this->add_item(item);
    }
    {
        Input_Item item("dos_edelta_ev");
        item.annotation = "delta energy for dos";
        read_sync_double(input.dos_edelta_ev);
        this->add_item(item);
    }
    {
        Input_Item item("dos_scale");
        item.annotation = "scale dos range by";
        read_sync_double(input.dos_scale);
        this->add_item(item);
    }
    {
        Input_Item item("dos_sigma");
        item.annotation = "gauss b coefficeinet(default=0.07)";
        read_sync_double(input.dos_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("dos_nche");
        item.annotation = "orders of Chebyshev expansions for dos";
        read_sync_int(input.dos_nche);
        this->add_item(item);
    }

    // Electronic Conductivity
    {
        Input_Item item("cal_cond");
        item.annotation = "calculate electronic conductivities";
        read_sync_bool(input.cal_cond);
        this->add_item(item);
    }
    {
        Input_Item item("cond_che_thr");
        item.annotation = "control the error of Chebyshev expansions for conductivities";
        read_sync_double(input.cond_che_thr);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dw");
        item.annotation = "frequency interval for conductivities";
        read_sync_double(input.cond_dw);
        this->add_item(item);
    }
    {
        Input_Item item("cond_wcut");
        item.annotation = "cutoff frequency (omega) for conductivities";
        read_sync_double(input.cond_wcut);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dt");
        item.annotation = "t interval to integrate Onsager coefficiencies";
        read_sync_double(input.cond_dt);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dtbatch");
        item.annotation = "exp(iH*dt*cond_dtbatch) is expanded with Chebyshev expansion";
        read_sync_int(input.cond_dtbatch);
        this->add_item(item);
    }
    {
        Input_Item item("cond_smear");
        item.annotation = "Smearing method for conductivities";
        read_sync_int(input.cond_smear);
        this->add_item(item);
    }
    {
        Input_Item item("cond_fwhm");
        item.annotation = "FWHM for conductivities";
        read_sync_double(input.cond_fwhm);
        this->add_item(item);
    }
    {
        Input_Item item("cond_nonlocal");
        item.annotation = "Nonlocal effects for conductivities";
        read_sync_bool(input.cond_nonlocal);
        this->add_item(item);
    }

    // berry_wannier
    {
        Input_Item item("berry_phase");
        item.annotation = "calculate berry phase or not";
        read_sync_bool(input.berry_phase);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.berry_phase)
            {
                if (para.input.basis_type != "pw" && para.input.basis_type != "lcao")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "calculate berry phase, please "
                                             "set basis_type = pw or lcao");
                }
                if (para.input.calculation != "nscf")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate berry phase, please set calculation = nscf");
                }
                if (!(para.input.gdir == 1 || para.input.gdir == 2 || para.input.gdir == 3))
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate berry phase, please set gdir = 1 or 2 or 3");
                }
                if (para.input.symmetry != "-1")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate berry phase, please set symmetry = -1");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("gdir");
        item.annotation = "calculate the polarization in the direction of the "
                          "lattice vector";
        read_sync_int(input.gdir);
        this->add_item(item);
    }
    {
        Input_Item item("towannier90");
        item.annotation = "use wannier90 code interface or not";
        read_sync_bool(input.towannier90);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.towannier90)
            {
                if (para.input.calculation != "nscf")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "to use towannier90, please set calculation = nscf");
                }
                if (para.input.nspin == 2)
                {
                    if (para.input.wannier_spin != "up" && para.input.wannier_spin != "down")
                    {
                        ModuleBase::WARNING_QUIT("ReadInput",
                                                 "to use towannier90, please set wannier_spin = up "
                                                 "or down");
                    }
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nnkpfile");
        item.annotation = "the wannier90 code nnkp file name";
        read_sync_string(input.nnkpfile);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_spin");
        item.annotation = "calculate spin in wannier90 code interface";
        read_sync_string(input.wannier_spin);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_method");
        item.annotation = "different implementation methods under Lcao basis set";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            /*
                       Developer's notes: on the repair of lcao_in_pw

                       lcao_in_pw is a special basis_type, for scf calculation,
                      it follows workflow of pw, but for nscf the toWannier90
                      calculation, the interface is in ESolver_KS_LCAO_elec,
                       therefore lcao_in_pw for towannier90 calculation follows
                      lcao.

                       In the future lcao_in_pw will have its own ESolver.

                       2023/12/22 use new psi_initializer to expand numerical
                      atomic orbitals, ykhuang
                   */
            if (para.input.towannier90 && para.input.basis_type == "lcao_in_pw")
            {
                para.input.wannier_method = 1;
            }
        };
        read_sync_int(input.wannier_method);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_mmn");
        item.annotation = "output .mmn file or not";
        read_sync_bool(input.out_wannier_mmn);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_amn");
        item.annotation = "output .amn file or not";
        read_sync_bool(input.out_wannier_amn);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_unk");
        item.annotation = "output UNK. file or not";
        read_sync_bool(input.out_wannier_unk);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_eig");
        item.annotation = "output .eig file or not";
        read_sync_bool(input.out_wannier_eig);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_wvfn_formatted");
        item.annotation = "output UNK. file in text format or in binary format";
        read_sync_bool(input.out_wannier_wvfn_formatted);
        this->add_item(item);
    }
}
} // namespace ModuleIO