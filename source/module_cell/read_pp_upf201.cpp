#include "read_pp.h"

// qianrui rewrite it 2021-5-10
// liuyu update 2023-09-17 add uspp support
int Pseudopot_upf::read_pseudo_upf201(std::ifstream &ifs, Atom_pseudo& pp)
{
    //--------------------------------------
    //-              PP_HEADER             -
    //--------------------------------------
    this->read_pseudo_upf201_header(ifs, pp);

    //--------------------------------------
    //-              PP_MESH               -
    //--------------------------------------
    this->read_pseudo_upf201_mesh(ifs, pp);

    //--------------------------------------
    //-              PP_NLCC               -
    //--------------------------------------
    if (pp.nlcc)
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NLCC", true, false))
        {
            ifs.ignore(150, '>'); // skip type, size, columns and so on.
        }
        else
        {
            ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NLCC>");
        }
        pp.rho_atc = std::vector<double>(pp.mesh, 0.0);
        for (int ir = 0; ir < pp.mesh; ir++)
        {
            ifs >> pp.rho_atc[ir];
        }
        this->skip_number(ifs, this->mesh_changed);
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NLCC>");
    }

    //--------------------------------------
    //-              PP_LOCAL              -
    //--------------------------------------
    if (!this->coulomb_potential)
    {
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_LOCAL", true, false))
        {
            ifs.ignore(150, '>'); // skip type, size, columns and so on.
        }
        else
        {
            ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_LOCAL>");
        }
        pp.vloc_at = std::vector<double>(pp.mesh, 0.0);
        for (int ir = 0; ir < pp.mesh; ir++)
        {
            ifs >> pp.vloc_at[ir];
        }
        this->skip_number(ifs, this->mesh_changed);
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_LOCAL>");
    }

    //--------------------------------------
    //-            PP_NONLOCAL             -
    //--------------------------------------
    this->read_pseudo_upf201_nonlocal(ifs, pp);

    //--------------------------------------
    //-            PP_PSWFC                -
    //--------------------------------------
    this->read_pseudo_upf201_pswfc(ifs, pp);

    //--------------------------------------
    //-            PP_FULL_WFC             -
    //--------------------------------------
    // if (has_wfc)
    // {
    //     this->read_pseudo_upf201_fullwfc(ifs);
    // }

    //--------------------------------------
    //-          PP_RHOATOM                -
    //--------------------------------------
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RHOATOM", true, false))
    {
        ifs.ignore(150, '>');
    }
    else
    {
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RHOATOM>");
    }
    pp.rho_at = std::vector<double>(pp.mesh, 0.0);
    for (int ir = 0; ir < pp.mesh; ir++)
    {
        ifs >> pp.rho_at[ir];
    }
    this->skip_number(ifs, this->mesh_changed);
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RHOATOM>");

    //--------------------------------------
    //-          PP_SPIN_ORB               -
    //--------------------------------------
    if (pp.has_so)
    {
        this->read_pseudo_upf201_so(ifs, pp);
    }

    ModuleBase::GlobalFunc::SCAN_END(ifs, "</UPF>", false);
    return 0;
}

void Pseudopot_upf::getnameval(std::ifstream& ifs, int& n, std::string* name, std::string* val)
{
    std::string txt, word;
    // get long txt
    ifs >> txt;
    while (ifs >> word)
    {
        size_t wl = word.length() - 1;
        txt = txt + " " + word;
        if (word.substr(wl, 1) == ">")
        {
            break;
        }
    }

    // count number of parameters according to "="
    size_t pos = 0;
    n = 0;
    while (1)
    {
        pos = txt.find("=", pos);
        if (pos == std::string::npos)
            break;
        pos++;
        n++;
    }

    // get name & value
    pos = 0;
    size_t pos2, ll;
    for (int i = 0; i < n; ++i)
    {
        pos2 = txt.find("=", pos);
        for (; pos2 > pos; --pos2) // There may be a space before "=";
        {
            if (txt.substr(pos2 - 1, 1) != " ")
                break;
        }
        ll = pos2 - pos;
        name[i] = txt.substr(pos, ll);
        // std::cout<<i<<" "<<name[i]<<std::endl;
        std::string mark;
        bool findmark = false;
        for (int j = 0; j < 100; ++j) // The mark can be ' or " or .
        {
            mark = txt.substr(pos2, 1);
            pos2++;
            if (mark == "\"" || mark == "\'" || mark == ".")
            {
                findmark = true;
                break;
            }
        }
        if (!findmark)
            ModuleBase::WARNING_QUIT(
                "Pseudopot_upf::getnameval",
                "The values are not in \' or \". Please improve the program in read_pp_upf201.cpp");
        pos = pos2;
        pos2 = txt.find(mark, pos);
        ll = pos2 - pos;
        std::string tmpval = txt.substr(pos, ll);
        tmpval = trim(tmpval);
        val[i] = tmpval;
        pos = pos2 + 1;
        for (int j = 0; j < 100; ++j)
        {
            if (txt.substr(pos, 1) == " " || txt.substr(pos, 1) == ",")
                pos++;
            else
                break;
        }
        // std::cout<<name[i]<<"=\""<<val[i]<<"\""<<std::endl;
    }
    return;
}

void Pseudopot_upf::read_pseudo_upf201_header(std::ifstream& ifs, Atom_pseudo& pp)
{
    if (!ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_HEADER"))
        ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_upf201_header", "Found no PP_HEADER");
    std::string name[50];
    std::string val[50];
    int nparameter;
    this->getnameval(ifs, nparameter, name, val);

    for (int ip = 0; ip < nparameter; ++ip)
    {
        if (name[ip] == "generated")
        {
            // add something//
        }
        else if (name[ip] == "author")
        {
        }
        else if (name[ip] == "date")
        {
        }
        else if (name[ip] == "comment")
        {
        }
        else if (name[ip] == "element")
        {
            pp.psd = val[ip];
        }
        else if (name[ip] == "pseudo_type")
        {
            pp.pp_type = val[ip];
            if (pp.pp_type == "SL")
            {
                ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_upf201_header",
                                         "SEMI-LOCAL PSEUDOPOTENTIAL IS NOT SUPPORTED");
            }
        }
        else if (name[ip] == "relativistic")
        {
            relativistic = val[ip];
        }
        else if (name[ip] == "is_ultrasoft")
        {
            if (val[ip] == "T" || val[ip] == "TRUE" || val[ip] == "True" || val[ip] == "true")
            {
                pp.tvanp = true;
            }
            else
            {
                pp.tvanp = false;
            }
        }
        else if (name[ip] == "is_paw")
        {
            if (val[ip] == "T" || val[ip] == "TRUE" || val[ip] == "True" || val[ip] == "true")
            {
                ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_upf201_header", "PAW POTENTIAL IS NOT SUPPORTED");
            }
        }
        else if (name[ip] == "is_coulomb")
        {
            if (val[ip] == "T" || val[ip] == "TRUE" || val[ip] == "True" || val[ip] == "true")
            {
                this->coulomb_potential = true;
            }
        }
        else if (name[ip] == "has_so")
        {
            if (val[ip] == "T" || val[ip] == "TRUE" || val[ip] == "True" || val[ip] == "true")
                pp.has_so = true;
            else
                pp.has_so = false;
        }
        else if (name[ip] == "has_wfc")
        {
            // if (val[ip] == "T" || val[ip] == "TRUE" || val[ip] == "True" || val[ip] == "true")
            // {
            //     has_wfc = true;
            // }
            // else
            // {
            //     has_wfc = false;
            // }
        }
        else if (name[ip] == "has_gipaw")
        {
        }
        else if (name[ip] == "paw_as_gipaw")
        {
        }
        else if (name[ip] == "core_correction")
        {
            if (val[ip] == "T" || val[ip] == "TRUE" || val[ip] == "True" || val[ip] == "true")
                pp.nlcc = true;
            else
                pp.nlcc = false;
        }
        else if (name[ip] == "functional")
        {
            pp.xc_func = val[ip];
        }
        else if (name[ip] == "z_valence")
        {
            pp.zv = std::stod(val[ip]);
        }
        else if (name[ip] == "total_psenergy")
        {
            pp.etotps = atof(val[ip].c_str());
        }
        else if (name[ip] == "wfc_cutoff")
        {
            pp.ecutwfc = atof(val[ip].c_str());
        }
        else if (name[ip] == "rho_cutoff")
        {
            pp.ecutrho = atof(val[ip].c_str());
        }
        else if (name[ip] == "l_max")
        {
            pp.lmax = atoi(val[ip].c_str());
        }
        else if (name[ip] == "l_max_rho")
        {
            this->lmax_rho = atoi(val[ip].c_str());
        }
        else if (name[ip] == "l_local")
        {
            this->lloc = atoi(val[ip].c_str());
        }
        else if (name[ip] == "mesh_size")
        {
            pp.mesh = atoi(val[ip].c_str());
            this->mesh_changed = false;
            if (pp.mesh % 2 == 0)
            {
                pp.mesh -= 1;
                this->mesh_changed = true;
            }
        }
        else if (name[ip] == "number_of_wfc")
        {
            pp.nchi = atoi(val[ip].c_str());
        }
        else if (name[ip] == "number_of_proj")
        {
            pp.nbeta = atoi(val[ip].c_str());
        }
        else
        {
            std::string warningstr
                = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
            ModuleBase::WARNING("PP_HEADRER reading", warningstr);
        }
    }
    if (this->coulomb_potential)
    {
        pp.nbeta = 0;
        pp.lmax = 0;
        this->lloc = 0;
    }
}

void Pseudopot_upf::read_pseudo_upf201_mesh(std::ifstream& ifs, Atom_pseudo& pp)
{
    std::string name[50];
    std::string val[50];
    int nparameter;
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_MESH", true, false))
    {
        this->getnameval(ifs, nparameter, name, val);
        for (int ip = 0; ip < nparameter; ++ip)
        {
            if (name[ip] == "dx")
            {
                dx = atof(val[ip].c_str());
            }
            else if (name[ip] == "mesh")
            {
                pp.mesh = atoi(val[ip].c_str());
                this->mesh_changed = false;
                if (pp.mesh % 2 == 0)
                {
                    pp.mesh -= 1;
                    this->mesh_changed = true;
                }
            }
            else if (name[ip] == "xmin")
            {
                xmin = atof(val[ip].c_str());
            }
            else if (name[ip] == "rmax")
            {
                rmax = atof(val[ip].c_str());
            }
            else if (name[ip] == "zmesh")
            {
                zmesh = atof(val[ip].c_str());
            }
            else
            {
                std::string warningstr
                    = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
                ModuleBase::WARNING("PP_MESH reading", warningstr);
            }
        }
    }
    else
    {
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_MESH>");
    }

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R", true, false))
    {
        ifs.ignore(150, '>');
    }
    else
    {
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R>");
    }
    assert(pp.mesh > 0);
    pp.r = std::vector<double>(pp.mesh, 0.0);
    pp.rab = std::vector<double>(pp.mesh, 0.0);
    for (int ir = 0; ir < pp.mesh; ir++)
    {
        ifs >> pp.r[ir];
    }
    this->skip_number(ifs, this->mesh_changed);
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_R>");

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB", true, false))
    {
        ifs.ignore(150, '>');
    }
    else
    {
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB>");
    }
    for (int ir = 0; ir < pp.mesh; ir++)
    {
        ifs >> pp.rab[ir];
    }
    this->skip_number(ifs, this->mesh_changed);
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RAB>");
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_MESH>");
}

void Pseudopot_upf::read_pseudo_upf201_nonlocal(std::ifstream& ifs, Atom_pseudo& pp)
{
    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
    if (pp.nbeta == 0)
    {
        return;
    }
    std::string word;
    std::string name[50];
    std::string val[50];
    int nparameter;
    this->kbeta = std::vector<int>(pp.nbeta);
    pp.lll = std::vector<int>(pp.nbeta);
    this->els_beta = std::vector<std::string>(pp.nbeta);
    this->rcut = std::vector<double>(pp.nbeta, 0.0);
    this->rcutus = std::vector<double>(pp.nbeta, 0.0);
    pp.betar.create(pp.nbeta, pp.mesh);
    pp.dion.create(pp.nbeta, pp.nbeta);
    for (int ib = 0; ib < pp.nbeta; ib++)
    {
        word = "<PP_BETA." + std::to_string(ib + 1);
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, word);
        this->getnameval(ifs, nparameter, name, val);
        // default value
        els_beta[ib] = "Xn";
        this->kbeta[ib] = pp.mesh;
        rcut[ib] = 0.0;
        rcutus[ib] = 0.0;
        for (int ip = 0; ip < nparameter; ++ip)
        {
            if (name[ip] == "type")
            {
            }
            else if (name[ip] == "size")
            {
            }
            else if (name[ip] == "columns")
            {
            }
            else if (name[ip] == "index")
            {
            }
            else if (name[ip] == "label")
            {
                els_beta[ib] = val[ip];
            }
            else if (name[ip] == "angular_momentum")
            {
                pp.lll[ib] = atoi(val[ip].c_str());
            }
            else if (name[ip] == "cutoff_radius_index")
            {
                this->kbeta[ib] = atoi(val[ip].c_str());
            }
            else if (name[ip] == "cutoff_radius")
            {
                rcut[ib] = atof(val[ip].c_str());
            }
            else if (name[ip] == "ultrasoft_cutoff_radius")
            {
                rcutus[ib] = atof(val[ip].c_str());
            }
            else
            {
                std::string warningstr
                    = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
                ModuleBase::WARNING("PP_BETA reading", warningstr);
            }
        }
        for (int ir = 0; ir < pp.mesh; ir++)
        {
            ifs >> pp.betar(ib, ir);
        }
        this->skip_number(ifs, this->mesh_changed);
        word = "</PP_BETA." + std::to_string(ib + 1) + ">";
        ModuleBase::GlobalFunc::SCAN_END(ifs, word);
    }

    // Read the hamiltonian terms D_ij
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_DIJ", true, false))
    {
        ifs.ignore(150, '>');
    }
    else
    {
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_DIJ>");
    }
    this->nd = pp.nbeta * pp.nbeta;
    for (int i = 0; i < pp.nbeta; i++)
    {
        for (int j = 0; j < pp.nbeta; j++)
        {
            ifs >> pp.dion(i, j);
        }
    }
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_DIJ>");

    // Read the augmentation charge section need by uspp
    if (pp.tvanp && ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_AUGMENTATION"))
    {
        this->getnameval(ifs, nparameter, name, val);
        // default value
        pp.nqlc = 2 * pp.lmax + 1;
        for (int ip = 0; ip < nparameter; ++ip)
        {
            if (name[ip] == "q_with_l")
            {
                if (val[ip] == "T" || val[ip] == "TRUE" || val[ip] == "True" || val[ip] == "true")
                {
                    q_with_l = true;
                }
                else
                {
                    q_with_l = false;
                }
            }
            else if (name[ip] == "nqf")
            {
                nqf = atoi(val[ip].c_str());
            }
            else if (name[ip] == "nqlc")
            {
                pp.nqlc = atoi(val[ip].c_str());
            }
            else
            {
                std::string warningstr
                    = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
                ModuleBase::WARNING("PP_AUGMENTATION reading", warningstr);
            }
        }

        // PP_Q
        if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_Q", true, false))
        {
            ifs.ignore(150, '>');
        }
        else
        {
            ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_Q>");
        }
        pp.qqq.create(pp.nbeta, pp.nbeta);
        for (int i = 0; i < pp.nbeta; i++)
        {
            for (int j = 0; j < pp.nbeta; j++)
            {
                ifs >> pp.qqq(i, j);
            }
        }
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_Q>");

        // Read polinomial coefficients for Q_ij expansion at small radius
        this->rinner = std::vector<double>(pp.nqlc, 0.0);
        if (nqf <= 0)
        {
            this->qfcoef.create(1, 1, 1, 1);
        }
        else
        {
            this->qfcoef.create(pp.nbeta, pp.nbeta, pp.nqlc, nqf);
            if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_QFCOEF>"))
            {
                for (int i = 0; i < pp.nbeta; i++)
                {
                    for (int j = 0; j < pp.nbeta; j++)
                    {
                        for (int k = 0; k < pp.nqlc; k++)
                        {
                            for (int l = 0; l < nqf; l++)
                            {
                                ifs >> qfcoef(i, j, k, l);
                            }
                        }
                    }
                }
            }
            ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_QFCOEF>");
            if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RINNER>"))
            {
                for (int i = 0; i < pp.nqlc; i++)
                {
                    ifs >> rinner[i];
                }
            }
            ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RINNER>");
        }

        // Read augmentation charge Q_ij
        if (q_with_l)
        {
            pp.qfuncl.create(2 * pp.lmax + 1, pp.nbeta * (pp.nbeta + 1) / 2, pp.mesh);
        }
        else
        {
            this->qfunc.create(pp.nbeta * (pp.nbeta + 1) / 2, pp.mesh);
        }

        for (int nb = 0; nb < pp.nbeta; nb++)
        {
            int ln = pp.lll[nb];
            for (int mb = nb; mb < pp.nbeta; mb++)
            {
                int lm = pp.lll[mb];
                int nmb = mb * (mb + 1) / 2 + nb;
                if (q_with_l)
                {
                    for (int l = std::abs(ln - lm); l <= ln + lm; l += 2)
                    {
                        word = "<PP_QIJL." + std::to_string(nb + 1) + "." + std::to_string(mb + 1) + "."
                               + std::to_string(l);
                        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, word);
                        ifs.ignore(150, '>');
                        for (int ir = 0; ir < pp.mesh; ir++)
                        {
                            ifs >> pp.qfuncl(l, nmb, ir);
                        }
                        this->skip_number(ifs, this->mesh_changed);
                        word = "</PP_QIJL." + std::to_string(nb + 1) + "." + std::to_string(mb + 1) + "."
                               + std::to_string(l) + ">";
                        ModuleBase::GlobalFunc::SCAN_END(ifs, word);
                    }
                }
                else
                {
                    word = "<PP_QIJ." + std::to_string(nb + 1) + "." + std::to_string(mb + 1);
                    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, word);
                    ifs.ignore(150, '>');
                    for (int ir = 0; ir < pp.mesh; ir++)
                    {
                        ifs >> this->qfunc(nmb, ir);
                    }
                    this->skip_number(ifs, this->mesh_changed);
                    word = "</PP_QIJ." + std::to_string(nb + 1) + "." + std::to_string(mb + 1) + ">";
                    ModuleBase::GlobalFunc::SCAN_END(ifs, word);
                }
            }
        }
        ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_AUGMENTATION>");
    }
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NONLOCAL>");

    pp.kkbeta = 0;
    for (int nb = 0; nb < pp.nbeta; nb++)
    {
        pp.kkbeta = (this->kbeta[nb] > pp.kkbeta) ? this->kbeta[nb] : pp.kkbeta;
    }
}

void Pseudopot_upf::read_pseudo_upf201_pswfc(std::ifstream& ifs, Atom_pseudo& pp)
{
    std::string word;
    std::string name[50];
    std::string val[50];
    int nparameter;
    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_PSWFC>");
    pp.els = std::vector<std::string>(pp.nchi, "");
    pp.lchi = std::vector<int>(pp.nchi, 0);
    this->nchi = std::vector<int>(pp.nchi, 0);
    pp.oc = std::vector<double>(pp.nchi, 0.0);
    this->epseu = std::vector<double>(pp.nchi, 0.0);
    this->rcut_chi = std::vector<double>(pp.nchi, 0.0);
    this->rcutus_chi = std::vector<double>(pp.nchi, 0.0);
    pp.chi.create(pp.nchi, pp.mesh);
    for (int iw = 0; iw < pp.nchi; iw++)
    {
        // default value
        pp.els[iw] = "Xn";
        nchi[iw] = -1;
        epseu[iw] = 0.0;
        rcut_chi[iw] = 0.0;
        rcutus_chi[iw] = 0.0;

        word = "<PP_CHI." + std::to_string(iw + 1);
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, word);
        this->getnameval(ifs, nparameter, name, val);
        for (int ip = 0; ip < nparameter; ++ip)
        {
            if (name[ip] == "type")
            {
            }
            else if (name[ip] == "size")
            {
            }
            else if (name[ip] == "columns")
            {
            }
            else if (name[ip] == "index")
            {
            }
            else if (name[ip] == "label")
            {
                pp.els[iw] = val[ip];
            }
            else if (name[ip] == "l")
            {
                pp.lchi[iw] = atoi(val[ip].c_str());
                if (nchi[iw] == -1)
                {
                    nchi[iw] = pp.lchi[iw] - 1;
                }
            }
            else if (name[ip] == "occupation")
            {
                pp.oc[iw] = atof(val[ip].c_str());
            }
            else if (name[ip] == "n")
            {
                nchi[iw] = atoi(val[ip].c_str());
            }
            else if (name[ip] == "pseudo_energy")
            {
                epseu[iw] = atof(val[ip].c_str());
            }
            else if (name[ip] == "cutoff_radius")
            {
                rcut_chi[iw] = atof(val[ip].c_str());
            }
            else if (name[ip] == "ultrasoft_cutoff_radius")
            {
                rcutus_chi[iw] = atof(val[ip].c_str());
            }
            else
            {
                std::string warningstr
                    = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
                ModuleBase::WARNING("PP_CHI reading", warningstr);
            }
        }
        for (int ir = 0; ir < pp.mesh; ir++)
        {
            ifs >> pp.chi(iw, ir);
        }
        for (int ir = 0; ir < pp.mesh; ir++)
        {
            assert(pp.chi.c[iw * pp.mesh + ir] == pp.chi(iw, ir));
        }
        this->skip_number(ifs, this->mesh_changed);
        word = "</PP_CHI." + std::to_string(iw + 1) + ">";
        ModuleBase::GlobalFunc::SCAN_END(ifs, word);
    }
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_PSWFC>");
}

/*
// only used in paw?
void Pseudopot_upf::read_pseudo_upf201_fullwfc(std::ifstream& ifs)
{
    std::string word;
    std::string name[50];
    std::string val[50];
    int nparameter;
    this->aewfc.create(this->nbeta, this->mesh);
    this->pswfc.create(this->nbeta, this->mesh);
    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_FULL_WFC");
    this->getnameval(ifs, nparameter, name, val);
    ifs.ignore(150, '>');
    for (int ib = 0; ib < nbeta; ib++)
    {
        // All-electron wavefunctions corresponding to beta functions
        word = "<PP_AEWFC." + std::to_string(ib + 1);
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, word);
        ifs.ignore(150, '>');
        for (int ir = 0; ir < mesh; ir++)
        {
            ifs >> this->aewfc(ib, ir);
        }
        word = "</PP_AEWFC." + std::to_string(ib + 1) + ">";
        ModuleBase::GlobalFunc::SCAN_END(ifs, word);

        // Pseudo wavefunctions corresponding to beta functions
        word = "<PP_PSWFC." + std::to_string(ib + 1);
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, word);
        ifs.ignore(150, '>');
        for (int ir = 0; ir < mesh; ir++)
        {
            ifs >> this->pswfc(ib, ir);
        }
        word = "</PP_PSWFC." + std::to_string(ib + 1) + ">";
        ModuleBase::GlobalFunc::SCAN_END(ifs, word);
    }
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_FULL_WFC>");
}*/

void Pseudopot_upf::read_pseudo_upf201_so(std::ifstream& ifs, Atom_pseudo& pp)
{
    std::string word;
    std::string name[50];
    std::string val[50];
    int nparameter;
    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_SPIN_ORB>");
    pp.jchi = std::vector<double>(pp.nchi, 0.0);
    pp.jjj = std::vector<double>(pp.nbeta, 0.0);
    pp.nn = std::vector<int>(pp.nchi, 0);

    for (int nw = 0; nw < pp.nchi; nw++)
    {
        word = "<PP_RELWFC." + std::to_string(nw + 1);
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, word);
        this->getnameval(ifs, nparameter, name, val);
        for (int ip = 0; ip < nparameter; ++ip)
        {
            if (name[ip] == "index")
            {
            }
            else if (name[ip] == "els")
            {
                pp.els[nw] = val[ip];
            }
            else if (name[ip] == "nn")
            {
                pp.nn[nw] = atoi(val[ip].c_str());
            }
            else if (name[ip] == "lchi")
            {
                pp.lchi[nw] = atoi(val[ip].c_str());
            }
            else if (name[ip] == "jchi")
            {
                pp.jchi[nw] = atof(val[ip].c_str());
            }
            else if (name[ip] == "oc")
            {
                pp.oc[nw] = atof(val[ip].c_str());
            }
            else
            {
                std::string warningstr
                    = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
                ModuleBase::WARNING("PP_RELWFC reading", warningstr);
            }
        }
        word = "</PP_RELWFC." + std::to_string(nw + 1) + ">";
        ModuleBase::GlobalFunc::SCAN_END(ifs, word);
    }

    for (int nb = 0; nb < pp.nbeta; nb++)
    {
        word = "<PP_RELBETA." + std::to_string(nb + 1);
        ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, word);
        this->getnameval(ifs, nparameter, name, val);
        for (int ip = 0; ip < nparameter; ++ip)
        {
            if (name[ip] == "index")
            {
            }
            else if (name[ip] == "lll")
            {
                pp.lll[nb] = atoi(val[ip].c_str());
            }
            else if (name[ip] == "jjj")
            {
                pp.jjj[nb] = atof(val[ip].c_str());
            }
            else
            {
                std::string warningstr
                    = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
                ModuleBase::WARNING("PP_RELBETA reading", warningstr);
            }
        }
        word = "</PP_RELBETA." + std::to_string(nb + 1) + ">";
        ModuleBase::GlobalFunc::SCAN_END(ifs, word);
    }
    ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_SPIN_ORB>");
}