#include "module_io/output_mulliken.h"

#include "module_base/formatter.h"
#include "module_base/name_angular.h"
#include "module_base/scalapack_connector.h"
#include "module_base/tool_quit.h"

#include <numeric>

namespace ModuleIO
{

template <typename TK>
Output_Mulliken<TK>::Output_Mulliken(Output_Sk<TK>* output_sk,
                                     Output_DMK<TK>* output_dmk,
                                     Parallel_Orbitals* ParaV,
                                     CellIndex* cell_index,
                                     const std::vector<int>& isk,
                                     int nspin)
    : output_sk_(output_sk), output_dmk_(output_dmk), ParaV_(ParaV), cell_index_(cell_index), isk_(isk), nspin_(nspin)
{
    this->set_nspin(nspin);
    this->set_ParaV(ParaV);
    this->cal_orbMulP();
}

template <typename TK>
void Output_Mulliken<TK>::write(int istep, std::string out_dir)
{
    std::vector<double> tot_chg = this->get_tot_chg();
    std::vector<std::vector<double>> atom_chg = this->get_atom_chg();
    std::map<std::vector<int>, double> orb_chg = this->get_orb_chg();
    std::stringstream as;
    as << out_dir << "mulliken.txt";
    std::ofstream os;
    if (istep == 0)
    {
        os.open(as.str(), std::ios::out);
    }
    else
    {
        os.open(as.str(), std::ios::app);
    }
    if (this->nspin_ == 1)
    {
        this->write_mulliken_nspin1(istep, tot_chg, atom_chg, orb_chg, os);
    }
    else if (this->nspin_ == 2)
    {
        this->write_mulliken_nspin2(istep, tot_chg, atom_chg, orb_chg, os);
    }
    else if (this->nspin_ == 4)
    {
        this->write_mulliken_nspin4(istep, tot_chg, atom_chg, orb_chg, os);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Output_Mulliken::write", "nspin must be 1, 2 or 4");
    }
    os.close();
}

template <typename TK>
void Output_Mulliken<TK>::write_mulliken_nspin1(int istep,
                                                const std::vector<double>& tot_chg,
                                                const std::vector<std::vector<double>>& atom_chg,
                                                std::map<std::vector<int>, double> orb_chg,
                                                std::ofstream& os)
{
    os << std::setprecision(4);
    /// step info
    os << "STEP: " << istep << std::endl;
    os << "CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM" << std::endl;
    os << " Total charge:\t" << tot_chg[0] << std::endl;
    /// orbital decomposed mulliken populations
    FmtCore fmt_of_chg("%20.4f");
    FmtCore fmt_of_label("%-20s");
    FmtCore fmt_of_Z("%20d");
    os << "Decomposed Mulliken populations" << std::endl;
    for (int iat = 0; iat < this->cell_index_->get_nat(); ++iat)
    {
        /// header of the table
        std::string atom_label = this->cell_index_->get_atom_label(iat);
        os << FmtCore::format("%-20d", iat) << FmtCore::format("%20s", std::string("Zeta of ") + atom_label)
           << FmtCore::format("%20s", std::string("Spin 1")) << std::endl;
        /// loop of L
        for (int L = 0; L <= this->cell_index_->get_maxL(iat); L++)
        {
            std::vector<double> sum_over_m_and_z(this->nspin_, 0.0);
            for (int Z = 0; Z < this->cell_index_->get_nchi(iat, L); Z++)
            {
                for (int M = 0; M < (2 * L + 1); M++)
                {
                    os << fmt_of_label.format(ModuleBase::Name_Angular[L][M]) << fmt_of_Z.format(Z)
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 0, L, Z, M}]) << std::endl;
                }
                // sum over m
                std::vector<double> sum_over_m(this->nspin_, 0.0);
                for (int is = 0; is < this->nspin_; is++)
                {
                    for (int M = 0; M < (2 * L + 1); M++)
                    {
                        sum_over_m[is] += orb_chg[std::vector<int>{iat, is, L, Z, M}];
                    }
                    sum_over_m_and_z[is] += sum_over_m[is];
                }
                if (L > 0)
                {
                    os << fmt_of_label.format(std::string("SUM OVER M")) << std::setw(20) << ""
                       << fmt_of_chg.format(sum_over_m[0]) << std::endl;
                }
            }
            os << fmt_of_label.format(std::string("SUM OVER M+Zeta")) << std::setw(20) << ""
               << fmt_of_chg.format(sum_over_m_and_z[0]) << std::endl;
            os << std::endl;
        }
        os << fmt_of_label.format(std::string("SUM OVER M+Zeta+L")) << std::setw(20) << ""
           << fmt_of_chg.format(atom_chg[iat][0]) << std::endl;
        os << std::endl;
        os << std::left << std::setw(30) << "Total Charge on atom:" << std::right << std::setw(10) << atom_label
           << fmt_of_chg.format(atom_chg[iat][0]) << std::endl;
        os << std::endl << std::endl;
    }
}

template <typename TK>
void Output_Mulliken<TK>::write_mulliken_nspin2(int istep,
                                                const std::vector<double>& tot_chg,
                                                const std::vector<std::vector<double>>& atom_chg,
                                                std::map<std::vector<int>, double> orb_chg,
                                                std::ofstream& os)
{
    os << std::setprecision(4);
    /// step info
    os << "STEP: " << istep << std::endl;
    os << "CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM" << std::endl;
    os << " Total charge of spin " << 1 << ":\t" << tot_chg[0] << std::endl;
    os << " Total charge of spin " << 2 << ":\t" << tot_chg[1] << std::endl;
    os << " Total charge:\t" << tot_chg[0] + tot_chg[1] << std::endl;
    /// orbital decomposed mulliken populations
    FmtCore fmt_of_chg("%20.4f");
    FmtCore fmt_of_label("%-20s");
    FmtCore fmt_of_Z("%20d");
    os << "Decomposed Mulliken populations" << std::endl;
    for (int iat = 0; iat < this->cell_index_->get_nat(); ++iat)
    {
        /// header of the table
        std::string atom_label = this->cell_index_->get_atom_label(iat);
        os << FmtCore::format("%-20d", iat) << FmtCore::format("%20s", std::string("Zeta of ") + atom_label)
           << FmtCore::format("%20s", std::string("Spin 1")) << FmtCore::format("%20s", std::string("Spin 2"))
           << FmtCore::format("%20s", std::string("Sum")) << FmtCore::format("%20s", std::string("Diff")) << std::endl;
        /// loop of L
        for (int L = 0; L <= this->cell_index_->get_maxL(iat); L++)
        {
            std::vector<double> sum_over_m_and_z(this->nspin_, 0.0);
            for (int Z = 0; Z < this->cell_index_->get_nchi(iat, L); Z++)
            {
                for (int M = 0; M < (2 * L + 1); M++)
                {
                    os << fmt_of_label.format(ModuleBase::Name_Angular[L][M]) << fmt_of_Z.format(Z)
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 0, L, Z, M}])
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 1, L, Z, M}])
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 0, L, Z, M}]
                                            + orb_chg[std::vector<int>{iat, 1, L, Z, M}])
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 0, L, Z, M}]
                                            - orb_chg[std::vector<int>{iat, 1, L, Z, M}])
                       << std::endl;
                }
                // sum over m
                std::vector<double> sum_over_m(this->nspin_, 0.0);
                for (int is = 0; is < this->nspin_; is++)
                {
                    for (int M = 0; M < (2 * L + 1); M++)
                    {
                        sum_over_m[is] += orb_chg[std::vector<int>{iat, is, L, Z, M}];
                    }
                    sum_over_m_and_z[is] += sum_over_m[is];
                }
                if (L > 0)
                {
                    os << fmt_of_label.format(std::string("SUM OVER M")) << std::setw(20) << ""
                       << fmt_of_chg.format(sum_over_m[0]) << fmt_of_chg.format(sum_over_m[1])
                       << fmt_of_chg.format(sum_over_m[0] + sum_over_m[1])
                       << fmt_of_chg.format(sum_over_m[0] - sum_over_m[1]) << std::endl;
                }
            }
            os << fmt_of_label.format(std::string("SUM OVER M+Zeta")) << std::setw(20) << ""
               << fmt_of_chg.format(sum_over_m_and_z[0]) << fmt_of_chg.format(sum_over_m_and_z[1])
               << fmt_of_chg.format(sum_over_m_and_z[0] + sum_over_m_and_z[1])
               << fmt_of_chg.format(sum_over_m_and_z[0] - sum_over_m_and_z[1]) << std::endl;
            os << std::endl;
        }
        os << fmt_of_label.format(std::string("SUM OVER M+Zeta+L")) << std::setw(20) << ""
           << fmt_of_chg.format(atom_chg[iat][0]) << fmt_of_chg.format(atom_chg[iat][1])
           << fmt_of_chg.format(atom_chg[iat][0] + atom_chg[iat][1])
           << fmt_of_chg.format(atom_chg[iat][0] - atom_chg[iat][1]) << std::endl;
        os << std::endl;
        os << std::left << std::setw(30) << "Total Charge on atom:" << std::right << std::setw(10) << atom_label
           << fmt_of_chg.format(atom_chg[iat][0] + atom_chg[iat][1]) << std::endl;
        os << std::left << std::setw(30) << "Total Magnetism on atom: " << std::right << std::setw(10) << atom_label
           << fmt_of_chg.format(atom_chg[iat][0] - atom_chg[iat][1]) << std::endl;
        os << std::endl << std::endl;
    }
}

template <typename TK>
void Output_Mulliken<TK>::write_mulliken_nspin4(int istep,
                                                const std::vector<double>& tot_chg,
                                                const std::vector<std::vector<double>>& atom_chg,
                                                std::map<std::vector<int>, double> orb_chg,
                                                std::ofstream& os)
{
    os << std::setprecision(4);
    /// step info
    os << "STEP: " << istep << std::endl;
    os << "CALCULATE THE MULLIkEN ANALYSIS FOR EACH ATOM" << std::endl;
    os << " Total charge:\t" << tot_chg[0] << std::endl;
    /// orbital decomposed mulliken populations
    FmtCore fmt_of_chg("%20.4f");
    FmtCore fmt_of_label("%-20s");
    FmtCore fmt_of_Z("%20d");
    os << "Decomposed Mulliken populations" << std::endl;
    for (int iat = 0; iat < this->cell_index_->get_nat(); ++iat)
    {
        /// header of the table
        std::string atom_label = this->cell_index_->get_atom_label(iat);
        os << FmtCore::format("%-20d", iat) << FmtCore::format("%20s", std::string("Zeta of ") + atom_label)
           << FmtCore::format("%20s", std::string("Spin 1")) << FmtCore::format("%20s", std::string("Spin 2"))
           << FmtCore::format("%20s", std::string("Spin 3")) << FmtCore::format("%20s", std::string("Spin 4"))
           << std::endl;
        /// loop of L
        for (int L = 0; L <= this->cell_index_->get_maxL(iat); L++)
        {
            std::vector<double> sum_over_m_and_z(this->nspin_, 0.0);
            for (int Z = 0; Z < this->cell_index_->get_nchi(iat, L); Z++)
            {
                for (int M = 0; M < (2 * L + 1); M++)
                {
                    os << fmt_of_label.format(ModuleBase::Name_Angular[L][M]) << fmt_of_Z.format(Z)
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 0, L, Z, M}])
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 1, L, Z, M}])
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 2, L, Z, M}])
                       << fmt_of_chg.format(orb_chg[std::vector<int>{iat, 3, L, Z, M}]) << std::endl;
                }
                // sum over m
                std::vector<double> sum_over_m(this->nspin_, 0.0);
                for (int is = 0; is < this->nspin_; is++)
                {
                    for (int M = 0; M < (2 * L + 1); M++)
                    {
                        sum_over_m[is] += orb_chg[std::vector<int>{iat, is, L, Z, M}];
                    }
                    sum_over_m_and_z[is] += sum_over_m[is];
                }
                if (L > 0)
                {
                    os << fmt_of_label.format(std::string("SUM OVER M")) << std::setw(20) << ""
                       << fmt_of_chg.format(sum_over_m[0]) << fmt_of_chg.format(sum_over_m[1])
                       << fmt_of_chg.format(sum_over_m[2]) << fmt_of_chg.format(sum_over_m[3]) << std::endl;
                }
            }
            os << fmt_of_label.format(std::string("SUM OVER M+Zeta")) << std::setw(20) << ""
               << fmt_of_chg.format(sum_over_m_and_z[0]) << fmt_of_chg.format(sum_over_m_and_z[1])
               << fmt_of_chg.format(sum_over_m_and_z[2]) << fmt_of_chg.format(sum_over_m_and_z[3]) << std::endl;
            os << std::endl;
        }
        os << fmt_of_label.format(std::string("SUM OVER M+Zeta+L")) << std::setw(20) << ""
           << fmt_of_chg.format(atom_chg[iat][0]) << fmt_of_chg.format(atom_chg[iat][1])
           << fmt_of_chg.format(atom_chg[iat][2]) << fmt_of_chg.format(atom_chg[iat][3]) << std::endl;
        os << std::endl;
        os << std::left << std::setw(30) << "Total Charge on atom:" << std::right << std::setw(10) << atom_label
           << fmt_of_chg.format(atom_chg[iat][0]) << std::endl;
        os << std::left << std::setw(30) << "Total Magnetism on atom: " << std::right << std::setw(10) << atom_label
           << fmt_of_chg.format(atom_chg[iat][1]) << fmt_of_chg.format(atom_chg[iat][2])
           << fmt_of_chg.format(atom_chg[iat][3]) << std::endl;
        os << std::endl << std::endl;
    }
}

/// set nspin
template <typename TK>
void Output_Mulliken<TK>::set_nspin(int nspin_in)
{
    if (nspin_in != 1 && nspin_in != 2 && nspin_in != 4)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_nspin", "nspin must be 1, 2 or 4");
    }
    this->nspin_ = nspin_in;
}

/// @brief  set ParaV
template <typename TK>
void Output_Mulliken<TK>::set_ParaV(Parallel_Orbitals* ParaV_in)
{
    this->ParaV_ = ParaV_in;
    int nloc = this->ParaV_->nloc;
    if (nloc <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_ParaV", "nloc <= 0");
    }
}

template <typename TK>
std::vector<double> Output_Mulliken<TK>::get_tot_chg()
{
    std::vector<double> tot_chg(this->nspin_, 0.0);
    int nw = this->cell_index_->get_nw();
    const int nlocal = (this->nspin_ == 4) ? nw / 2 : nw;
    for (int is = 0; is != this->nspin_; ++is)
    {
        for (size_t iw = 0; iw != nlocal; ++iw)
        {
            tot_chg[is] += this->orbMulP_(is, iw);
        }
    }
    return tot_chg;
}

template <typename TK>
std::vector<std::vector<double>> Output_Mulliken<TK>::get_atom_chg()
{
    int nat = this->cell_index_->get_nat();
    std::vector<std::vector<double>> atom_chg(nat, std::vector<double>(this->nspin_, 0.0));
    int num = 0;
    for (int iat = 0; iat < nat; iat++)
    {
        int nw_it = this->cell_index_->get_nw(iat);
        for (int iw = 0; iw < nw_it; iw++)
        {
            for (int is = 0; is < this->nspin_; is++)
            {
                atom_chg[iat][is] += std::abs(this->orbMulP_(is, num)) < 1e-10 ? 0.0 : this->orbMulP_(is, num);
            }
            num++;
        }
    }
    return atom_chg;
}

template <typename TK>
std::map<std::vector<int>, double> Output_Mulliken<TK>::get_orb_chg()
{
    int nat = this->cell_index_->get_nat();
    std::map<std::vector<int>, double> orb_chg;
    for (int is = 0; is < this->nspin_; is++)
    {
        int num = 0;
        for (int iat = 0; iat < nat; iat++)
        {
            int maxL = this->cell_index_->get_maxL(iat);
            for (int L = 0; L <= maxL; L++)
            {
                int nchi = this->cell_index_->get_nchi(iat, L);
                for (int Z = 0; Z < this->cell_index_->get_nchi(iat, L); Z++)
                {
                    for (int M = 0; M < (2 * L + 1); M++)
                    {
                        orb_chg[std::vector<int>{iat, is, L, Z, M}]
                            = std::abs(this->orbMulP_(is, num)) < 1e-10 ? 0.0 : this->orbMulP_(is, num);
                        num++;
                    }
                }
            }
        }
    }
    return orb_chg;
}

template <typename TK>
void Output_Mulliken<TK>::collect_MW(ModuleBase::matrix& MecMulP, const ModuleBase::ComplexMatrix& mud, int nw, int isk)
{
    if (this->nspin_ == 1 || this->nspin_ == 2)
    {
        for (size_t i = 0; i < nw; ++i)
        {
            if (this->ParaV_->in_this_processor(i, i))
            {
                const int ir = this->ParaV_->global2local_row(i);
                const int ic = this->ParaV_->global2local_col(i);
                MecMulP(isk, i) += mud(ic, ir).real();
            }
        }
    }
    else if (this->nspin_ == 4)
    {
        for (size_t i = 0; i < nw; ++i)
        {
            const int index = i % 2;
            if (!index)
            {
                const int j = i / 2;
                const int k1 = 2 * j;
                const int k2 = 2 * j + 1;
                if (this->ParaV_->in_this_processor(k1, k1))
                {
                    const int ir = this->ParaV_->global2local_row(k1);
                    const int ic = this->ParaV_->global2local_col(k1);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) += mud(ic, ir).real();
                }
                if (this->ParaV_->in_this_processor(k1, k2))
                {
                    const int ir = this->ParaV_->global2local_row(k1);
                    const int ic = this->ParaV_->global2local_col(k2);
                    // note that mud is column major
                    MecMulP(1, j) += mud(ic, ir).real();
                    // M_y = i(M_{up,down} - M_{down,up}) = -(M_{up,down} - M_{down,up}).imag()
                    MecMulP(2, j) -= mud(ic, ir).imag();
                }
                if (this->ParaV_->in_this_processor(k2, k1))
                {
                    const int ir = this->ParaV_->global2local_row(k2);
                    const int ic = this->ParaV_->global2local_col(k1);
                    MecMulP(1, j) += mud(ic, ir).real();
                    // M_y = i(M_{up,down} - M_{down,up}) = -(M_{up,down} - M_{down,up}).imag()
                    MecMulP(2, j) += mud(ic, ir).imag();
                }
                if (this->ParaV_->in_this_processor(k2, k2))
                {
                    const int ir = this->ParaV_->global2local_row(k2);
                    const int ic = this->ParaV_->global2local_col(k2);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) -= mud(ic, ir).real();
                }
            }
        }
    }
}

template <typename TK>
void Output_Mulliken<TK>::print_atom_mag(const std::vector<std::vector<double>>& atom_chg, std::ostream& os)
{
    int nat = this->cell_index_->get_nat();
    std::vector<std::string> atom_label;
    std::vector<double> mag_x(nat, 0.0);
    std::vector<double> mag_y(nat, 0.0);
    std::vector<double> mag_z(nat, 0.0);
    if (this->nspin_ == 2)
    {
        const std::vector<std::string> title = {"Total Magnetism (uB)", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f"};
        FmtTable table(title, nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            atom_label.push_back(this->cell_index_->get_atom_label(iat, true));
            mag_z[iat] = atom_chg[iat][0] - atom_chg[iat][1];
        }
        table << atom_label << mag_z;
        os << table.str() << std::endl;
    }
    else if (this->nspin_ == 4)
    {
        std::vector<double> magnitude(nat, 0.0);
        std::vector<double> polar(nat, 0.0);
        std::vector<double> azimuth(nat, 0.0);
        const std::vector<std::string> title = {"Total Magnetism (uB)", "x", "y", "z"};
        const std::vector<std::string> fmts = {"%26s", "%20.10f", "%20.10f", "%20.10f"};
        FmtTable table(title, nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::RIGHT});
        for (int iat = 0; iat < nat; ++iat)
        {
            atom_label.push_back(this->cell_index_->get_atom_label(iat, true));
            mag_x[iat] = atom_chg[iat][1];
            mag_y[iat] = atom_chg[iat][2];
            mag_z[iat] = atom_chg[iat][3];
            magnitude[iat] = std::sqrt(mag_x[iat] * mag_x[iat] + mag_y[iat] * mag_y[iat] + mag_z[iat] * mag_z[iat]);
            polar[iat] = std::acos(mag_z[iat] / magnitude[iat]) * 180.0 / ModuleBase::PI;
            azimuth[iat] = std::atan2(mag_y[iat], mag_x[iat]) * 180.0 / ModuleBase::PI;
        }
        table << atom_label << mag_x << mag_y << mag_z;
        os << table.str() << std::endl;
        /// output mag in polar coordinates
        const std::vector<std::string> title_polar = {"Total Magnetism (uB)", "Magnitude (uB)", "Polar (degree)", "Azimuth (degree)"};
        const std::vector<std::string> fmts_polar = {"%26s", "%20.10f", "%20.10f", "%20.10f"};
        FmtTable table_polar(title_polar, nat, fmts_polar, {FmtTable::Align::RIGHT, FmtTable::Align::RIGHT});
        table_polar << atom_label << magnitude << polar << azimuth;
        os << table_polar.str() << std::endl;
    }
    else if (this->nspin_ == 1)
    {
        /// do nothing due to no mag info available
    }
    else
    {
        ModuleBase::WARNING_QUIT("Output_Mulliken::print_atom_mag", "nspin must be 1, 2 or 4");
    }
}

template <typename TK>
std::vector<std::vector<double>> Output_Mulliken<TK>::get_atom_mulliken(std::vector<std::vector<double>>& atom_chg)
{
    int nat = this->cell_index_->get_nat();
    std::vector<std::vector<double>> atom_mulliken(nat, std::vector<double>(this->nspin_, 0.0));
    for (int iat = 0; iat < nat; iat++)
    {
        if (this->nspin_ == 1)
        {
            atom_mulliken[iat][0] = atom_chg[iat][0];
        }
        else if (this->nspin_ == 2)
        {
            atom_mulliken[iat][0] = atom_chg[iat][0] + atom_chg[iat][1];
            atom_mulliken[iat][1] = atom_chg[iat][0] - atom_chg[iat][1];
        }
        else if (this->nspin_ == 4)
        {
            atom_mulliken[iat][0] = atom_chg[iat][0];
            atom_mulliken[iat][1] = atom_chg[iat][1];
            atom_mulliken[iat][2] = atom_chg[iat][2];
            atom_mulliken[iat][3] = atom_chg[iat][3];
        }
        else
        {
            ModuleBase::WARNING_QUIT("Output_Mulliken::get_atom_mulliken", "nspin must be 1, 2 or 4");
        }
    }
    return atom_mulliken;
}

template <>
void Output_Mulliken<std::complex<double>>::cal_orbMulP()
{
    ModuleBase::TITLE("module_deltaspin", "cal_MW_k");
    int nw = this->cell_index_->get_nw();
    const int nlocal = (this->nspin_ == 4) ? nw / 2 : nw;
    ModuleBase::matrix MecMulP(this->nspin_, nlocal, true);
    this->orbMulP_.create(this->nspin_, nlocal, true);
    for (size_t ik = 0; ik != this->isk_.size(); ++ik)
    {
        auto p_Sk = this->output_sk_->get_Sk(ik);
        auto p_DMk = this->output_dmk_->get_DMK(ik);
        ModuleBase::ComplexMatrix mud(this->ParaV_->ncol, this->ParaV_->nrow, true);
#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const std::complex<double> one_float = {1.0, 0.0}, zero_float = {0.0, 0.0};
        pzgemm_(&N_char,
                &T_char,
                &nw,
                &nw,
                &nw,
                &one_float,
                p_DMk,
                &one_int,
                &one_int,
                this->ParaV_->desc,
                p_Sk,
                &one_int,
                &one_int,
                this->ParaV_->desc,
                &zero_float,
                mud.c,
                &one_int,
                &one_int,
                this->ParaV_->desc);
        this->collect_MW(MecMulP, mud, nw, this->isk_[ik]);
#endif
    }
#ifdef __MPI
    MPI_Allreduce(MecMulP.c, this->orbMulP_.c, this->nspin_ * nlocal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}

template <>
void Output_Mulliken<double>::cal_orbMulP()
{
    ModuleBase::TITLE("Mulliken_Charge", "cal_mulliken");
    int nw = this->cell_index_->get_nw();
    const int nspin = (this->nspin_ == 2) ? 2 : 1;
    const int nlocal = (this->nspin_ == 4) ? nw / 2 : nw;
    // std::vector<std::vector<double>> MecMulP(GlobalV::NSPIN, std::vector<double>(nlocal, 0));
    // std::vector<std::vector<double>> orbMulP(GlobalV::NSPIN, std::vector<double>(nlocal, 0));
    ModuleBase::matrix MecMulP(this->nspin_, nlocal, true);
    this->orbMulP_.create(this->nspin_, nlocal, true);

    for (size_t is = 0; is != nspin; ++is)
    {
        ModuleBase::matrix mud;
        auto p_Sk = this->output_sk_->get_Sk(is);
        auto p_DMk = this->output_dmk_->get_DMK(is);
        mud.create(this->ParaV_->ncol, this->ParaV_->nrow);
#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const double one_float = 1.0, zero_float = 0.0;
        pdgemm_(&N_char,
                &T_char,
                &nw,
                &nw,
                &nw,
                &one_float,
                p_DMk,
                &one_int,
                &one_int,
                this->ParaV_->desc,
                p_Sk,
                &one_int,
                &one_int,
                this->ParaV_->desc,
                &zero_float,
                mud.c,
                &one_int,
                &one_int,
                this->ParaV_->desc);
        if (this->nspin_ == 1 || this->nspin_ == 2)
        {
            for (size_t i = 0; i != nw; ++i)
                if (this->ParaV_->in_this_processor(i, i))
                {
                    const int ir = this->ParaV_->global2local_row(i);
                    const int ic = this->ParaV_->global2local_col(i);
                    MecMulP(is, i) += mud(ic, ir);
                }
        }
#endif
    }
#ifdef __MPI
    MPI_Reduce(MecMulP.c, this->orbMulP_.c, this->nspin_ * nlocal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}

// case for nspin<4, gamma-only k-point
template class Output_Mulliken<double>;
// case for nspin<4, multi-k-points
template class Output_Mulliken<std::complex<double>>;

} // namespace ModuleIO
