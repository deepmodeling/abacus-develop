#include "source_pw/module_pwdft/onsite_proj_print.h"
#include "source_base/formatter.h"

namespace projectors {
namespace print {

void print_orb_chg(
    const UnitCell* ucell,
    const std::vector<std::complex<double>>& occs,
    const std::vector<int>& iat_nh,
    const std::vector<std::string>& atom_labels)
{
    // parameters for orbital charge output
    FmtCore fmt_of_chg("%15.4f");
    FmtCore fmt_of_label("%-15s");
    GlobalV::ofs_running << std::endl;
    GlobalV::ofs_running << "-------------------------------------------------------------------------------------------" << std::endl;
    GlobalV::ofs_running << "Orbital Charge Analysis      Charge         Mag(x)         Mag(y)         Mag(z)" << std::endl;
    GlobalV::ofs_running << "-------------------------------------------------------------------------------------------" << std::endl;

    // parameters for mag output
    std::vector<double> mag_x(ucell->nat, 0.0);
    std::vector<double> mag_y(ucell->nat, 0.0);
    std::vector<double> mag_z(ucell->nat, 0.0);
    const std::vector<std::string> title = {"Total Magnetism (uB)", "", "", ""};
    const std::vector<std::string> fmts = {"%-26s", "%20.10f", "%20.10f", "%20.10f"};
    const std::vector<std::string> orb_names = {"s", "p", "d", "f", "g"};
    FmtTable table(/*titles=*/title, 
                   /*nrows=*/ucell->nat, 
                   /*formats=*/fmts, 
                   /*indent=*/0, 
                   /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
    // parameters for mag output
    int occ_index = 0;
    for(int iat=0; iat<ucell->nat; iat++)
    {
        const int it = ucell->iat2it[iat];
        std::string atom_label = atom_labels[it];
        int ia = ucell->iat2ia[iat];
        GlobalV::ofs_running << FmtCore::format("%-20s", atom_label+std::to_string(ia+1)) << std::endl;
        std::vector<double> sum(4, 0.0);
        int current_l = 1;
        std::vector<double> charge_mag(4, 0.0);
        for(int ih=0; ih<iat_nh[iat]; ih++)
        {
            charge_mag[3] += (occs[occ_index] - occs[occ_index + 3]).real();
            charge_mag[1] += (occs[occ_index + 1] + occs[occ_index + 2]).real();
            charge_mag[2] += (occs[occ_index + 1] - occs[occ_index + 2]).imag();
            charge_mag[0] += (occs[occ_index] + occs[occ_index + 3]).real();
            if(ih == current_l * current_l - 1)
            {
                sum[0] += charge_mag[0];
                sum[1] += charge_mag[1];
                sum[2] += charge_mag[2];
                sum[3] += charge_mag[3];
                GlobalV::ofs_running << FmtCore::format("%20s", orb_names[current_l-1])
                    << fmt_of_chg.format(charge_mag[0]) << fmt_of_chg.format(charge_mag[1])
                    << fmt_of_chg.format(charge_mag[2]) << fmt_of_chg.format(charge_mag[3]) << std::endl;
                current_l++;
                charge_mag.assign(4, 0.0);
            }
            occ_index += 4;
        }
        mag_x[iat] = sum[1];
        mag_y[iat] = sum[2];
        mag_z[iat] = sum[3];
        GlobalV::ofs_running << FmtCore::format("%20s", std::string("Sum")) << ""
                    << fmt_of_chg.format(sum[0]) << fmt_of_chg.format(sum[1])
                    << fmt_of_chg.format(sum[2]) << fmt_of_chg.format(sum[3]) << std::endl;
    }
    GlobalV::ofs_running << "-------------------------------------------------------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::endl;
    
    // Print magnetism table
    print_mag_table(atom_labels, mag_x, mag_y, mag_z);
}

void print_mag_table(
    const std::vector<std::string>& atom_labels,
    const std::vector<double>& mag_x,
    const std::vector<double>& mag_y,
    const std::vector<double>& mag_z)
{
    const std::vector<std::string> title = {"Total Magnetism (uB)", "", "", ""};
    const std::vector<std::string> fmts = {"%-26s", "%20.10f", "%20.10f", "%20.10f"};
    FmtTable table(/*titles=*/title, 
                   /*nrows=*/mag_x.size(), 
                   /*formats=*/fmts, 
                   /*indent=*/0, 
                   /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
    
    table << atom_labels << mag_x << mag_y << mag_z;
    GlobalV::ofs_running << table.str() << std::endl;
}

void print_proj_status(int it, int nproj_it)
{
    if(nproj_it == 0)
    {
        std::cout << "BECP_PW >> No projectors defined for type " << it << std::endl;
    }
}

} // namespace print
} // namespace projectors
