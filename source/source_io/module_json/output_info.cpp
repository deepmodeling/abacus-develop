#include "output_info.h"

#include "abacusjson.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

#ifdef __JSON
#include <nlohmann/json.hpp>
#endif

#include <cmath>
#include <utility>

namespace Json
{

#ifdef __JSON

void init_output_array_obj()
{
    AbacusJson::append_json({"output"},
                            {{"e_fermi", nullptr},
                             {"energy", nullptr},
                             {"scf_converge", nullptr},
                             {"force", nullptr},
                             {"stress", nullptr},
                             {"coordinate", jsonValue::array()},
                             {"mag", jsonValue::array()},
                             {"cell", jsonValue::array()}});
}

void add_output_cell_coo_stress_force(const UnitCell& ucell,
                                      const ModuleBase::matrix& force,
                                      const double fac,
                                      const ModuleBase::matrix& stress,
                                      const double unit_transform,
                                      const bool cal_force,
                                      const bool cal_stress)
{
    const double output_acc = 1.0e-8;
    if (cal_force)
    {
        jsonValue force_array = jsonValue::array();
        int iat = 0;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
            {
                const double fx = std::abs(force(iat, 0)) > output_acc ? force(iat, 0) * fac : 0.0;
                const double fy = std::abs(force(iat, 1)) > output_acc ? force(iat, 1) * fac : 0.0;
                const double fz = std::abs(force(iat, 2)) > output_acc ? force(iat, 2) * fac : 0.0;
                force_array.push_back(jsonValue::array({fx, fy, fz}));
                ++iat;
            }
        }
        AbacusJson::set_json({"output", -1, "force"}, std::move(force_array));
    }

    if (cal_stress)
    {
        jsonValue stress_array = jsonValue::array();
        for (int i = 0; i < 3; ++i)
        {
            stress_array.push_back(jsonValue::array({stress(i, 0) * unit_transform,
                                                     stress(i, 1) * unit_transform,
                                                     stress(i, 2) * unit_transform}));
        }
        AbacusJson::set_json({"output", -1, "stress"}, std::move(stress_array));
    }

    const double lat0_angstrom = ucell.lat0_angstrom;
    jsonValue coordinates = jsonValue::array();
    jsonValue mag = jsonValue::array();
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            const ModuleBase::Vector3<double>& tau = ucell.atoms[it].tau[ia];
            coordinates.push_back(jsonValue::array({tau[0] * lat0_angstrom,
                                                    tau[1] * lat0_angstrom,
                                                    tau[2] * lat0_angstrom}));
            mag.push_back(ucell.atoms[it].mag[ia]);
        }
    }
    AbacusJson::set_json({"output", -1, "coordinate"}, std::move(coordinates));
    AbacusJson::set_json({"output", -1, "mag"}, std::move(mag));
    AbacusJson::set_json({"output", -1, "cell"},
                         {{ucell.latvec.e11 * lat0_angstrom,
                           ucell.latvec.e12 * lat0_angstrom,
                           ucell.latvec.e13 * lat0_angstrom},
                          {ucell.latvec.e21 * lat0_angstrom,
                           ucell.latvec.e22 * lat0_angstrom,
                           ucell.latvec.e23 * lat0_angstrom},
                          {ucell.latvec.e31 * lat0_angstrom,
                           ucell.latvec.e32 * lat0_angstrom,
                           ucell.latvec.e33 * lat0_angstrom}});
}

void add_output_efermi_converge(const double efermi, const bool scf_converge)
{
    AbacusJson::set_json({"output", -1, "e_fermi"}, efermi);
    AbacusJson::set_json({"output", -1, "scf_converge"}, scf_converge);
}

void add_output_energy(const double energy)
{
    AbacusJson::set_json({"output", -1, "energy"}, energy);
}

void add_output_scf_mag(const double total_mag,
                        const double absolute_mag,
                        const double energy,
                        const double ediff,
                        const double drho,
                        const double time)
{
    AbacusJson::set_json({"output", -1, "total_mag"}, total_mag);
    AbacusJson::set_json({"output", -1, "absolute_mag"}, absolute_mag);
    AbacusJson::append_json({"output", -1, "scf"},
                            {{"energy", energy}, {"ediff", ediff}, {"drho", drho}, {"time", time}});
}

#endif // __JSON
} // namespace Json
