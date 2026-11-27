#include "output_info.h"

#include "abacusjson.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

#ifdef __JSON
#include <nlohmann/json.hpp>
#include <stdexcept>
#endif

#include <cmath>
#include <utility>

namespace Json
{

#ifdef __JSON

namespace
{
jsonValue& current_output()
{
    jsonValue& root = AbacusJson::document();
    const jsonValue::iterator output = root.find("output");
    if (output == root.end() || !output->is_array())
    {
        throw std::invalid_argument("JSON output records must be initialized as an array");
    }
    if (output->empty())
    {
        throw std::out_of_range("JSON output record is not initialized");
    }
    jsonValue& record = output->back();
    if (!record.is_object())
    {
        throw std::invalid_argument("JSON output record must be an object");
    }
    return record;
}
} // namespace

void init_output_array_obj()
{
    jsonValue& output = *AbacusJson::document().emplace("output", jsonValue::array()).first;
    if (!output.is_array())
    {
        throw std::invalid_argument("JSON output must be an array");
    }
    output.push_back({{"e_fermi", nullptr},
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
    jsonValue& output = current_output();
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
        output["force"] = std::move(force_array);
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
        output["stress"] = std::move(stress_array);
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
    output["coordinate"] = std::move(coordinates);
    output["mag"] = std::move(mag);
    output["cell"] = {{ucell.latvec.e11 * lat0_angstrom,
                       ucell.latvec.e12 * lat0_angstrom,
                       ucell.latvec.e13 * lat0_angstrom},
                      {ucell.latvec.e21 * lat0_angstrom,
                       ucell.latvec.e22 * lat0_angstrom,
                       ucell.latvec.e23 * lat0_angstrom},
                      {ucell.latvec.e31 * lat0_angstrom,
                       ucell.latvec.e32 * lat0_angstrom,
                       ucell.latvec.e33 * lat0_angstrom}};
}

void add_output_efermi_converge(const double efermi, const bool scf_converge)
{
    jsonValue& output = current_output();
    output["e_fermi"] = efermi;
    output["scf_converge"] = scf_converge;
}

void add_output_energy(const double energy)
{
    current_output()["energy"] = energy;
}

void add_output_scf_mag(const double total_mag,
                        const double absolute_mag,
                        const double energy,
                        const double ediff,
                        const double drho,
                        const double time)
{
    jsonValue& output = current_output();
    output["total_mag"] = total_mag;
    output["absolute_mag"] = absolute_mag;
    // Acquire the history only after inserting other fields: ordered_json may reallocate them.
    jsonValue& scf = *output.emplace("scf", jsonValue::array()).first;
    if (!scf.is_array())
    {
        throw std::invalid_argument("JSON SCF history must be an array");
    }
    scf.push_back({{"energy", energy}, {"ediff", ediff}, {"drho", drho}, {"time", time}});
}

#endif // __JSON
} // namespace Json
