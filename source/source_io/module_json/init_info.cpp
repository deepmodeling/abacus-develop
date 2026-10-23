#include "init_info.h"

#include "abacusjson.h"
#include "source_cell/atom_spec.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"

#ifdef __JSON
#include <nlohmann/json.hpp>
#include <stdexcept>

namespace Json
{
namespace
{
jsonValue& init_section()
{
    jsonValue& init = *AbacusJson::document().emplace("init", jsonValue::object()).first;
    if (!init.is_object())
    {
        throw std::invalid_argument("JSON init section must be an object");
    }
    return init;
}
} // namespace

void gen_init(UnitCell* ucell, const Input_para& inp)
{
    jsonValue info = {{"point_group", ucell->symm.pgname},
                      {"point_group_in_space", ucell->symm.spgname},
                      {"natom", ucell->nat},
                      {"nband", inp.nbands}};

    int nelec_total = 0;
    for (int it = 0; it < ucell->ntype; ++it)
    {
        const Atom& atom = ucell->atoms[it];
        nelec_total += atom.ncpp.zv * atom.na;
        info["natom_each_type"][atom.label] = atom.na;
        info["nelectron_each_type"][atom.label] = atom.ncpp.zv;
    }

    info["nelectron"] = nelec_total;
    info["ecutwfc"] = inp.ecutwfc;
    info["ecutwfc_unit"] = "Ry";
    info["smearing_method"] = inp.smearing_method;
    info["smearing_sigma"] = inp.smearing_sigma;
    info["smearing_sigma_unit"] = "Ry";
    info["kmesh_type"] = inp.kmesh_type;
    info["kspacing"] = jsonValue::array({inp.kspacing[0], inp.kspacing[1], inp.kspacing[2]});
    info["koffset"] = jsonValue::array({inp.koffset[0], inp.koffset[1], inp.koffset[2]});
    // Shallow update: preserve other generators' fields, replace this generator's containers.
    init_section().update(info);
}

void add_nkstot(int nkstot)
{
    init_section()["nkstot"] = nkstot;
}

void gen_stru(UnitCell* ucell, const Input_para& inp)
{
    AbacusJson::document()["comment"] =
        "Unless otherwise specified, the unit of energy is eV and the unit of length is Angstrom";

    jsonValue info = jsonValue::object();
    for (int it = 0; it < ucell->ntype; ++it)
    {
        const Atom& atom = ucell->atoms[it];
        info["element"][atom.label] = atom.ncpp.psd;
        const std::string orbital = inp.orbital_dir + ucell->orbital_fn[it];
        info["orb"][atom.label] = orbital.empty() ? jsonValue(nullptr) : jsonValue(orbital);
        info["pp"][atom.label] = ucell->pseudo_fn[it];
    }

    const double lat0_angstrom = ucell->lat0_angstrom;
    for (int it = 0; it < ucell->ntype; ++it)
    {
        const Atom& atom = ucell->atoms[it];
        for (int ia = 0; ia < atom.na; ++ia)
        {
            const ModuleBase::Vector3<double>& tau = atom.tau[ia];
            info["coordinate"].push_back(jsonValue::array({tau[0] * lat0_angstrom,
                                                         tau[1] * lat0_angstrom,
                                                         tau[2] * lat0_angstrom}));
            info["mag"].push_back(atom.mag[ia]);
            info["label"].push_back(atom.label);
        }
    }
    info["cell"] = {{ucell->latvec.e11 * lat0_angstrom,
                     ucell->latvec.e12 * lat0_angstrom,
                     ucell->latvec.e13 * lat0_angstrom},
                    {ucell->latvec.e21 * lat0_angstrom,
                     ucell->latvec.e22 * lat0_angstrom,
                     ucell->latvec.e23 * lat0_angstrom},
                    {ucell->latvec.e31 * lat0_angstrom,
                     ucell->latvec.e32 * lat0_angstrom,
                     ucell->latvec.e33 * lat0_angstrom}};
    init_section().update(info);
}

} // namespace Json
#endif // __JSON
