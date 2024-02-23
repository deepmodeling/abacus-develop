#include "init_info.h"
#include "module_io/input.h"
#include "../para_json.h"
#include "abacusjson.h"


//Add json objects to init
namespace Json
{

#ifdef __RAPIDJSON


// void gen_init(ModuleSymmetry::Symmetry *symm,Atom *atoms){
//     std::string pgname = symm->pgname;
//     std::string spgname = symm->spgname;
//     AbacusJson::add_json({"init", "pgname"}, pgname,false);
//     AbacusJson::add_json({"init", "spgname"}, spgname,false);


//     int  numAtoms = atoms->na;
//     AbacusJson::add_json({"init", "natom"}, numAtoms,false);
//     AbacusJson::add_json({"init", "nband"}, INPUT.nbands,false);


// }

void gen_init(UnitCell *ucell){
    std::string pgname = ucell->symm.pgname;
    std::string spgname = ucell->symm.spgname;
    AbacusJson::add_json({"init", "point_group"}, pgname,false);
    AbacusJson::add_json({"init", "point_group_in_space"}, spgname,false);


    int  numAtoms = ucell->atoms->na;
    AbacusJson::add_json({"init", "natom"}, numAtoms,false);
    AbacusJson::add_json({"init", "nband"}, GlobalV::NBANDS,false);

    int ntype = ucell->ntype;
    AbacusJson::add_json({"init", "nelectron"}, ntype,false);

    for (int it = 0; it < ntype; it++)
    {
        std::string label = ucell->atoms[it].label;
        int number = ucell->atoms[it].ncpp.zv;
        AbacusJson::add_json({"init", "nelectron_each_type",label}, number,false);
    }    

}


void add_nkstot(int nkstot,int nkstot_ibz){
    Json::AbacusJson::add_json({"init", "nkstot"}, nkstot,false);
    Json::AbacusJson::add_json({"init", "nkstot_ibz"}, nkstot_ibz,false);
}

#endif
} // namespace Json