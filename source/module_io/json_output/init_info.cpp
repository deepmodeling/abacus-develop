#include "init_info.h"
#include "module_io/input.h"
#include "../para_json.h"
#include "abacusjson.h"


//Add json objects to init
namespace Json
{

#ifdef __RAPIDJSON



void gen_init(UnitCell *ucell){
    std::string pgname = ucell->symm.pgname;
    std::string spgname = ucell->symm.spgname;
    AbacusJson::add_json({"init", "point_group"}, pgname,false);
    AbacusJson::add_json({"init", "point_group_in_space"}, spgname,false);


    int  numAtoms = ucell->atoms->na;
    AbacusJson::add_json({"init", "natom"}, numAtoms,false);
    AbacusJson::add_json({"init", "nband"}, GlobalV::NBANDS,false);

    int ntype = ucell->ntype,nelec_total=0;
    AbacusJson::add_json({"init", "nelectron"}, 0,false);
    for (int it = 0; it < ntype; it++)
    {
        std::string label = ucell->atoms[it].label;
        int number = ucell->atoms[it].ncpp.zv;

        nelec_total+=ucell->atoms[it].ncpp.zv * ucell->atoms[it].na;
        AbacusJson::add_json({"init", "nelectron_each_type",label}, number,false);
    }    
    
    AbacusJson::add_json({"init", "nelectron"}, nelec_total,false);


}


void add_nkstot(int nkstot,int nkstot_ibz){
    Json::AbacusJson::add_json({"init", "nkstot"}, nkstot,false);
    Json::AbacusJson::add_json({"init", "nkstot_ibz"}, nkstot_ibz,false);
}


void gen_stru(UnitCell *ucell){

    int ntype = ucell->ntype;


    
    //array of pseudopotential file
    std::string* pseudo_fn = ucell->pseudo_fn;

    //array of orbital file
    std::string* orbital_fn = ucell->orbital_fn;


    //add atom element,orbital file and pseudopotential file
    for(int i=0;i<ntype;i++){
        std::string atom_label = ucell->atoms[i].label;

        std::string atom_element = ucell->atoms[i].ncpp.psd;

        Json::jsonValue element_obj(JobjectType);
        element_obj.JaddStringKV(atom_label,atom_element);
        Json::AbacusJson::add_json({"init","element"}, element_obj,false);
    

        std::string orbital_str = GlobalV::global_orbital_dir + orbital_fn[i];
        if(!orbital_str.compare("")){
            Json::jsonValue nullValue;
            nullValue.SetNull();
            Json::AbacusJson::add_json({"init","orb",atom_label}, nullValue,false);
        }else {
            Json::AbacusJson::add_json({"init","orb",atom_label}, orbital_str,false);
        }    
        std::string pseudo_str = pseudo_fn[i];
        Json::AbacusJson::add_json({"init","pp",atom_label}, pseudo_str,false);
    
 
    }
    

    //atom coordinate, mag and label
    double lat0 = ucell->lat0;
    std::string* label = ucell->atom_label;
    for(int i=0;i<ntype;i++){
        ModuleBase::Vector3<double>* tau = ucell->atoms[i].tau;
        int na = ucell->atoms[i].na;
        for(int j=0;j<na;j++){
            Json::jsonValue coordinateArray(JarrayType);
            coordinateArray.JPushBack(tau[j][0]*lat0);
            coordinateArray.JPushBack(tau[j][1]*lat0);
            coordinateArray.JPushBack(tau[j][2]*lat0);
            Json::AbacusJson::add_json({"init","coordinate"}, coordinateArray,true);
        
            Json::AbacusJson::add_json({"init","mag"}, ucell->atoms[i].mag[j],true);

    std::string str = label[i];
    Json::AbacusJson::add_json({"init","label"}, str,true);        
        }

    }


    //cell 
    {
        Json::jsonValue cellArray1(JarrayType);
        Json::jsonValue cellArray2(JarrayType);
        Json::jsonValue cellArray3(JarrayType);
        cellArray1.JPushBack(ucell->latvec.e11);
        cellArray1.JPushBack(ucell->latvec.e12);
        cellArray1.JPushBack(ucell->latvec.e13);
        cellArray2.JPushBack(ucell->latvec.e21);
        cellArray2.JPushBack(ucell->latvec.e22);
        cellArray2.JPushBack(ucell->latvec.e23);
        cellArray3.JPushBack(ucell->latvec.e31);
        cellArray3.JPushBack(ucell->latvec.e32);
        cellArray3.JPushBack(ucell->latvec.e33);
        Json::AbacusJson::add_json({"init","cell"}, cellArray1,true);
        Json::AbacusJson::add_json({"init","cell"}, cellArray2,true);
        Json::AbacusJson::add_json({"init","cell"}, cellArray3,true);
    }
    return;
}


#endif
} // namespace Json