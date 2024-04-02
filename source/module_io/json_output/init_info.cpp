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

    // Json::AbacusJson::add_Json(pgname,false,"init", "point_group");
    // Json::AbacusJson::add_Json(spgname,false,"init","point_group_in_space");

    int  numAtoms = ucell->nat;
    AbacusJson::add_json({"init", "natom"}, numAtoms,false);
    AbacusJson::add_json({"init", "nband"}, GlobalV::NBANDS,false);

    // Json::AbacusJson::add_Json(numAtoms,false,"init", "natom");
    // Json::AbacusJson::add_Json(GlobalV::NBANDS,false,"init", "nband");

    int ntype = ucell->ntype,nelec_total=0;
    for (int it = 0; it < ntype; it++)
    {
        std::string label = ucell->atoms[it].label;
        int atom_number = ucell->atoms[it].na;
        int number = ucell->atoms[it].ncpp.zv;

        nelec_total+=ucell->atoms[it].ncpp.zv * ucell->atoms[it].na;
        AbacusJson::add_json({"init", "natom_each_type",label}, atom_number,false);
        AbacusJson::add_json({"init", "nelectron_each_type",label}, number,false);


        //Json::AbacusJson::add_Json(number,false,"init", "nelectron_each_type",label);
    }    
    
    AbacusJson::add_json({"init", "nelectron"}, nelec_total,false);

    // Json::AbacusJson::add_Json(nelec_total,false,"init", "nelectron");
}


void add_nkstot(int nkstot,int nkstot_ibz){
    Json::AbacusJson::add_json({"init", "nkstot"}, nkstot,false);
    Json::AbacusJson::add_json({"init", "nkstot_ibz"}, nkstot_ibz,false);

    // Json::AbacusJson::add_Json(nkstot,false,"init", "nkstot");
    // Json::AbacusJson::add_Json(nkstot_ibz,false,"init", "nkstot_ibz");
}


void gen_stru(UnitCell *ucell){
    AbacusJson::add_json({"comment"},"Unless otherwise specified, the unit of energy is eV and the unit of length is Angstrom",false);

    int ntype = ucell->ntype;


    
    //array of pseudopotential file
    std::string* pseudo_fn = ucell->pseudo_fn;

    //array of orbital file
    std::string* orbital_fn = ucell->orbital_fn;


    //add atom element,orbital file and pseudopotential file
    for(int i=0;i<ntype;i++){
        std::string atom_label = ucell->atoms[i].label;

        std::string atom_element = ucell->atoms[i].ncpp.psd;

        Json::AbacusJson::add_json({"init","element",atom_label}, atom_element,false);


        std::string orbital_str = GlobalV::global_orbital_dir + orbital_fn[i];
        if(!orbital_str.compare("")){
            Json::jsonValue nullValue;
            nullValue.SetNull();
            Json::AbacusJson::add_json({"init","orb",atom_label}, nullValue,false);
        
            // Json::AbacusJson::add_Json(nullValue,false,"init","orb",atom_label);

        }else {
            Json::AbacusJson::add_json({"init","orb",atom_label}, orbital_str,false);
            // Json::AbacusJson::add_Json(orbital_str,false,"init","orb",atom_label);
        }    
        std::string pseudo_str = pseudo_fn[i];
        Json::AbacusJson::add_json({"init","pp",atom_label}, pseudo_str,false);

        // Json::AbacusJson::add_Json(pseudo_str,false,"init","pp",atom_label);
 
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
            // Json::AbacusJson::add_Json(coordinateArray,true,"init","coordinate");
 
            Json::AbacusJson::add_json({"init","mag"}, ucell->atoms[i].mag[j],true);

            // Json::AbacusJson::add_Json(ucell->atoms[i].mag[j],true,"init","mag");
 
            std::string str = label[i];
            Json::AbacusJson::add_json({"init","label"}, str,true);  
            // Json::AbacusJson::add_Json(str,true,"init","label");      
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

        // Json::AbacusJson::add_Json(cellArray1,true,"init","cell"); 
        // Json::AbacusJson::add_Json(cellArray2,true,"init","cell");  
        // Json::AbacusJson::add_Json(cellArray3,true,"init","cell");   
    }
    return;
}


#endif
} // namespace Json