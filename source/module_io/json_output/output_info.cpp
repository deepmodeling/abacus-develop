#include "output_info.h"
#include "module_io/input.h"
#include "../para_json.h"
#include "abacusjson.h"


//Add json objects to init
namespace Json
{

#ifdef __RAPIDJSON


    // Adjust the position of the json object and set the initial value
    void init_output_array_obj(){


        
        jsonValue scf_obj(JobjectType);

        Json::jsonValue nullValue;
        nullValue.SetNull();
        scf_obj.JaddNormal("e_fermi",nullValue);
        scf_obj.JaddNormal("energy",nullValue);
        scf_obj.JaddNormal("scf_converge",nullValue);

        jsonValue force(JobjectType);
        jsonValue stress(JobjectType);
        jsonValue coordinate(JarrayType);
        jsonValue mag(JarrayType);
        jsonValue cell(JarrayType);

        scf_obj.JaddNormal("force",nullValue);
        scf_obj.JaddNormal("stress",nullValue);
        scf_obj.JaddNormal("coordinate",coordinate);
        scf_obj.JaddNormal("mag",mag);
        scf_obj.JaddNormal("cell",cell);


        AbacusJson::add_json({"output"},scf_obj,true);
    }

    void add_output_cell_coo_stress_force(
        UnitCell *ucell,
        ModuleBase::matrix force, double fac,
        ModuleBase::matrix stress, double unit_transform
    ) {
        int iat = 0;
        const double output_acc = 1.0e-8;

        if (GlobalV::CAL_FORCE){
            //add force
            Json::jsonValue force_array(JarrayType);
            for (int it = 0; it < ucell->ntype; it++)
            {
                for (int ia = 0; ia < ucell->atoms[it].na; ia++)
                {
                    Json::jsonValue force_subarray(JarrayType);
                    double fx = std::abs(force(iat, 0)) > output_acc ? force(iat, 0) * fac : 0.0;
                    double fy = std::abs(force(iat, 1)) > output_acc ? force(iat, 1) * fac : 0.0;
                    double fz = std::abs(force(iat, 2)) > output_acc ? force(iat, 2) * fac : 0.0;

                    force_subarray.JPushBack(fx);
                    force_subarray.JPushBack(fy);
                    force_subarray.JPushBack(fz);
                    force_array.JPushBack(force_subarray);
                    iat++;
                }
            }
            Json::AbacusJson::add_json({"output",-1,"force"}, force_array,false);

            // AbacusJson::add_Json(force_array,false,"output",-1,"force");
        }

        if (GlobalV::CAL_STRESS){
        //add stress
            Json::jsonValue stress_array(JarrayType);
            for (int i = 0; i < 3; i++)
            {
                Json::jsonValue stress_subarray(JarrayType);
                double sx = stress(i, 0) * unit_transform;
                double sy = stress(i, 1) * unit_transform;
                double sz = stress(i, 2) * unit_transform;
                stress_subarray.JPushBack(sx);
                stress_subarray.JPushBack(sy);
                stress_subarray.JPushBack(sz);
                stress_array.JPushBack(stress_subarray);
            }            
            Json::AbacusJson::add_json({"output",-1,"stress"}, stress_array,false);

            // AbacusJson::add_Json(stress_array,false,"output",-1,"stress");
        }
        //add coordinate
        int ntype = ucell->ntype;
        double lat0 = ucell->lat0;
        for(int i=0;i<ntype;i++){
            ModuleBase::Vector3<double>* tau = ucell->atoms[i].tau;
            int na = ucell->atoms[i].na;
            for(int j=0;j<na;j++){
                Json::jsonValue coordinateArray(JarrayType);
                coordinateArray.JPushBack(tau[j][0]*lat0);
                coordinateArray.JPushBack(tau[j][1]*lat0);
                coordinateArray.JPushBack(tau[j][2]*lat0);
                Json::AbacusJson::add_json({"output",-1,"coordinate"}, coordinateArray,true);
                Json::AbacusJson::add_json( {"output",-1,"mag"},ucell->atoms[i].mag[j],true);
            }
        }        

        //add cell 
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
            Json::AbacusJson::add_json({"output",-1,"cell"}, cellArray1,true);
            Json::AbacusJson::add_json({"output",-1,"cell"}, cellArray2,true);
            Json::AbacusJson::add_json({"output",-1,"cell"}, cellArray3,true);

            // Json::AbacusJson::add_Json(cellArray1,true,"output",-1,"cell");
            // Json::AbacusJson::add_Json(cellArray2,true,"output",-1,"cell");
            // Json::AbacusJson::add_Json(cellArray2,true,"output",-1,"cell");
        }

    }

    void add_output_efermi_energy_converge(double efermi, double energy ,bool scf_converge ){
        Json::AbacusJson::add_json({"output",-1,"e_fermi"}, efermi,false);
        Json::AbacusJson::add_json({"output",-1,"energy"}, energy,false);
        Json::AbacusJson::add_json({"output",-1,"scf_converge"}, scf_converge,false);
    
        // Json::AbacusJson::add_Json(efermi,false,"output",-1,"e_fermi");
        // Json::AbacusJson::add_Json(energy,false,"output",-1,"energy");
        // Json::AbacusJson::add_Json(scf_converge,false,"output",-1,"scf_converge");
    
    }

    void add_output_scf_mag(
        double total_mag, double absolute_mag,
        double energy, double ediff, double drho,double time
    ){
        Json::AbacusJson::add_json({"output",-1,"total_mag"}, total_mag,false);
        Json::AbacusJson::add_json({"output",-1,"absolute_mag"}, absolute_mag,false);

        // Json::AbacusJson::add_Json(total_mag,false,"output",-1,"total_mag");
        // Json::AbacusJson::add_Json(absolute_mag,false,"output",-1,"absolute_mag");

        Json::jsonValue scf_obj(JobjectType);
        scf_obj.JaddNormal("energy",energy);
        scf_obj.JaddNormal("ediff",ediff);    
        scf_obj.JaddNormal("drho",drho);    
        scf_obj.JaddNormal("time",time);          
        Json::AbacusJson::add_json({"output",-1,"scf"}, scf_obj,true);

        // Json::AbacusJson::add_Json(scf_obj,true,"output",-1,"scf");
    }

#endif
} // namespace Json