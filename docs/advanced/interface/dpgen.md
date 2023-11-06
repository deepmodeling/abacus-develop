# DP-GEN


[DP-GEN](https://github.com/deepmodeling/dpgen), the deep potential generator, is a package designed to generate deep learning based model of interatomic potential energy and force fields (Yuzhi Zhang, Haidi Wang, Weijie Chen, Jinzhe Zeng, Linfeng Zhang, Han Wang, and Weinan E, DP-GEN: A concurrent learning platform for the generation of reliable deep learning based potential energy models, Computer Physics Communications, 2020, 107206). ABACUS can now interface with DP-GEN to generate deep potentials and perform autotests. The minimum recommended version is ABACUS 3.0, dpdata 0.2.8, and dpgen 0.10.7 . In the following part, we take the FCC aluminum as an example.

## init_bulk and run

This example can be found in examples/dpgen-example/init_and_run directory.

Firstly, one needs to prepare input files for ABACUS calculation, e.g., “INPUT”, "INPUT.md", "KPT", "Al.STRU", "Al_ONCV_PBE-1.0.upf", which are the main input file containing input tags, k-point mesh, crystal structure and pseudoptential, respectively. "INPUT" is for scf calculation, and "INPUT.md" is for AIMD (ab-initio molecular dynamic) calculation.

Secondly, for the "dpgen init_bulk" step, an `init.json` file should be provided:


```json
{
    "init_fp_style":    "ABACUS",   # abacus interface
    "stages":           [1,2,3,4],
    "cell_type":        "fcc",
    "super_cell":       [2, 1, 1],
    "elements":         ["Al"],
    "from_poscar":      true,
    "from_poscar_path": "./Al.STRU",
    "potcars":          ["Al_ONCV_PBE-1.0.upf"],
    "relax_incar":      "INPUT",
    "relax_kpt":        "KPT",
    "md_incar" :        "INPUT.md",
    "md_kpt" :          "KPT",
    "skip_relax":       false,
    "scale":            [1.00],
    "pert_numb":        10,
    "pert_box":         0.01,
    "pert_atom":        0.01,
    "coll_ndata":       10,
    "_comment":         "that's all"
}
```

Next, for the "dpgen run" step, the following `run_param.json` should be provided.
```json
{
    "type_map": [
        "Al"
    ],
    "mass_map": [
        26.9815
    ],
    "init_data_prefix": "./",
    "init_data_sys": [
        "Al.STRU.01x01x01/02.md/sys-0004/deepmd"
    ],
    "sys_format": "abacus/stru",   # the initial structures are in ABACUS/STRU formate
    "sys_configs_prefix": "./",
    "sys_configs": [
        [
            "Al.STRU.01x01x01/01.scale_pert/sys-0004/scale*/00000*/STRU" 
        ],
        [
            "Al.STRU.01x01x01/01.scale_pert/sys-0004/scale*/000010/STRU" 
        ]
    ],
    "_comment": " 00.train ",
    "numb_models": 4,
    "default_training_param": {
        "model": {
            "type_map": [
                "Al"
            ],
            "descriptor": {
                "type": "se_e2_a",
                "sel": [
                    16
                ],
                "rcut_smth": 0.5,
                "rcut": 5.0,
                "neuron": [
                    10,
                    20,
                    40
                ],
                "resnet_dt": true,
                "axis_neuron": 12,
                "seed": 1
            },
            "fitting_net": {
                "neuron": [
                    25,
                    50,
                    100
                ],
                "resnet_dt": false,
                "seed": 1
            }
        },
        "learning_rate": {
            "type": "exp",
            "start_lr": 0.001,
            "decay_steps": 100
        },
        "loss": {
            "start_pref_e": 0.02,
            "limit_pref_e": 2,
            "start_pref_f": 1000,
            "limit_pref_f": 1,
            "start_pref_v": 0.0,
            "limit_pref_v": 0.0
        },
        "training": {
            "stop_batch": 20000,
            "disp_file": "lcurve.out",
            "disp_freq": 1000,
            "numb_test": 4,
            "save_freq": 1000,
            "save_ckpt": "model.ckpt",
            "disp_training": true,
            "time_training": true,
            "profiling": false,
            "profiling_file": "timeline.json",
            "_comment": "that's all"
        }
    },
    "_comment": "01.model_devi ",
    "model_devi_dt": 0.002,
    "model_devi_skip": 0,
    "model_devi_f_trust_lo": 0.05,
    "model_devi_f_trust_hi": 0.15,
    "model_devi_clean_traj": false,
    "model_devi_jobs": [
        {
            "sys_idx": [
                0
            ],
            "temps": [
                50
            ],
            "press": [
                1.0
            ],
            "trj_freq": 10,
            "nsteps": 300,
            "ensemble": "nvt",
            "_idx": "00"
        },
        {
            "sys_idx": [
                1
            ],
            "temps": [
                50
            ],
            "press": [
                1.0
            ],
            "trj_freq": 10,
            "nsteps": 3000,
            "ensemble": "nvt",
            "_idx": "01"
        }
    ],
    "fp_style": "abacus/scf",
    "shuffle_poscar": false,
    "fp_task_max": 20,
    "fp_task_min": 5,
    "fp_pp_path": "./",
    "fp_pp_files": ["Al_ONCV_PBE-1.0.upf"],   # the pseudopotential file
    "fp_orb_files": ["Al_gga_9au_100Ry_4s4p1d.orb"],  # the orbital file (use only in LCAO calculation)
    "k_points":[2, 2, 2, 0, 0, 0],  # k-mesh setting
    "user_fp_params":{  # All the ABACUS input paramters are defined here
    "ntype": 1,         # defining input parameters from INPUT files is not supported yet.
    "ecutwfc": 80,      
    "mixing_type": "broyden",
    "mixing_beta": 0.8,
    "symmetry": 1,
    "nspin": 1,
    "ks_solver": "cg",
    "smearing_method": "mp",
    "smearing_sigma": 0.002,
    "scf_thr":1e-8,
    "cal_force":1,        # calculate force must be set to 1 in dpgen calculation
    "kspacing": 0.01  # when KSPACING is set, the above k_points setting becomes invalid.
    }
}
```

## autotest

This example can be found in examples/dpgen-example/autotest directory.

`dpgen autotest` supports to perform `relaxation`,`eos` (equation of state),`elastic`,`surface`,`vacancy`, and `interstitial` calculations with ABACUS. A `property.json` and `machine.json` file need to be provided. For example,

`property.json`:
```json

{
    "structures":    ["confs/"],
    "interaction": {
        "type":         "abacus",
        "incar":        "./INPUT",
        "potcar_prefix":"./",
        "potcars":      {"Al": "Al.PD04.PBE.UPF"},
        "orb_files": {"Al":"Al_gga_10au_100Ry_3s3p2d.orb"}
    },
    "_relaxation": {
            "cal_type": "relaxation",
            "cal_setting":{
                    "input_prop": "./INPUT.rlx"
            }
     },
    "properties": [
        {
         "type":         "eos",
         "vol_start":    0.85,
         "vol_end":      1.15,
         "vol_step":     0.01,
         "cal_setting": {
                         "relax_pos": true,
                         "relax_shape": true,
                         "relax_vol": false,
                         "overwrite_interaction":{
                                     "type": "abacus",
                                     "incar": "./INPUT",
                                     "potcar_prefix":"./",
                                     "orb_files": {"Al":"Al_gga_10au_100Ry_3s3p2d.orb"},
                                     "potcars": {"Al": "Al.PD04.PBE.UPF"} }
                        }
        },
         {
         "type":         "elastic",
         "skip":         false,
         "norm_deform":   1e-2,
         "shear_deform":  1e-2
        },
        {
         "type":         "vacancy",
         "skip":         false,
         "supercell":    [2, 2, 2]
        },
        {
         "type":           "surface",
         "skip":         true,
         "min_slab_size":  15,
         "min_vacuum_size":11,
         "pert_xz":        0.01,
         "max_miller":     3,
         "cal_type":       "static"
        }
        ]
}
```

`machine.json`

```json
{
  "api_version": "1.0",
  "deepmd_version": "2.1.0",
  "train" :[
    {
      "command": "dp",
      "machine": {
        "batch_type": "DpCloudServer",
        "context_type": "DpCloudServerContext",
        "local_root" : "./",
        "remote_profile":{
          "email": "xxx@xxx.xxx",
          "password": "xxx",
          "program_id": 000,
            "input_data":{
                "api_version":2,
                "job_type": "indicate",
                "log_file": "00*/train.log",
                "grouped":true,
                "job_name": "Al-train-VASP",
                "disk_size": 100,
                "scass_type":"c8_m32_1 * NVIDIA V100",
                "platform": "ali",
                "image_name":"LBG_DeePMD-kit_2.1.0_v1",
                "on_demand":0
            }
        }
      },
      "resources": {
        "number_node":123473334635,
        "local_root":"./",
        "cpu_per_node": 4,
        "gpu_per_node": 1,
        "queue_name": "GPU",
        "group_size": 1
      }
    }],
  "model_devi":
    [{
      "command": "lmp -i input.lammps -v restart 0",
      "machine": {
        "batch_type": "DpCloudServer",
        "context_type": "DpCloudServerContext",
        "local_root" : "./",
        "remote_profile":{
          "email": "xxx@xxx.xxx",
          "password": "xxx",
          "program_id": 000,
            "input_data":{
              "api_version":2,
              "job_type": "indicate",
              "log_file": "*/model_devi.log",
              "grouped":true,
              "job_name": "Al-devia-ABACUS",
              "disk_size": 200,
              "scass_type":"c8_m32_1 * NVIDIA V100",
              "platform": "ali",
              "image_name":"LBG_DeePMD-kit_2.1.0_v1",
              "on_demand":0
            }
        }
      },
      "resources": {
        "number_node": 28348383,
        "local_root":"./",
        "cpu_per_node": 4,
        "gpu_per_node": 1,
        "queue_name": "GPU",
        "group_size": 100
      }
    }],
  "fp":
    [{
      "command": "OMP_NUM_THREADS=1 mpirun -np 16 abacus",
      "machine": {
        "batch_type": "DpCloudServer",
        "context_type": "DpCloudServerContext",
        "local_root" : "./",
        "remote_profile":{
          "email": "xxx@xxx.xxx",
          "password": "xxx",
         "program_id": 000,
            "input_data":{
              "api_version":2,
              "job_type": "indicate",
              "log_file": "task*/fp.log",
              "grouped":true,
              "job_name": "al-DFT-test",
              "disk_size": 100,
              "scass_type":"c32_m128_cpu",
              "platform": "ali",
              "image_name":"XXXXX",
              "on_demand":0
            }
        }
      },
      "resources": {
        "number_node": 712254638889,
        "cpu_per_node": 32,
        "gpu_per_node": 0,
        "queue_name": "CPU",
        "group_size": 2,
        "local_root":"./",
        "source_list": ["/opt/intel/oneapi/setvars.sh"]
      }
    }
  ]
}

```

For each property, the command `dpgen autotest make property.json` will generate the input files, `dpgen autotest run property.json machine.json` will run the corresponding tasks, and `dpgen autotest post property.json` will collect the final results. 

Notes:
- The ABACUS-DPGEN interface can be used in both pw and lcao basis.

