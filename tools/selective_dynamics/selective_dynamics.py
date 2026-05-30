#!/usr/bin/env python3
import os
import re
import sys
import glob
import time
import yaml
import shutil
import subprocess
from ase import Atoms
from ase.io import read, write


def split_stru_with_index(input_file, selected_indices, output1, output2):

    atoms = read(input_file, verbose=True)
    pp = atoms.info["pp"]
    basis = atoms.info["basis"]

    for i, atom in enumerate(atoms):
        atom.tag = i  

    atoms1 = atoms[selected_indices]


    remaining_indices = [i for i in range(len(atoms)) if i not in selected_indices]
    atoms2 = atoms[remaining_indices]

    
    write(output1, atoms1, format='abacus', pp=pp, basis=basis)
    write(output2, atoms2, format='abacus', pp=pp, basis=basis)
    
    print(f" Success! \n A.stru: {len(atoms1)} atoms\n B.stru: {len(atoms2)} atoms\n")
    return selected_indices, remaining_indices, pp, basis


def merge_stru_by_index(file1, file2, indices1, indices2, pp, basis, output_file):

    atoms1 = read(file1)
    atoms2 = read(file2)
    

    total_atoms = len(atoms1) + len(atoms2)
    merged_atoms = [None] * total_atoms
    
    for i, idx in enumerate(indices1):
        merged_atoms[idx] = atoms1[i]
    
    for i, idx in enumerate(indices2):
        merged_atoms[idx] = atoms2[i]
    

    cell = atoms2.get_cell()
    pbc = atoms2.get_pbc()

    positions = [atom.position for atom in merged_atoms]
    symbols = [atom.symbol for atom in merged_atoms]
    
    merged = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=pbc)

    write(output_file, merged, format='abacus', pp=pp, basis=basis)
    print(f"success! total {len(merged)} atoms")


def parse_forces_from_file(indices, file_path="running_scf.log"):

    with open(file_path, 'r', encoding='utf-8') as f:
        contents = f.read()

    lines = contents.split('\n')

    start_index = -1
    for i, line in enumerate(lines):
        if "#TOTAL-FORCE (eV/Angstrom)#" in line:
            start_index = i+4
            break

    results = []
    for i in indices:
        results.append(lines[start_index+i])

    result_text = '\n'.join(results)
    return result_text


def load_config(config_file="config.yaml"):

    with open(config_file, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def submit_jobs():

    # loading settings
    config = load_config()
    
    # split structure into A and B
    origin_structure = config['origin_structure']
    selected_indices = config['selected_indices']
    indices1, indices2, pp, basis = split_stru_with_index(origin_structure, selected_indices, 
                                                'A.stru', 'B.stru')


    # setting.conf for phonopy
    setting_conf_content = config['setting.conf']
    with open("setting.conf", "w", encoding="utf-8") as f:
        f.write(setting_conf_content)


    # generate perturbed structures
    result = subprocess.run("phonopy setting.conf --abacus -d -c A.stru", 
                            shell=True, capture_output=True, text=True)
    print(result.stdout)


    all_files = sorted(glob.glob("STRU-*"))
    TOTAL_TASKS = len(all_files)
    TASKS_PER_BATCH = config['tasks_per_batch']
    WAIT_TIME = config['wait_time']
    TOTAL_BATCHES = (TOTAL_TASKS + TASKS_PER_BATCH - 1) // TASKS_PER_BATCH  
    print(f' TOTAL_TASKS = {TOTAL_TASKS}\n TASKS_PER_BATCH = {TASKS_PER_BATCH}\n TOTAL_BATCHES = {TOTAL_BATCHES}\n')


    # job script
    job_script_content = config['job_script']
    with open("job.sh", "w", encoding="utf-8") as f:
        f.write(job_script_content)
    os.chmod("job.sh", 0o755)

    # KPT
    kpt_content = config['kpt']
    with open("KPT", "w", encoding="utf-8") as f:
        f.write(kpt_content)

    input_content = config['input']

    num_digits = len(str(TOTAL_TASKS))
    for batch in range(1, TOTAL_BATCHES + 1):
        print("-" * 30)
        print(f"Current batch: {batch}/{TOTAL_BATCHES}")
        print("-" * 30)
        
        for task in range(1, TASKS_PER_BATCH + 1):
            index = task - 1 + (batch - 1) * TASKS_PER_BATCH
            
            if index >= TOTAL_TASKS:
                break
            

            new_index = f"{index+1:0{num_digits+1}d}"

            dir_name = f"job_{new_index}"
            os.makedirs(dir_name, exist_ok=True)
            os.chdir(dir_name)


            # STRU
            perturbed_stru = f"../{all_files[index]}"
            merge_stru_by_index(perturbed_stru, '../B.stru', indices1, indices2, 
                                    pp, basis, 'STRU')

            # INPUT
            with open("INPUT", "w") as f:
                f.write(input_content)


            # submit job
            try:
                result = subprocess.run(["sbatch", "../job.sh"], 
                                      check=True, capture_output=True, text=True)
                print(f"submit job: {new_index}")
            except subprocess.CalledProcessError as e:
                print(f"Failed - {e}")
            
            os.chdir("..")
        
        # sleep if not the last batch
        if batch < TOTAL_BATCHES:
            print(f"wait {WAIT_TIME} seconds...")
            time.sleep(WAIT_TIME)
        
        print("-" * 30)
    
    print("All job submitted!")
    print("-" * 30)

def postprocess():
    
    # loading settings
    config = load_config()

    selected_indices = config['selected_indices']
    natom = len(selected_indices)

    all_files = sorted(glob.glob("STRU-*"))
    TOTAL_TASKS = len(all_files)
    num_digits = len(str(TOTAL_TASKS))

    for task in range(1, TOTAL_TASKS + 1):

        new_index = f"{task:0{num_digits+1}d}"
        dir_name = f"job_{new_index}"
        os.chdir(dir_name)

        out_folders = glob.glob("OUT.*")
        out_folder = out_folders[0]
        os.chdir(out_folder)

        atom_forces = parse_forces_from_file(selected_indices)

        # phonon.log
        log = f"""TOTAL ATOM NUMBER = {natom}

 #TOTAL-FORCE (eV/Angstrom)#
 -------------------------------------------------------------------------
     Atoms              Force_x              Force_y              Force_z 
 -------------------------------------------------------------------------\n"""

        with open("phonon.log", "w") as f:
            f.write(log)
            f.write(str(atom_forces))
            f.write("\n -------------------------------------------------------------------------")


        os.chdir("../..")

    result = subprocess.run("phonopy -f job_*/OUT*/phonon.log",
                           shell=True, check=True, capture_output=True, text=True)

    
    mesh_content = config['mesh.conf']
    with open("mesh.conf", "w") as f:
                f.write(mesh_content)

    result = subprocess.run("phonopy -t mesh.conf",
                           shell=True, check=True, capture_output=True, text=True)




def main():
    if len(sys.argv) != 2:
        print("usage: python3 selective_dynamics.py --submit")
        print("usage: python3 selective_dynamics.py --post")
        sys.exit(1)

    if sys.argv[1] == "--submit":
        submit_jobs()
    elif sys.argv[1] == "--post":
        postprocess()
    else:
        print(f"No such parameter: {sys.argv[1]}")

if __name__ == "__main__":
    main()