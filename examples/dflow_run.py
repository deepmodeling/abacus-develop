#!/usr/bin/env python3

import sys,copy,re
import argparse
import json
from typing import List
from dflow import (
    Workflow,
    Step,
    argo_range,
    SlurmRemoteExecutor,
    upload_artifact,
    download_artifact,
    InputArtifact,
    OutputArtifact,
    ShellOPTemplate
)
from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices
)
import time
import subprocess, os, shutil, glob
from pathlib import Path
from typing import List
from dflow.plugins.lebesgue import LebesgueContext
from dflow import config, s3_config
config["host"] = "http://xxx.xxx"
config["k8s_api_server"] = "https://xxx.xxx"
s3_config["endpoint"] = "xxx.xxx"

class AbacusExample(OP):
    """
    class for AbacusExample
    """
    def __init__(self, infomode=1):
        self.infomode = infomode

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "env_cmd": str,
            "exepath": str,
            "nproc": int,
            "threads": int,
            "input": Artifact(Path),
            "tests": Artifact(Path),
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "output": Artifact(Path)
        })

    @OP.exec_sign_check
    def execute(self, op_in: OPIO) -> OPIO:
        cwd = os.getcwd()
        os.chdir(op_in["input"])
        cmd= op_in["env_cmd"] + " && " \
             "echo 'ABACUS_PATH=" + op_in["exepath"] + "' >> ../SETENV && " + \
             "echo 'ABACUS_NPROCS=%d"%op_in["nproc"] + "' >> ../SETENV && " + \
             "echo 'ABACUS_THREADS=%d"%op_in["threads"] + "' >> ../SETENV && " + \
             "bash runall.sh" 
        subprocess.call(cmd, shell=True)
        os.chdir(cwd)
        op_out = OPIO({
            "output": op_in["input"]
        })
        return op_out


### find all directories where runall.sh has been prepared
def FindRunDirs(path):
    os.chdir(path)
    _dir = glob.glob('*')
    _dir.sort()
    find_dir = []
    for idir in range(len(_dir)):
        if os.path.isdir(_dir[idir]):
            find_dir.append(_dir[idir])
    run_dir = []
    for idir in range(len(find_dir)):
        os.chdir(find_dir[idir])
        _files = glob.glob('*')
        _files.sort()
        if 'runall.sh' in _files:
            run_dir.append(find_dir[idir])
        os.chdir('../')
    return run_dir

### find all pp and orb files needed
def FindPPORB(path,pp_orb):
    os.chdir(path)
    _stru = glob.glob('*/STRU')
    _stru.sort()
    #print(_stru)
    pporb_list = []
    for istru in range(len(_stru)):
        with open(_stru[istru],'r') as fp:
            ilines = fp.read()
            for ipp in range(len(pp_orb)):
                if pp_orb[ipp] in ilines:
                    pporb_list.append(pp_orb[ipp])
    os.chdir('../')
    if len(pporb_list)==0:
        raise RuntimeError("Can't find pp_orb file in ../tests/PP_ORB")
    else:
        #print(pporb_list)
        pporb_list_tmp = list(set(pporb_list))
        return pporb_list_tmp

# funcs related to post process
def FindDirs(path):
    os.chdir(path)
    _dir = glob.glob('*')
    _dir.sort()
    find_dir = []
    for idir in range(len(_dir)):
        if os.path.isdir(_dir[idir]):
            find_dir.append(_dir[idir])
    os.chdir('../')
    return find_dir

def CheckRunningLogs(path):
    os.chdir(path)
    _files = glob.glob('running*')
    file_status = {}
    for ifile in range(len(_files)):
        abspath = os.path.join(os.path.abspath('./'),_files[ifile])
        with open(_files[ifile],'r') as ff:
            final_line = ff.read().split('\n')[-2]
            if (re.search('Total  Time',final_line)):
                file_status[abspath] = "Ok\t"
            else:
                file_status[abspath] = "Failed\t"
    os.chdir('../')
    return file_status

def CheckJobStatus(path):
    sys_list = FindDirs(path)
    os.chdir(path)
    for isys in range(len(sys_list)):
        #print(sys_list[isys])
        out_dir = FindDirs(sys_list[isys])
        if (len(out_dir)==0):
            print("Failed\t No OUT.ABACUS directory in",path+"/"+sys_list[isys])
        else:
            os.chdir(sys_list[isys])
            files_status = CheckRunningLogs(out_dir[0])
            for key in files_status:
                print(files_status[key],key)
            os.chdir('../')
    os.chdir('../')


def main(run_params,run_dir):
    jcwd = os.getcwd()
    pp_dir = os.path.join(jcwd,'../tests/PP_ORB/')
    os.chdir(pp_dir)
    pp_orb = glob.glob('*')
    os.chdir(jcwd)
    abacus = PythonOPTemplate(AbacusExample,image=run_params["LBG_IMAGE"] ,command=['python3'])
    job_list = []
    os.makedirs('PP_ORB', exist_ok = True)
    jobs_dict={}
    for idir in range(len(run_dir)):
        #print("")
        #print(run_dir[idir])
        pporb_list = FindPPORB(run_dir[idir],pp_orb)
        #print(pporb_list)
        pporb_files = []
        for ipp in range(len(pporb_list)):
            shutil.copy2(os.path.join(pp_dir,pporb_list[ipp]),os.path.join(jcwd,'PP_ORB'))
            pporb_files.append(os.path.join(jcwd,'PP_ORB',pporb_list[ipp]))
        jpath = os.path.join(jcwd,run_dir[idir])
        jobs = [jpath]
        step_name = "ABACUS-EXAMPLE"+str(len(job_list))
        jobs_dict[step_name] = run_dir[idir]
        abacus_example = Step(name=step_name,
                template=abacus,
                parameters={
                    "exepath":run_params["EXE_PATH"],
                    "nproc"  :run_params["NPROCS"],
                    "threads":run_params["NTHREADS"],
                    "env_cmd":run_params["ENV_CMD"]},
                artifacts={
                    "input":upload_artifact(jobs),
                    "tests":upload_artifact(pporb_files)})
        job_list.append(abacus_example)
    
    wf = Workflow(name="abacus-functions", context=lebesgue_context, host="http://xxx.xxx")
    wf.add(job_list)
    wf.submit()

    f1 = open('check.log','w')
    sys.stdout = f1
    print("Jobs submmited!")
    print("Starts waiting...")
    for key in jobs_dict:
        print(key, jobs_dict[key])
    sys.stdout.flush()

    step_status = [0]*len(job_list)
    while wf.query_status() in ["Pending","Running"]:
        for ii in range(len(job_list)):
            step_name = "ABACUS-EXAMPLE"+str(ii)
            if len(wf.query_step(name=step_name))==0:
                continue
            else:
                step = wf.query_step(name=step_name)[0]
                if step.phase == 'Succeeded' and step_status[ii] == 0:
                    download_artifact(step.outputs.artifacts["output"])
                    step_status[ii] = 1
                    print(step_name, jobs_dict[step_name],' Finished!')
                    CheckJobStatus(os.path.join(jcwd,jobs_dict[step_name]))
                    sys.stdout.flush()
        time.sleep(4)
    assert(wf.query_status() == 'Succeeded')
    for ii in range(len(job_list)):
        step_name = "ABACUS-EXAMPLE"+str(ii)
        step = wf.query_step(name=step_name)[0]
        if step.phase == 'Succeeded' and step_status[ii] == 0:
            download_artifact(step.outputs.artifacts["output"])
            step_status[ii] = 1
            print(step_name, jobs_dict[step_name],' Finished!')
            CheckJobStatus(os.path.join(jcwd,jobs_dict[step_name]))
            sys.stdout.flush()
    print("All done!!")
    f1.close()
    shutil.rmtree('PP_ORB')

def RandomDisturbParser():
    parser = argparse.ArgumentParser(
        description="Script to automatically run ABACUS examples with dflow")
    parser.add_argument('-run', '--run', type=int,
            default=0, help='run dflow: 1 run, default 0')
    parser.add_argument('-post', '--post', type=int,
            default=0, help='checkout job status: 1 check, default 0')
    return parser.parse_args()


if __name__ == "__main__":
    args = RandomDisturbParser()
    run = args.run
    post = args.post
    lebesgue_context = LebesgueContext(
        username="xxx@xxx.xxx",
        password="xxxxxx",
        executor="lebesgue_v2",
        extra='{"scass_type":"c16_m32_cpu","program_id":xxx}',
        app_name='Default',
        org_id='123',
        user_id='456',
        tag='')
    run_params={}
    run_params["LBG_IMAGE"] = "ABACUS"
    run_params["EXE_PATH"]  = 'abacus'
    run_params["NPROCS"]    = 8
    run_params["NTHREADS"]  = 1
    run_params["ENV_CMD"]   = "ulimit -s unlimited && " + \
                              ". /opt/intel/oneapi/setvars.sh"
    cwd = os.getcwd()
    ### find all directories where runall.sh has been prepared
    run_dir = FindRunDirs(cwd)
    ### or add examples directory manually
    #run_dir = ["electrostatic_potential","scf","wfc"]
    if run == 1:
        main(run_params,run_dir)
    if post == 1:
        for idir in range(len(run_dir)):
            CheckJobStatus(run_dir[idir])
