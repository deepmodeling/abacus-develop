#!/usr/bin/env python3

import sys,copy,re
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
            "exepath": str,
            "nproc":int,
            "threads":int,
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
        cmd= "ulimit -s unlimited && " + \
             "echo 'ABACUS_PATH=" + op_in["exepath"] + "' >> ../SETENV && " + \
             "echo 'ABACUS_NPROCS=%d"%op_in["nproc"] + "' >> ../SETENV && " + \
             "echo 'ABACUS_THREADS=%d"%op_in["threads"] + "' >> ../SETENV && " + \
             "export OMPI_ALLOW_RUN_AS_ROOT=1 && " + \
             "export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 && " + \
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
                if (re.search(".+?"+pp_orb[ipp],ilines)):
                    pporb_list.append(pp_orb[ipp])
    os.chdir('../')
    if len(pporb_list)==0:
        raise RuntimeError("Can't find pp_orb file in ../tests/PP_ORB")
    else:
        #print(pporb_list)
        pporb_list_tmp = list(set(pporb_list))
        return pporb_list_tmp


def main(run_params):
    cwd = os.getcwd()
    pp_dir = os.path.join(cwd,'../tests/PP_ORB/')
    os.chdir(pp_dir)
    pp_orb = glob.glob('*')
    os.chdir(cwd)
    abacus = PythonOPTemplate(AbacusExample,image='ABACUS_GNU',command=['python3'])
    jcwd = os.getcwd()
    job_list = []
    run_dir = FindRunDirs(cwd)
    os.makedirs('PP_ORB', exist_ok = True)
    for idir in range(len(run_dir)):
        #print("")
        #print(run_dir[idir])
        pporb_list = FindPPORB(run_dir[idir],pp_orb)
        pporb_files = []
        for ipp in range(len(pporb_list)):
            shutil.copy2(os.path.join(pp_dir,pporb_list[ipp]),os.path.join(jcwd,'PP_ORB'))
            pporb_files.append(os.path.join(jcwd,'PP_ORB',pporb_list[ipp]))
        jpath = os.path.join(jcwd,run_dir[idir])
        jobs = [jpath]
        abacus_example = Step(name="ABACUS-EXAMPLE"+str(len(job_list)),
                parameters={"exepath":run_params[0],"nproc":run_params[1],"threads":run_params[2]},
                template=abacus,
                artifacts={"input": upload_artifact(jobs),"tests":upload_artifact(pporb_files)})
        job_list.append(abacus_example)
    
    wf = Workflow(name="abacus-test", context=lebesgue_context, host="http://xxx.xxx")
    wf.add(job_list)
    wf.submit()

    step_status = [0]*len(job_list)
    while wf.query_status() in ["Pending","Running"]:
        for ii in range(len(job_list)):
            if len(wf.query_step(name="ABACUS-EXAMPLE"+str(ii)))==0:
                continue
            else:
                step = wf.query_step(name="ABACUS-EXAMPLE"+str(ii))[0]
                if step.phase == 'Succeeded' and step_status[ii] == 0:
                    download_artifact(step.outputs.artifacts["output"])
                    step_status[ii] = 1
        time.sleep(4)
    assert(wf.query_status() == 'Succeeded')
    for ii in range(len(job_list)):
        step = wf.query_step(name="ABACUS-EXAMPLE"+str(ii))[0]
        if step.phase == 'Succeeded' and step_status[ii] == 0:
            download_artifact(step.outputs.artifacts["output"])
            step_status[ii] = 1
    shutil.rmtree('PP_ORB')

if __name__ == "__main__":
    lebesgue_context = LebesgueContext(
    username="xxx@xxx.xxx",
    password="xxxxxx",
    executor="lebesgue_v2",
    extra='{"scass_type":"c16_m32_cpu","program_id":xxx}',
    app_name='Default',
    org_id='123',
    user_id='456',
    tag='',
    )
    ABACUS_PATH='abacus'
    ABACUS_NPROCS=8
    ABACUS_THREADS=1
    run_params = [ABACUS_PATH,ABACUS_NPROCS,ABACUS_THREADS]
    main(run_params)
