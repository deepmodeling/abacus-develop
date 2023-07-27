This is a test of dpgen+abacus.

This test will use bohrium machine to run ABACUS, you need to have bohrium
account. The ABACUS bohrium image can be found in bohrium website, or you can 
use "registry.dp.tech/deepmodeling/abacus-intel:latest", which is the latest
ABACUS compiled by intel.

Before run the tests, you need to install dpgen and dpdata:
`pip install dpgen`
`pip install dpdata`


There are two ways to run this test:
1. Autotest by runall.sh:
	a. Set enviroment variables:
		`export BOHRIUM_USERNAME=<replace by your bohrium username>`
		`export BOHRIUM_PASSWORD=<replace by your bohrium password>`
		`export BOHRIUM_PROJECT_ID=<replace by your bohrium project id>`
		`export ABACUS_IMAGE=<replace by ABACUS bohrium image>`	
	b. Run runall.sh:
		`bash runall.sh`
	c. After finishing the running, there are two files will be produced in sub-path: 
		- run.time
		- run.result # if the run is finished normally, the run.result should be 0.

2. Manually run the tests:
	a. Modify "machine.json" file in each test folder, and fill in the values 
           of "email", "password", "program_id", and "image_name" in `fp` section.
	b. Do the `init_bulk` and `run` tests by `bash run.sh` or run before commands:
		`cd init_and_run`
 		`dpgen init_bulk init.json machine.json`
  		`dpgen run run_param.json machine.json`
		`cd ..`
	c. Do the `autotest` test by `bash run.sh` or run below commands:
		`cd autotest`
		`dpgen autotest make relax.json`
 		`dpgen autotest run relax.json machine.json`
		`dpgen autotest post relax.json`
		`dpgen autotest make property.json`
		`dpgen autotest run property.json machine.json`
		`dpgen autotest post property.json`

