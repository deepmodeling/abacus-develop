This is a test of dpgen+abacus.

The usage of dpgen can be found in
https://docs.deepmodeling.com/projects/dpgen/en/latest/overview/overview.html

Before run the test, you need install dpgen and dpdata:
`pip install dpgen`
`pip install dpdata`

You need to modify the machine.json file in each sub-path, and a detail of machine.json 
can be found in https://docs.deepmodeling.com/projects/dpgen/en/latest/init/init-bulk-mdata.html.
Here, the machine.json file is a format based on bohrium platform.

There has two ways to run this test:
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
		`cd ..`

