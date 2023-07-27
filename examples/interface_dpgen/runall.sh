# install latest dpgen and dpdata
latest_dpgen=1
if [ $latest_dpgen == 1 ];then
	pip install -U dpdispatcher bohrium-sdk setuptools_scm oss2
	pip install git+https://gitee.com/deepmodeling/dpdata.git@devel -i https://pypi.tuna.tsinghua.edu.cn/simple
	pip install git+https://gitee.com/deepmodeling/dpgen.git@devel -i https://pypi.tuna.tsinghua.edu.cn/simple
fi
dpgen -h

for ifolder in "init_and_run" \
	"autotest"
do
	cd $ifolder
	
	# modify machine.json by Environment Variables
	python ../modify_machine.py

	start_time=$(date +%s)
	bash run.sh
	echo $? > run.result
	end_time=$(date +%s)
	elapsed_time=$((end_time - start_time))
	echo ${elapsed_time} > run.time

	cd ..
done

