rm -r build

cmake -B build -DCMAKE_INSTALL_PREFIX=/public1/home/t6s000394/jghan/software/abacus-develop/rdmft-abacus/ -DCMAKE_CXX_COMPILER=icpx -DMPI_CXX_COMPILER=mpiicpc -DELPA_DIR=/public1/home/t6s000394/jghan/software/elpa-2021.11/ -DLibxc_DIR=/public1/home/t6s000394/jghan/software/libxc-6.2.2/ -DLIBRI_DIR=/public1/home/t6s000394/jghan/software/LibRI -DLIBCOMM_DIR=/public1/home/t6s000394/jghan/software/LibComm -DCEREAL_INCLUDE_DIR=/public1/home/t6s000394/jghan/software/cereal/cereal-1.3.2/include

cmake --build build -j 24 2>job.err

cmake --install build
