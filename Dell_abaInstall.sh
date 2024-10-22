rm -r build

cmake -B build -DCMAKE_INSTALL_PREFIX=/home/jghan/software/abacus/ -DCMAKE_CXX_COMPILER=mpiicpc -DELPA_DIR=/home/jghan/resource/elpa-2021.11/ -DCMAKE_INSTALL_PREFIX=/home/jghan/software/abacus-develop/ -DLibxc_DIR=/home/jghan/resource/libxc-6.2.2/ -DLIBRI_DIR=/home/jghan/resource/LibRI -DLIBCOMM_DIR=/home/jghan/resource/LibComm -DCEREAL_INCLUDE_DIR=/home/jghan/resource/cereal/cereal-1.3.2/include

cmake --build build -j 12 2>job.err

cmake --install build
