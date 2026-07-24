sed -i "s/^lambda.*/lambda = $1/g" input
mpirun -np 16 pw.x -i input | tee out.log && rm -rf pwscf.save