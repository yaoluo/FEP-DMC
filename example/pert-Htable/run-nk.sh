#!/bin/bash 



work=$PWD
echo $work

nsv=20
for nk in  20 40 60 
do

echo $nk
mkdir gkq-$nk
cd gkq-$nk

ln -sf ../../pert-svd/gkq/eph_S.dat
ln -sf ../../pert-svd/gkq/eph_U.dat
ln -sf ../../pert-svd/gkq/eph_V.dat
ln -sf ../../pert-svd/gkq/g_Rp.dat
cd ..

echo "&perturbo 
 prefix='lif-sp3'
 calc_mode = 'tabulate-H'
 nk_svd    = $nk

 !!! svd setting !!! 
 read_svd = .true. 
 apply_svd = .true.
 svd_dir = './gkq-${nk}'
 nsvd = $nsv
/">pert.in

export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread

nprocs=1
mpirun -np ${nprocs} perturbo.x -npools ${nprocs} -i pert.in 

python3 alias.py gkq-${nk}

done 

