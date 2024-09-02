#!/bin/bash 

Htable=../../pert-Htable

work=$PWD
echo $work


T=50
ik=0
nsv=20

for imc in 1 
do
for nk in  20 40 60    
do

cd $work

mkdir $nk
cd $nk

ln -sf ../../lif-sp3_epwan.h5 

echo "1 T 
${T}.00 11.5 1.0E+18" > temper.in
ln -sf $epw 
echo "&perturbo
 prefix='lif-sp3'
 calc_mode = 'diagmc-EZ'

 DMC_Method = 0
 band_min = 1
 band_max = 1
 zeroTMC = .true.
 read_H = .true.
 nsvd = $nsv 
 nk_svd = ${nk}
 tauMin = 0.d0 
 alpha_frohlich = 1.0
 svd_dir = '${Htable}/gkq-${nk}'
 
 ftemper = 'temper.in'
 phfreq_cutoff = 1.0 !mev
 print_sign = .false. 

/">pert.in


echo "Nmcmc(1e4)      Px  Py   Pz      mu(eV)  se_check   nq_se_checkmaxOrder 
1000            $ik  $ik   0     0.0          0    50000 500
types of update
7
prob for each update, for green-function, there are 7 updates  
0.1 0.1 0.1 0.1 0.10 0.10 0.1"          >diagMC.in

export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread

nprocs=1
mpirun -np 1 perturbo.x -npools 1 -i pert.in 

done 
done 



