## First-principles DMC
 
Our First-prinicples DMC(FEP-DMC) is a module of the open-source package of Perturbo. For more information of Perturbo, please visit our github website https://perturbo-code.github.io/ . 

## System Requirements
### Hardware requirement 
The minimal RAM is 64GB to run the full example. 
### OS reuqirement 
This package is supported for Linux only. It has been tested on the following systems. 

Linux: Ubuntu 20.04

### External package dependencies 
PERTURBO currently uses a small number of subroutines from the Quantum Espresso (QE) packages. Therefore, it needs to be compiled on top of QE. We assume that the users have already compiled QE successfully.

FEP-DMC is built on top of QE version 6.5. See the downloading link [https://gitlab.com/QEF/q-e/-/tags/qe-6.5](https://gitlab.com/QEF/q-e/-/tags/qe-6.5). Other versions of QE are not supported currently. Please note that QE should be compiled with the Intel Fortran compiler and HDF5 (recommended version >= 1.8.5). HDF5 should be compiled with the same Intel Fortran compiler. HDF5 should be configured with `--enable-fortran --enable-fortran2003`. The website [https://perturbo-code.github.io/MB-group-resource-sharing/docs/compile-codes/](https://perturbo-code.github.io/MB-group-resource-sharing/docs/compile-codes/) is useful for installation of HDF5 and QE on your own HPC or Ubuntu, although it is for our own HPC.

In addition, the installation of QE also requires the Intel Fortran compiler, `BLAS`, and `LAPACK`. We recommend `Intel oneAPI` [https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html), which includes `BLAS` and `LAPACK` via the MKL in oneAPI.
 

## Installation Guide 
On my ubuntu workstation, this whole installation, including compiling QE and perturbo, takes around 15 mins. 
Assume that HDF5 and Intel oneAPI have been installed and activated.
Go to the directory of QE and copy the `perturbo-fep-dmc` there.
```bash 
$ cd QE_direcroty/
$ cp -rf /some_directory/perturbo-fep-dmc  .
$ ls 
archive          Doc                    KS_Solvers  perturbo-fep-dmc  TDDFPT
atomic           environment_variables  LAXlib      PHonon            test-suite
clib             EPW                    License     PlotPhon          upftools
configure        FFTXlib                logo.jpg    PP                UtilXlib
CONTRIBUTING.md  GUI                    LR_Modules  pseudo            XSpectra
COUPLE           GWW                    Makefile    PW
CPV              HP                     make.inc    PWCOND
dev-tools        include                Modules     QHA
dft-d3           install                NEB         README.md
$

```


In an Ubuntu 20.04 OS system, I compile QE with the following commands:
```bash
# my HDF is installed in /home/yaoluo/software/hdf5-1.8.18-ifort
$ ./configure MPIF90=mpiifort --with-hdf5=/home/yaoluo/software/hdf5-1.8.18-ifort
$ make pw ph pp 
```

Perturbo uses the config file make.inc of QE for most of the compiler options.
In the perturbo-fep-dmc/ directory, we provide you with a make.sys for setting up the compilation.

```bash
$ cd perturbo-fep-dmc

# modify the file make.sys
>> vim make.sys
------------------------------------------------------------------------------------------
......
# modify the paths for HDF5 library
# my HDF is installed in /home/yaoluo/software/hdf5-1.8.18-ifort
......
IFLAGS += -I/home/yaoluo/software/hdf5-1.8.18-ifort/include
HDF5_LIB = -L/home/yaoluo/software/hdf5-1.8.18-ifort/lib -lhdf5 -lhdf5_fortran
......
------------------------------------------------------------------------------------------
```

and make the necessary changes to it. You will then be ready to compile PERTURBO:

```bash
$ make perturbo 

#add perturbo bin to the path for the example
$ export PATH="<QE-directory>/perturbo-fep-dmc/bin:$PATH"
```

After the compiling, a directory called _"bin"_ is generated, which contains one executable, `perturbo.x`

## example 
This whole process takes around 20 mins on my workstation.  
We are going to calculate the formation energy of electron polarons in LiF with compressed electron-phonon interactions.
This example is dedicated to reproducing some points of Fig.2(b) in the main text.
Before calculation, one should download the "lif-sp3_epwan.h5" file via this link " https://www.dropbox.com/scl/fi/40djsz04bkrdyrr8pmzf2/DMC-Demo.zip?rlkey=kyil1w6j7f5tt5ii84ii5oi4o&dl=0 ". 
This file contains the electron-phonon interactions. 
There are three steps.
### 1. Compress e-ph 
```bash 
$ cd example/pert-svd
$ ln -sf ../lif-sp3_epwan.h5
$ ./run.sh 
```
### 2. Construct Hamiltonian table on uniform grid
```bash 
$ cd example/pert-Htable
$ ln -sf ../lif-sp3_epwan.h5
$ ./run-nk.sh 
```
### 3. run FEP-DMC
```bash 
$ cd example/E-nk

# caluclate polaron binding energy for different size of grid 
# 20^3, 40^3, 60^3, 80^3 
$ ./loop-nk.sh
$ ls 
20  40  60  loop-nk.sh
# 20 for 20^3 grid; 40 for 40^3 grid; 60 for 60^3 grid
$ cd 20 
$ ls 
chemical_pot.dat  orderStat.dat-1  pert.in           sign.dat-1
diagMC.in         parameter.dat    PolaronLat.dat-1  temper.in
lif-sp3_epwan.h5  parameter.dat-1  PolaronWF.dat-1   Znph.dat-1
# parameter.dat-1 contains polaron ground state energy 
# expected output for 20^3 grid 
$ cat parameter.1
    #mu = 
    0.00000E+00
 #alpha = 
    0.10000E+01
    #Ep = 
    0.93549E+01
 #Etrue = 
    0.9107056695E+01
 #<g/Z> = 
    0.9400030978E+00
  #<Z0> = 
    0.2238440137E+00
# Here, Ep = the band electron energy, Etrue = polaron ground state energy; 
# The polaron binding energy equals `Etrue - Ep = -0.25 eV` for 20^3 grid. 
# Znph.dat-1 contains the n phonon amplitude, the last two lines are polaron ground state energy 
# PolaronWF.dat-1 contains the reduced density matrix in electron basis 
# 
```


## References

For more details on the code structure of Perturbo, we refer the users to the manuscript accompying the source code: 

- Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Ivan Maliyov, Xiao Tong, Marco Bernardi, <i>"PERTURBO: A software package for ab initio electron-phonon interactions, charge transport and ultrafast dynamics"</i>. Preprint: <a href="https://arxiv.org/abs/2002.02045" target="_blank"> arXiv: 2002.02045</a> (2020).

**When using results from PERTURBO in your publications, please cite the PERTURBO paper given above and acknowledge the use of the code.**
