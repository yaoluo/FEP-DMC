## Perturbo

Perturbo is an open-source software package for first-principles calculations of charge transport and ultrafast carrier dynamics in solid state materials, including metals, semiconductors, insulators, and 2D materials. In the current version, Perturbo mainly computes electron-phonon (e-ph) interactions and phonon limited transport properties in the framework of the Boltzmann transport equation (BTE). These include the carrier mobility, electrical conductivity, and Seebeck coefficient. Perturbo can also compute the ultrafast carrier dynamics (for now, with fixed phonon occupations) by explicitly evolving in time the electron BTE. 


Perturbo is written in Fortran with hybrid parallelization (MPI plus OpenMP). The main output format is HDF5, which is easily portable from one machine to another and is convenient for postprocessing using high-level languauges (e.g., Python and Julia).  Perturbo has a core software, called `perturbo.x`, for electron dynamics calculations and an interface software, called `qe2pert.x`, to read output files of Quantum Espresso (version 6.4.1) and Wannier90 (W90, version 3.0 or higher). The `qe2pert.x` interface software generates a HDF5 file, which is then read from the core `perturbo.x` software. In principle, any other third-party density functional theory (DFT) codes (e.g., VASP) can use Perturbo as long as the interface of the DFT codes can prepare a HDF5 output format for Perturbo to read.


## Features

Perturbo has the following stable features:

* Phonon-limited carrier mobility, electrical conductivity and Seebeck coefficient. 
* Phonon-limited carrier mean free path and relaxation times
* Imaginary part of e-ph self-energy and e-ph scattering rates
* e-ph matrix elements for nonpolar and polar materials, and their Wannier interpolation
* Interpolated electronic band structure and phonon dispersion
* Ultrafast carrier dynamics with fixed phonon occupation


## Installation

PERTURBO currently uses a small number of subroutines from the Quantum Espresso (QE) packages. Therefore, it needs to be compiled on top of QE. We assume that the users have already compiled QE successfully.

Clone from GitHub (or extract .tar.gz) into the QE directory, which creates a _"perturbo"_ directory containing the source files, utitlities, and user's manual. 

Perturbo uses the config file _make.inc_ of QE for most of the compiler options. 
The config file _make.sys_ inside the directory _"perturbo"_ specifies additional options required by Perturbo. 

Modify _make.sys_ to make it suitable for your system, such as the OpenMP options and path to the HDF5 library (not needed if HDF5 library is already specified in _make.inc_ of QE).

Once the file _make.sys_ has been modified, you are ready to compile PERTURBO.

```bash
$ make
```

After the compiling, a directory called _"bin"_ is generated, which contains two executables, `perturbo.x` and `qe2pert.x`.



Full documentation and tutorial examples are available in the user's manual and at <https://perturbo-code.github.io>.



## References

For more details on the code structure of Perturbo, we refer the users to the manuscript accompying the source code: 

- Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Ivan Maliyov, Xiao Tong, Marco Bernardi, <i>"PERTURBO: A software package for ab initio electron-phonon interactions, charge transport and ultrafast dynamics"</i>. Preprint: <a href="https://arxiv.org/abs/2002.02045" target="_blank"> arXiv: 2002.02045</a> (2020).

**When using results from PERTURBO in your publications, please cite the PERTURBO paper given above and acknowledge the use of the code.**
