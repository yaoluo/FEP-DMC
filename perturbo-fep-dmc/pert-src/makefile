 
include ../make.sys


PERTMOD = \
MyLAPACK.f90 \
HDF5_utils.f90 \
wigner_seitz_cell.f90 \
polar_correction.f90 \
electron_wannier.f90 \
force_constant.f90 \
elph_matrix_wannier.f90 \
epwan_hdf5_io.f90 \
pert_const.f90 \
qe_mpi_mod.f90 \
pert_utils.f90 \
pert_param.f90  \
pert_data.f90   \
pert_output.f90 \
band_structure.f90 \
phonon_dispersion.f90 \
elphon_coupling_matrix.f90 \
boltz_utils.f90 \
vector_list.f90 \
boltz_grid.f90  \
imaginary_selfenergy.f90 \
boltz_trans_output.f90 \
boltz_trans_mod.f90 \
boltz_scatter.f90 \
boltz_scatter_integral.f90 \
boltz_grid_neighbors.f90 \
boltz_external_field.f90 \
boltz_dynamics_solver.f90 \
boltz_dynamics_mod.f90 \
selfenergy.f90 \
cumulant_expansion.f90 \
cumulant_utility.f90  \
eph_svd.f90 \
random_tool.f90 \
diagMC_pre.f90 \
diagMC.f90 \
diagMC_debug.f90 \
diagMC_JJ_updates.f90 \
diagMC_updates.f90 

PERTSRC = \
environment_setup.f90 \
weight_dos.f90 \
load_data_pert.f90 \
calc_bands.f90 \
calc_ephmat.f90 \
calc_electron_mfp.f90 \
real_selfenergy.f90 \
electron_imsigma.f90 \
electron_resigma.f90 \
boltz_setup.f90 \
transport.f90 \
trans_postproc.f90 \
carrier_dynamics_run.f90 \
carrier_dynamics_postproc.f90 \
calc_selfenergy.f90 \
process_spectral.f90 \
calc_spectral_cumulant.f90 \
calc_trans_cumulant.f90 \
calc_trans_cumulant_opcond.f90 \
diagMC_gt.f90 \
diagMC_JJ.f90 \
tabulate_svd.f90 \
perturbo.f90

OBJ = $(PERTMOD:.f90=.o) $(PERTSRC:.f90=.o)

Include_FFT=-I/usr/include

EXE = perturbo.x

debug=-fcheck=all -fbacktrace -g

default: $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(PHMODS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS) -lstdc++
se:
	$(LD) $(debug)  $(MODFLAG) -c selfenergy.f90

diagmc:
	$(LD) $(F90FLAGS)  $(MODFLAG) -c diagMC.f90 diagMC_updates.f90 diagMC_JJ.f90 

eph_svd:
	$(LD) $(debug)  $(MODFLAG) -c eph_svd.f90

%.o: %.f90 
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<
clean:
	rm -f *.o *~ *.d *.mod *.i

