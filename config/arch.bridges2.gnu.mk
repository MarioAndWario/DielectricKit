# arch.mk for BerkeleyGW codes
#Currently Loaded Modules:
#  1) allocations/1.0   3) openmpi/4.0.5-gcc10.2.0               5) mkl/2020.4.304   7) anaconda3/2020.11
#  2) gcc/10.2.0        4) phdf5/1.10.7-openmpi4.0.5-gcc10.2.0   6) fftw/3.3.8
# suitable for Bridges-2@PSC

COMPFLAG  = -DGNU
PARAFLAG  = -DMPI #-DOMP
MATHFLAG  = -DUSESCALAPACK -DHDF5 -DUSEFFTW3 -DUNPACKED -DUSEELPA
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

VEC_FLAGS = -ffast-math -march=znver1 -mtune=znver1 -mfma -mavx2 -m3dnow -fomit-frame-pointer
FOPT_FLAGS = $(VEC_FLAGS) -funroll-loops  -funsafe-loop-optimizations -funsafe-math-optimizations -ftree-vect-loop-version -ftree-vectorize  -ffree-form -ffree-line-length-none -ffpe-summary=none
COPT_FLAGS = $(VEC_FLAGS) -funroll-loops  -funsafe-loop-optimizations -funsafe-math-optimizations -ftree-vect-loop-version -ftree-vectorize

MPI_DIR = $(MPI_GNU_DIR)
FCPP    = cpp -C -nostdinc
F90free = $(MPI_DIR)/mpif90
LINK    = $(MPI_DIR)/mpif90
FOPTS   = -O3 -fallow-argument-mismatch $(FOPT_FLAGS)

#FOPTS   = -O0 -fallow-argument-mismatch -Wall  -fcheck=all  -pedantic  -fbacktrace
FNOOPTS = $(FOPTS)
MOD_OPT = -J
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = $(MPI_DIR)/mpiCC
C_COMP  = $(MPI_DIR)/mpicc
C_LINK  = $(MPI_DIR)/mpiCC
C_OPTS  = -O2 $(COPT_FLAGS)
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
FFTWPATH     = $(FFTW3_DIR)
FFTWLIB = -L$(FFTWPATH) -lfftw3 -lfftw3_mpi
FFTWINCLUDE  = $(FFTWPATH)/include

LAPACKLIB    = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 -lpthread -lm -ldl

SCALAPACKLIB = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64

HDF5LIB      = -L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
HDF5INCLUDE  = $(HDF5_DIR)/include

ELPAINCLUDE = $(ELPA_DIR)/include/elpa-2017.11.001/modules -I$(ELPA_DIR)/include/elpa-2017.11.001/elpa
ELPALIB = $(ELPA_DIR)/lib/libelpa.a

SLATECLIB = $(SLATEC_DIR)/libslatec.a

TESTSCRIPT = sbatch bridges2.scr
