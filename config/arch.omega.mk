# arch.mk for BerkeleyGW codes
#Currently Loaded Modulefiles:
#  1) emacs/25.1            3) python/3.6            5) hdf5/1.8.18-intel-p   7) fftw/3.3.6-intel
#  2) intel/2016.4.072      4) openmpi/2.0.2-intel   6) mkl/2016.4.072
# suitable for Omega at LBNL with HDF5

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI #-DOMP
MATHFLAG  = -DUSESCALAPACK -DHDF5 -DUSEFFTW3 -DUNPACKED -DUSEELPA
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG

VEC_FLAGS = -xCORE-AVX512

FCPP    = cpp -C -nostdinc
F90free = mpif90 -free
LINK    = mpif90 -static-intel -static-libgcc -static-libstdc++  -Wl,--as-needed -liomp5 -Wl,--no-as-needed
FOPTS   = -O3 $(VEC_FLAGS) -no-ipo -ip -align array64byte -threads -heap-arrays 4096 -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -assume byterecl -assume nobuffered_io -assume nobuffered_stdout -no-prec-div

#FOPTS   = -g -O0 -check all -Warn all -traceback
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiCC
C_COMP  = mpicc
C_LINK  = mpiCC -Wl,--as-needed -liomp5 -Wl,--no-as-needed
C_OPTS  = -O3 $(VEC_FLAGS) -no-ipo -ip -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -fno-alias -ansi-alias
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
MKLPATH      = $(MKL_DIR)
FFTWPATH     = $(FFTW3_DIR)
FFTWLIB = -L$(FFTWPATH) -lfftw3 -lfftw3_mpi
FFTWINCLUDE  = $(FFTWPATH)/include

LAPACKLIB    = -L$(MKLPATH)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64  -lpthread -lm -ldl

SCALAPACKLIB = -L$(MKLPATH)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64

HDF5LIB      = -L$(HDF5_DIR)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
HDF5INCLUDE  = $(HDF5_DIR)/include

ELPAINCLUDE = $(ELPA_DIR)/include/elpa-2017.11.001/modules -I$(ELPA_DIR)/include/elpa-2017.11.001/elpa
ELPALIB = $(ELPA_DIR)/lib/libelpa.a

SLATECLIB = $(SLATEC_DIR)/lib/libslatec.a

TESTSCRIPT = sbatch hbar.scr
