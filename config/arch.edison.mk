# arch.mk for BerkeleyGW codes
#
# Suitable for Edison at NERSC.
# This is the recommended arch.mk for Edsion.
#
# Supports MPI+OpenMP parallelism and HDF5.
# FFTs are provided by FFTW, and linear algebra by libSci.
#
# To compile BerkeleyGW, you'll need to:
# module swap PrgEnv-intel PrgEnv-gnu && module load cray-hdf5-parallel && module load fftw 
#
# JRD
# 2013, NERSC
#
# 2016-01-20: All tests pass @ r6960 (FHJ)

# Precompiler options
#
COMPFLAG  = -DGNU
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG -DVERBOSE

FCPP    = /usr/bin/cpp -C
F90free = ftn -ffree-form -ffree-line-length-none -fopenmp -fno-second-underscore
LINK    = ftn -fopenmp
# FHJ: -funsafe-math-optimizations breaks Haydock and doesn't give any significant speedup
FOPTS   = -O3 -funroll-loops #-funsafe-math-optimizations -ffast-math
# FHJ: Useful flags for debugging:
#FOPTS   = -O0 -g -fbounds-check -fbacktrace -Wall -finit-real=nan -finit-integer=-100 -ffpe-trap=invalid -fopt-info-vec-all=vec
FNOOPTS = $(FOPTS)
MOD_OPT = -J 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = CC
C_COMP  = cc
C_LINK  = CC
C_OPTS  = -O3 -ffast-math
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
#
FFTWLIB      = -L$(FFTW_LIB) -lfftw3_omp -lfftw3 
FFTWINCLUDE  = $(FFTW_INC)
PERFORMANCE  = 

HDF5PATH     = $(HDF5_DIR)
HDF5LIB      = $(HDF5PATH)/lib/libhdf5hl_fortran.a \
               $(HDF5PATH)/lib/libhdf5_hl.a \
               $(HDF5PATH)/lib/libhdf5_fortran.a \
               $(HDF5PATH)/lib/libhdf5.a -lz -ldl
HDF5INCLUDE  = $(HDF5PATH)/include

LAPACKLIB = 

TESTSCRIPT = sbatch edison.scr
