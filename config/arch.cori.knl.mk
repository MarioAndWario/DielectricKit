# arch.mk for BerkeleyGW codes using Intel compiler.
# Features: MPI, OpenMP, LibSci, HDF5
#
# Suitable for Cori Phase 1 at NERSC
#
# *** This is the recommended arch.mk for Cori Phase 1 ***
#
# Do:
# module swap intel intel/15.0.1.133
# module load fftw cray-hdf5-parallel
#
# Felipe H. da Jornada (2015)
#
# Status:
# Passes most of the testsuite, except for the Haydock code, which has some
# instabilities issues in its algorithm. One of the graphene tests also fails
# because it finds a different coarse-grid k-point, but the answer doesn't
# change much. This is likely just an ambiguity in the algorithm.

# Precompiler options
#
COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG -DVERBOSE

#FCPP    = fpp -free -C # -P
FCPP    = /usr/bin/cpp -C -nostdinc
F90free = ftn -xMIC-AVX512 -free -qopenmp
LINK    = ftn -xMIC-AVX512 -qopenmp
#FOPTS   =  -O3 -xHost -ip -no-ipo -no-prec-div -fp-model fast=2 -align array64byte
#FNOOPTS =  -O2 -xHost -no-ip -no-ipo -no-prec-div -fp-model fast=2 -align array64byte
#FOPTS =  -g -O0 -check all -Warn all -traceback
#FOPTS =  -g -O0 -check none -Warn all -traceback
FOPTS   = -fast -no-ip -no-ipo -align array64byte
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = CC -xMIC-AVX512 -qopenmp
C_COMP  = cc -xMIC-AVX512 -qopenmp
C_LINK  = CC -xMIC-AVX512 -qopenmp
C_OPTS  = -fast -no-ip -no-ipo -align
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

FFTWPATH     =
FFTWLIB      = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a \
               $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl

LAPACKLIB = $(FFTWLIB)

# Math Libraries
#
# FHJ: LibSci already provides FFT routines. But are they threaded?!
#FFTWLIB      = -L$(FFTW_DIR)/lib -lfftw3_omp -lfftw3

FFTWINCLUDE  = $(MKLROOT)/include/fftw/ #$(FFTW_INC)
PERFORMANCE  = 

HDF5_LDIR    =  $(HDF5_DIR)/lib
# We use LibSci instead of intel MKL (uncomment -mkl )
HDF5LIB      =  $(HDF5_LDIR)/libhdf5hl_fortran.a \
                $(HDF5_LDIR)/libhdf5_hl.a \
                $(HDF5_LDIR)/libhdf5_fortran.a \
                $(HDF5_LDIR)/libhdf5.a -lz -ldl # -mkl
HDF5INCLUDE  = $(HDF5_DIR)/include

TESTSCRIPT = sbatch cori2.scr
