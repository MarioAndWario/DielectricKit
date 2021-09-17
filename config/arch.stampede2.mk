# arch.mk for BerkeleyGW codes
# suitable for Stampede2 at TACC
#   1) git/2.24.1      3) cmake/3.16.1   5) TACC           7) libfabric/1.7.0   9) phdf5/1.8.16  11) python3/3.6.3
#   2) autotools/1.1   4) xalt/2.8       6) intel/17.0.4   8) impi/17.0.3      10) fftw3/3.3.6

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEELPA
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG -DVERBOSE

TACC_VEC_FLAGS = -xCORE-AVX512 # for skx only
FCPP    = cpp -C -nostdinc
F90free = mpiifort -free
LINK    = mpiifort -static-intel -static-libgcc -static-libstdc++  -Wl,--as-needed -liomp5 -Wl,--no-as-needed
# We need the -fp-model precise to pass the testsuite.

FOPTS   = -O3 $(TACC_VEC_FLAGS) -no-ipo -ip -align array64byte -threads -heap-arrays 4096 -fp-model fast=2 -complex-limited-range -assume byterecl -qopenmp -no-prec-div
#FOPTS   = -g -O0 -check all -Warn all -traceback $(TACC_VEC_FLAGS) -align array64byte -threads -heap-arrays 4096 -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -assume byterecl -qopenmp
FNOOPTS = $(FOPTS)

MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiicc
C_COMP  = mpiicc
C_LINK  = mpiicc -Wl,--as-needed -liomp5 -Wl,--no-as-needed
C_OPTS  = -O3 $(TACC_VEC_FLAGS) -no-ipo -ip -fp-model fast=2 -complex-limited-range -fno-alias -ansi-alias -qopenmp
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
MKLPATH      = $(MKLROOT)/lib/intel64
FFTWINCLUDE  = $(MKLROOT)/include/fftw

FFTWLIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lm -ldl

LAPACKLIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

SCALAPACKLIB = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a

HDF5PATH     = $(TACC_HDF5_LIB)
HDF5INCLUDE  = $(HDF5PATH)/../include
HDF5LIB      =  -Wl,--start-group $(HDF5PATH)/libhdf5hl_fortran.a $(HDF5PATH)/libhdf5_hl.a $(HDF5PATH)/libhdf5_fortran.a $(HDF5PATH)/libhdf5.a $(HDF5PATH)/libsz.a -Wl,--end-group -lz

ELPAINCLUDE = $(ELPA_DIR)/include/elpa/modules -I$(ELPA_DIR)/elpa/elpa
ELPALIB = -Wl,--start-group $(ELPA_DIR)/lib/libelpa_openmp.a -Wl,--end-group

SLATECLIB = $(SLATEC_DIR)/lib/libslatec.a

TESTSCRIPT = sbatch stampede2.scr
