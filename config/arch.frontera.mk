# arch.mk for BerkeleyGW codes
#
# suitable for Frontera at TACC
# NOTE: epsilon sometimes freezes when compiled with HDF5.
# We are working on this issue.
#
  # 1) intel/19.0.5   3) autotools/1.2   5) pmix/3.1.4      7) xalt/2.8.5   9) impi/19.0.7    11) phdf5/1.8.16
  # 2) git/2.24.1     4) cmake/3.16.1    6) hwloc/1.11.12   8) TACC        10) python3/3.7.0  12) fftw3/3.3.8

COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEELPA
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG #-DVERBOSE

TACC_VEC_FLAGS = -xCORE-AVX512 # for skx only
FCPP    = cpp -C -nostdinc
F90free = mpiifort -free
LINK    = mpiifort -static-intel -static-libgcc -static-libstdc++  -Wl,--as-needed -liomp5 -Wl,--no-as-needed
FOPTS   = -O3 $(TACC_VEC_FLAGS) -no-ipo -ip -align array64byte -fp-model fast=2 -complex-limited-range -assume byterecl -assume nobuffered_io -assume nobuffered_stdout  -no-prec-div -qopenmp -qopt-zmm-usage=high -threads -heap-arrays 4096
#FOPTS   = -g -O0 -check all -Warn all -traceback -qopenmp -assume nobuffered_io -assume nobuffered_stdout
#https://stackoverflow.com/questions/9751996/how-can-i-find-the-cause-for-a-memory-leak-in-fortran-2003-program?answertab=oldest#tab-top
#FOPTS = -O2 -qopenmp -stand f03 -assume realloc_lhs  -check all  -traceback  -warn all  -fstack-protector  -assume protect_parens  -implicitnone
FNOOPTS = $(FOPTS)

MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = mpiicc
C_COMP  = mpiicc
C_LINK  = mpiicc -Wl,--as-needed -liomp5 -Wl,--no-as-needed
C_OPTS  = -O3 $(TACC_VEC_FLAGS) -no-ipo -ip -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -fno-alias -ansi-alias -qopenmp
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

# Math Libraries
MKLPATH      = $(MKLROOT)/lib/intel64
FFTWINCLUDE  = $(MKLROOT)/include/fftw

FFTWLIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lm -ldl

LAPACKLIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group

SCALAPACKLIB = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a

HDF5PATH     = $(TACC_HDF5_LIB)
HDF5INCLUDE  = $(HDF5PATH)/../include
HDF5LIB      =  -Wl,--start-group $(HDF5PATH)/libhdf5hl_fortran.a $(HDF5PATH)/libhdf5_hl.a $(HDF5PATH)/libhdf5_fortran.a $(HDF5PATH)/libhdf5.a $(HDF5PATH)/libsz.a -Wl,--end-group -lz

ELPA_DIR = /home1/03355/tg826544/software/ELPA/elpa_2017/elpa/intel-skx-omp
ELPAINCLUDE = ${ELPA_DIR}/include/elpa/modules -I${ELPA_DIR}/elpa/elpa
ELPALIB = -Wl,--start-group ${ELPA_DIR}/lib/libelpa_openmp.a -Wl,--end-group

SLATECLIB = $(SLATEC_DIR)/lib/libslatec.a

TESTSCRIPT = sbatch frontera.scr
