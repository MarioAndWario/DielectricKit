# arch.mk for BerkeleyGW codes using Intel compiler.
# Features: MPI, OpenMP, LibSci, HDF5
#
# Suitable for Cori Phase 1 at NERSC
#
 #  1) modules/3.2.11.4                                 13) dvs/2.12_2.2.156-7.0.1.1_8.9__g5aab709e
 #  2) altd/2.0                                         14) alps/6.6.58-7.0.1.1_6.4__g437d88db.ari
 #  3) craype-network-aries                             15) rca/2.2.20-7.0.1.1_4.46__g8e3fb5b.ari
 #  4) craype/2.6.2                                     16) atp/2.1.3
 #  5) cray-libsci/19.06.1                              17) PrgEnv-intel/6.0.5
 #  6) udreg/2.3.2-7.0.1.1_3.31__g8175d3d.ari           18) craype-haswell
 #  7) ugni/6.0.14.0-7.0.1.1_7.33__ge78e5b0.ari         19) cray-mpich/7.7.10
 #  8) pmi/5.0.14                                       20) craype-hugepages2M
 #  9) dmapp/7.1.1-7.0.1.1_4.48__g38cf134.ari           21) intel/19.0.3.199
 # 10) gni-headers/5.0.12.0-7.0.1.1_6.28__g3b1768f.ari  22) cray-hdf5-parallel/1.10.5.2
 # 11) xpmem/2.2.20-7.0.1.1_4.10__g0475745.ari          23) cray-fftw/3.3.8.4
 # 12) job/2.2.4-7.0.1.1_3.36__g36b56f4.ari             24) python/3.7-anaconda-2019.10

# Precompiler options
# export CRAYPE_LINK_TYPE=static
COMPFLAG  = -DINTEL
PARAFLAG  = -DMPI -DOMP
MATHFLAG  = -DUSESCALAPACK -DUNPACKED -DUSEFFTW3 -DHDF5 -DUSEELPA
#-DSHEPPACK
# Only uncomment DEBUGFLAG if you need to develop/debug BerkeleyGW.
# The output will be much more verbose, and the code will slow down by ~20%.
#DEBUGFLAG = -DDEBUG -DVERBOSE

FCPP    = /usr/bin/cpp -C -nostdinc
F90free = ftn -free
LINK    = ftn -static-intel -static-libgcc -static-libstdc++  -Wl,--as-needed -liomp5 -Wl,--no-as-needed
FOPTS   = -O3 -xCORE-AVX2 -no-ipo -ip -align array64byte -heap-arrays 4096 -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -assume byterecl -qopenmp
#FOPTS   =  -g -O0 -check none -Warn all -traceback -xCORE-AVX2 -no-ipo -ip -align array64byte -heap-arrays 4096 -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -assume byterecl -qopenmp
FNOOPTS = $(FOPTS)
MOD_OPT = -module 
INCFLAG = -I

C_PARAFLAG  = -DPARA -DMPICH_IGNORE_CXX_SEEK
CC_COMP = CC
C_COMP  = cc
C_LINK  = CC -Wl,--as-needed -liomp5 -Wl,--no-as-needed
C_OPTS = -O3 -xCORE-AVX2 -no-ipo -ip -qopt-zmm-usage=high -fp-model fast=2 -complex-limited-range -fno-alias -ansi-alias -qopenmp
C_DEBUGFLAG =

REMOVE  = /bin/rm -f

LAPACKLIB = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl

FFTWLIB = 
FFTWINCLUDE  = $(MKLROOT)/include/fftw/
PERFORMANCE  = 

HDF5_LDIR    =  $(HDF5_DIR)/lib
HDF5LIB      =  -L$(HDF5_LDIR) -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -ldl

HDF5INCLUDE  = $(HDF5_DIR)/include

ELPAINCLUDE = ${ELPA_DIR}/include/elpa_openmp-2017.11.001/modules -I${ELPA_DIR}/include/elpa_openmp-2017.11.001/elpa
ELPALIB = ${ELPA_DIR}/lib/libelpa_openmp.a

# SHEPPACK_DIR = /global/homes/m/mwu/softwares/SHEPPACK/bulid_BGW
# SHEPPACKINCLUDE = /global/homes/m/mwu/softwares/SHEPPACK/bulid_BGW

SLATECLIB = $(SLATEC_DIR)/lib/libslatec.a

TESTSCRIPT = sbatch cori2.scr
