PREFIX=..
include $(PREFIX)/Common/common-rules.mk

# object files from $(COMMON) directory
COMOBJ = misc.o gmap.o check_inversion.o inversion.o vcoul_generator.o fftw.o fft_parallel.o \
         minibzaverage.o symmetries.o fullbz.o createpools.o input_utils.o \
         wfn_rho_vxc_io.o sort.o blas.o cells.o \
         io_utils.o epsread_hdf5.o hdf5_io.o so32su2.o wfn_io_hdf5.o epswrite_hdf5.o write_matrix.o irrbz.o \
         complex_incomplete_gamma.o scalapack.o

ifeq ($(findstring -DUSESCALAPACK,$(MATHFLAG)),-DUSESCALAPACK)

endif

clean-flavored:
	$(MAKE_CLEAN) *.mod

COMMONOBJS = $(addprefix $(COMMON)/,$(COMOBJ)) $(GLOBALOBJS) $(TILE_LIBS)

SRC_CM	= inread_calc_mdat.f90 distrib_cm.f90 input_cm.f90 mtxel_cc.f90 genwf_bse.f90 genwf_kernel.f90 calc_mdat.f90 
OBJ_CM	= $(SRC_CM:.f90=.o)

SRC_CC	= inread_calc_chi.f90 distrib_cc.f90 input_cc.f90 genwf_bse.f90 genwf_cc.f90 mtxel_cc.f90 calc_chi.f90
OBJ_CC	= $(SRC_CC:.f90=.o)

SRC_EPSINV = EpsInv.f90
OBJ_EPSINV = $(SRC_EPSINV:.f90=.o)

TOOLOBJ = eigvec_k.o

# put truly external libs here, i.e. not part of our buildsystem at all
LIBS = $(SCALAPACKLIB) $(FFTWLIB) $(LAPACKLIB) $(SLATECLIB) $(HDF5LIB) $(TILE_LIBS_EXT)

default: EpsInv
all: calc_mdat calc_chi EpsInv eigvec_k

calc_mdat: calc_mdat$(FLAVOR).x
calc_chi: calc_chi$(FLAVOR).x
EpsInv: EpsInv$(FLAVOR).x
eigvec_k: eigvec_k$(FLAVOR).x
tools: eigvec_k$(FLAVOR).x

calc_mdat$(FLAVOR).x: $(COMMONOBJS) $(OBJ_CM) $(COMMON)/scalapack.o $(COMMON)/lapack.o $(SPGLIB)/libsymspg.a
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)

calc_chi$(FLAVOR).x: $(COMMONOBJS) $(OBJ_CC) $(COMMON)/scalapack.o $(COMMON)/lapack.o $(SPGLIB)/libsymspg.a
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)

EpsInv$(FLAVOR).x: $(COMMONOBJS) $(SPGLIB)/libsymspg.a $(OBJ_EPSINV) $(COMMON)/scalapack.o $(COMMON)/lapack.o $(SPGLIB)/libsymspg.a
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)

eigvec_k$(FLAVOR).x: eigvec_k.o $(GLOBALOBJS)
	$(LINK) $(FOPTS) -o $@ $^

$(OBJ_K) $(OBJ_ABS) $(OBJ_INT) $(TOOLOBJ): $(GLOBALMODS)
input_cm.o :  $(COMMON)/sort_m.mod $(COMMON)/misc_m.mod $(COMMON)/input_utils_m.mod $(COMMON)/wfn_rho_vxc_io_m.mod $(COMMON)/read_rho_vxc_m.mod $(COMMON)/fullbz_m.mod $(COMMON)/checkbz_m.mod
input_cc.o :  $(COMMON)/sort_m.mod $(COMMON)/misc_m.mod $(COMMON)/input_utils_m.mod $(COMMON)/wfn_rho_vxc_io_m.mod $(COMMON)/read_rho_vxc_m.mod $(COMMON)/fullbz_m.mod $(COMMON)/checkbz_m.mod distrib_cc.o
distrib_cm.o distrib_cc.o : $(COMMON)/misc_m.mod
genwf_cc.o : genwf_cc.f90 $(COMMON)/blas_m.mod $(COMMON)/gmap_m.mod $(COMMON)/input_utils_m.mod $(COMMON)/misc_m.mod $(COMMON)/so32su2_m.mod $(COMMON)/sort_m.mod
mtxel_cc.o mtxel_cc_m.mod : $(COMMON)/fftw_m.mod $(COMMON)/misc_m.mod 
genwf_kernel.o : $(COMMON)/blas_m.mod $(COMMON)/gmap_m.mod $(COMMON)/input_utils_m.mod $(COMMON)/misc_m.mod $(COMMON)/so32su2_m.mod $(COMMON)/sort_m.mod
inread_calc_mdat.o inread_calc_mdat_m.mod : $(GLOBALMODS)
inread_calc_chi.o inread_calc_chi_m.mod : $(GLOBALMODS)
