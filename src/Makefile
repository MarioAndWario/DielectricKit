PREFIX=..
include $(PREFIX)/Common/common-rules.mk

# object files from $(COMMON) directory
COMOBJ = misc.o check_inversion.o inversion.o vcoul_generator.o \
          symmetries.o input_utils.o sort.o blas.o \
          io_utils.o epsread_hdf5.o hdf5_io.o so32su2.o wfn_io_hdf5.o epswrite_hdf5.o write_matrix.o scalapack.o


ifeq ($(findstring -DUSESCALAPACK,$(MATHFLAG)),-DUSESCALAPACK)

endif

clean-flavored:
	$(MAKE_CLEAN) *.mod

COMMONOBJS = $(addprefix $(COMMON)/,$(COMOBJ)) $(GLOBALOBJS)

SRC_EPSINV = EpsInv.f90
OBJ_EPSINV = $(SRC_EPSINV:.f90=.o)

# External libs
LIBS = $(SCALAPACKLIB) $(LAPACKLIB) $(SLATECLIB) $(HDF5LIB)

default: EpsInv
EpsInv: EpsInv$(FLAVOR).x
EpsInv$(FLAVOR).x: $(COMMONOBJS) $(SPGLIB)/libsymspg.a $(OBJ_EPSINV) $(COMMON)/scalapack.o $(COMMON)/lapack.o $(SPGLIB)/libsymspg.a
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)
