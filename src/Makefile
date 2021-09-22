PREFIX=..
include $(PREFIX)/Common/common-rules.mk
# include ../arch.mk
# include common-rules.mk

# object files from $(COMMON) directory
COMOBJ = misc.o check_inversion.o inversion.o vcoul_generator.o \
          symmetries.o input_utils.o sort.o blas.o \
          io_utils.o epsread_hdf5.o hdf5_io.o so32su2.o \
          wfn_io_hdf5.o epswrite_hdf5.o write_matrix.o scalapack.o

GLOBOBJ = global.o typedefs.o nrtype.o push_pop.o message.o peinfo.o timing.o intrinsics.o scalapack_aux.o

clean-flavored:
	$(MAKE_CLEAN) *.mod

GLOBALOBJS = $(addprefix ,$(GLOBOBJ))
GLOBALMODS = $(GLOBALOBJS:.o=_m.mod)
COMMONOBJS = $(addprefix ,$(COMOBJ)) $(GLOBALOBJS)

SRC_EPSINV = EpsInv.f90
OBJ_EPSINV = $(SRC_EPSINV:.f90=.o)

# External libs
LIBS = $(SCALAPACKLIB) $(LAPACKLIB) $(HDF5LIB) # $(SLATECLIB) 

all: EpsInv
default: EpsInv
EpsInv: EpsInv$(FLAVOR).x
EpsInv$(FLAVOR).x: $(COMMONOBJS) $(SPGLIB)/libsymspg.a $(OBJ_EPSINV) scalapack.o lapack.o $(SPGLIB)/libsymspg.a
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)

f_defs.h :
#$(ALL_COMMON_OBJ) : $(GLOBALMODS) f_defs.h compiler.h
global.o global_m.mod: scalapack_aux_m.mod nrtype_m.mod timing_m.mod typedefs_m.mod peinfo_m.mod push_pop_m.mod message_m.mod f_defs.h #compiler.h
misc.o misc_m.mod : blas_m.mod global_m.mod f_defs.h scalapack_m.mod
peinfo.o peinfo_m.mod : nrtype_m.mod f_defs.h intrinsics_m.mod
scalapack_aux.o scalapack_aux_m.mod : push_pop_m.mod f_defs.h
timing.o timing_m.mod : push_pop_m.mod nrtype_m.mod f_defs.h intrinsics_m.mod peinfo_m.mod
input_utils.o input_utils_m.mod : blas_m.mod global_m.mod f_defs.h
blas.o blas_m.mod : global_m.mod f_defs.h
intrinsics.o intrinsics_m.mod : f_defs.h #compiler.h
lapack.o lapack_m.mod : global_m.mod f_defs.h
scalapack.o scalapack_m.mod : global_m.mod f_defs.h
inversion.o inversion_m.mod: lapack_m.mod scalapack_m.mod #undef.h
sort.o sort_m.mod : global_m.mod f_defs.h
typedefs.o typedefs_m.mod : nrtype_m.mod f_defs.h
nrtype.o nrtype_m.mod : f_defs.h
push_pop.o push_pop_m.mod : peinfo_m.mod message_m.mod nrtype_m.mod f_defs.h
message.o message_m.mod : peinfo_m.mod nrtype_m.mod f_defs.h
check_inversion.o check_inversion_m.mod : global_m.mod
io_utils.o io_utils_m.mod : global_m.mod
write_matrix.o : scalapack_m.mod
write_matrix.o write_matrix_m.mod : global_m.mod scalapack_m.mod io_utils_m.mod f_defs.h hdf5_io_m.mod
symmetries.o symmetries_m.mod : global_m.mod misc_m.mod sort_m.mod
hdf5_io.o hdf5_io_m.mod : global_m.mod
wfn_io_hdf5.o wfn_io_hdf5_m.mod : hdf5_io_m.mod global_m.mod
epsread_hdf5.o epsread_hdf5_m.mod : global_m.mod hdf5_io_m.mod
epswrite_hdf5.o epswrite_hdf5_m.mod : global_m.mod hdf5_io_m.mod wfn_io_hdf5_m.mod
so32su2.o so32su2_m.mod : so32su2.f90 global_m.mod f_defs.h
