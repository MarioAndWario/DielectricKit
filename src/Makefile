include ../arch.mk

ifeq ($(findstring -DUSESCALAPACK,$(MATHFLAG)),-DUSESCALAPACK)
  ifneq ($(findstring -DMPI,$(PARAFLAG)),-DMPI)
    $(error -DUSESCALAPACK flag requires -DMPI; please check your arch.mk file)
  endif
endif

# Object files for module files
COMMONOBJS = misc.o check_inversion.o inversion.o vcoul_generator.o \
          symmetries.o input_utils.o sort.o blas.o \
          io_utils.o epsread_hdf5.o hdf5_io.o so32su2.o \
          wfn_io_hdf5.o epswrite_hdf5.o write_matrix.o scalapack.o lapack.o

# Global object files for all subroutines
GLOBALOBJS = global.o intrinsics.o nrtype.o typedefs.o peinfo.o message.o push_pop.o scalapack_aux.o timing.o 
GLOBALMODS = $(GLOBALOBJS:.o=_m.mod)

# External libs
LIBS = $(SCALAPACKLIB) $(LAPACKLIB) $(HDF5LIB) # $(SLATECLIB) 
SPGLIB = ../spglib-1.0.9

#Default values for the C++ library, stored in $(CC_LINK_LIBS), when linking
#C++ code with Fortran linker. NOTE: we don`t have a default value for xlf!
ifeq ($(CC_LINK_LIBS),)
  ifeq ($(findstring -DINTEL,$(COMPFLAG)),-DINTEL)
    CC_LINK_LIBS = -cxxlib
  else ifeq ($(findstring -DPGI,$(COMPFLAG)),-DPGI)
    CC_LINK_LIBS = -pgcpplibs
  else ifeq ($(findstring -DPATH,$(COMPFLAG)),-DPATH)
    CC_LINK_LIBS = -lstl -lcxxrt -lpthread -ldl -lgcc -leh -lc -lmv -lmpath -lm
  else ifeq ($(findstring -DSUN,$(COMPFLAG)),-DSUN)
    CC_LINK_LIBS = -lCstd -lCrun -lc
  else
#Default case: works with GNU,G95,CRAY,OPEN64 plus many combinations of
# compilers, such as NAG+g++, ABSOFT+g++, etc.
    CC_LINK_LIBS = -lstdc++
  endif
endif

# Fortran preprocessing
CPPOPT  = $(TYPEFLAG) $(COMPFLAG) $(PARAFLAG) $(MATHFLAG) $(DEBUGFLAG)

# C/C++
C_CPPOPT  = $(TYPEFLAG) $(COMPFLAG) $(C_PARAFLAG) $(C_DEBUGFLAG)

INCLUDE = -I.
FTNINC += -I. -I$(HDF5INCLUDE)
f90_CPP = $(FCPP) $(INCLUDE) $(CPPOPT) $< > $(basename $<).p.f
F90_CMD = $(F90free) $(FTNINC) -c $(FOPTS) $(basename $<).p.f -o $(basename $<).o $(MOD_OPT)$(dir $<)
RM_P_F_CMD = @$(REMOVE) $(basename $<).p.f

# Rules:
# Remove all implicit suffix rules
.SUFFIXES:

%.o %_m.mod : %.f90
	$(f90_CPP)
	$(F90_CMD)
	$(RM_P_F_CMD)

%.o : %.cpp
	$(CC_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@

%.o : %.cc
	$(CC_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@

%.o : %.c
	$(C_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@

default: EpsInv
EpsInv: EpsInv.x
EpsInv.x: $(GLOBALOBJS) $(COMMONOBJS) $(SPGLIB)/libsymspg.a EpsInv.o 
	$(LINK) $(FOPTS) -o $@ $^ $(LIBS)

# Dependence
$(COMMONOBJS) : $(GLOBALMODS) f_defs.h

# Basic modules
typedefs.o typedefs_m.mod : nrtype_m.mod
peinfo.o peinfo_m.mod : nrtype_m.mod intrinsics_m.mod 
message.o message_m.mod : peinfo_m.mod nrtype_m.mod
push_pop.o push_pop_m.mod : peinfo_m.mod message_m.mod nrtype_m.mod
scalapack_aux.o scalapack_aux_m.mod : push_pop_m.mod
timing.o timing_m.mod : push_pop_m.mod nrtype_m.mod intrinsics_m.mod peinfo_m.mod
global.o global_m.mod: intrinsics_m.mod scalapack_aux_m.mod nrtype_m.mod timing_m.mod typedefs_m.mod peinfo_m.mod push_pop_m.mod message_m.mod

# Utilitarian modules
input_utils.o input_utils_m.mod : blas_m.mod
misc.o misc_m.mod : blas_m.mod scalapack_m.mod
symmetries.o symmetries_m.mod : misc_m.mod sort_m.mod
inversion.o inversion_m.mod: lapack_m.mod scalapack_m.mod
write_matrix.o write_matrix_m.mod : scalapack_m.mod io_utils_m.mod hdf5_io_m.mod
wfn_io_hdf5.o wfn_io_hdf5_m.mod : hdf5_io_m.mod
epsread_hdf5.o epsread_hdf5_m.mod : hdf5_io_m.mod
epswrite_hdf5.o epswrite_hdf5_m.mod : hdf5_io_m.mod wfn_io_hdf5_m.mod

# Main object
EpsInv.o: $(GLOBALOBJS) $(COMMONOBJS) $(SPGLIB)/libsymspg.a 

# Space group library
SPGLIB_SRC = $(addprefix $(SPGLIB)/, cell.c debug.c hall_symbol.c kpoint.c lattice.c mathfunc.c pointgroup.c primitive.c \
refinement.c site_symmetry.c sitesym_database.c spacegroup.c spg_database.c \
spglib.c symmetry.c spglib_f.c spglib_f_meta.c)
SPGLIB_OBJ = $(SPGLIB_SRC:.c=.o)

SPGLIB_HEADERS = $(addprefix $(SPGLIB)/, cell.h debug.h hall_symbol.h kpoint.h lattice.h mathfunc.h \
pointgroup.h primitive.h refinement.h site_symmetry.h sitesym_database.h \
spacegroup.h spg_database.h spglib.h symmetry.h)

spglib: $(SPGLIB)/libsymspg.a
clean-spglib:
	cd $(SPGLIB) && $(MAKE) clean
cleanall-spglib:
	cd $(SPGLIB) && $(MAKE) cleanall

SPGLIB_OBJ : $(SPGLIB_HEADERS)

ifndef AR
AR = /usr/bin/ar
endif

$(SPGLIB)/libsymspg.a : $(SPGLIB_OBJ)
	$(AR) ru $@ $^

clean:
	rm -f *.o *.p.f *~ core *__genmod.f90 *.a *.mod *.x

.PHONY: all default clean spglib clean-spglib cleanall-spglib
