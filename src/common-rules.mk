#include ../arch.mk

# it is ok not to be present for some targets such as 'clean' and 'all-flavors'
# if it is needed, it will cause an error according to the flavor.mk target below.

# ifeq ($(findstring -DUSESCALAPACK,$(MATHFLAG)),-DUSESCALAPACK)
#   ifneq ($(findstring -DMPI,$(PARAFLAG)),-DMPI)
#     $(error -DUSESCALAPACK flag requires -DMPI; please check your arch.mk file)
#   endif
# endif

# #Default values for the C++ library, stored in $(CC_LINK_LIBS), when linking
# #C++ code with Fortran linker. NOTE: we don`t have a default value for xlf!
# ifeq ($(CC_LINK_LIBS),)
#   ifeq ($(findstring -DINTEL,$(COMPFLAG)),-DINTEL)
#     CC_LINK_LIBS = -cxxlib
#   else ifeq ($(findstring -DPGI,$(COMPFLAG)),-DPGI)
#     CC_LINK_LIBS = -pgcpplibs
#   else ifeq ($(findstring -DPATH,$(COMPFLAG)),-DPATH)
#     CC_LINK_LIBS = -lstl -lcxxrt -lpthread -ldl -lgcc -leh -lc -lmv -lmpath -lm
#   else ifeq ($(findstring -DSUN,$(COMPFLAG)),-DSUN)
#     CC_LINK_LIBS = -lCstd -lCrun -lc
#   else
# #Default case: works with GNU,G95,CRAY,OPEN64 plus many combinations of
# # compilers, such as NAG+g++, ABSOFT+g++, etc.
#     CC_LINK_LIBS = -lstdc++
#   endif
# endif

# ifeq ($(findstring USESCALAPACK,$(MATHFLAG)),)
#   SCALAPACKLIB =
# endif

# # Fortran preprocessing
# CPPOPT  = $(TYPEFLAG) $(COMPFLAG) $(PARAFLAG) $(MATHFLAG) $(DEBUGFLAG)

# # C/C++
# C_CPPOPT  = $(TYPEFLAG) $(COMPFLAG) $(C_PARAFLAG) $(C_DEBUGFLAG)

#COMMON = .
#COMMON = $(PREFIX)/Common
#SPGLIB = ..//spglib-1.0.9

# # this one is for C, C++, and cpp step for Fortran
# INCLUDE = -I$(COMMON)
# # this one is for Fortran compilation
# FTNINC = $(INCFLAG) $(COMMON) $(INCFLAG)
# ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5)
#   FTNINC += $(INCFLAG) $(HDF5INCLUDE)
# endif
# ifeq ($(findstring -DUSEELPA,$(MATHFLAG)),-DUSEELPA)
#   FTNINC += $(INCFLAG) $(ELPAINCLUDE)
# endif
# ifeq ($(findstring -DSHEPPACK,$(MATHFLAG)),-DSHEPPACK)
#   FTNINC += $(INCFLAG) $(SHEPPACKINCLUDE)
# endif

#MAKE_CLEAN = -$(REMOVE) *.o *.p.f *~ core *__genmod.f90 *.a

# ifort -warn all creates __genmod.f90 files as pseudo-modules
#F90_CMD = $(F90free) $(FTNINC) -c $(FOPTS) $(basename $<).f90 -o $(basename $<).o $(MOD_OPT)$(dir $<)
#F90_CMD = $(F90free) $(FTNINC) -c $(FOPTS) $(basename $<).p.f -o $(basename $<).o $(MOD_OPT)$(dir $<)

# #clang C-preprocessing treats files incorrectly if they have .F90 extension
# ifeq ($(findstring clang,$(FCPP)),clang)
#   f90_CPP = cp $(basename $<).f90 $(basename $<)_cp.F90; $(FCPP) $(INCLUDE) $(CPPOPT) $(basename $<)_cp.F90 > $(basename $<).p.f; $(REMOVE) $(basename $<)_cp.F90
# else
#   f90_CPP = $(FCPP) $(INCLUDE) $(CPPOPT) $< > $(basename $<).p.f
# endif
# F90_CPP = $(FCPP) -P $(INCLUDE) $(CPPOPT) $< > $(basename $<).90

# ifneq (,$(filter $(COMPFLAG),-DOPEN64 -DPATH -DABSOFT -DCRAY))
# # these compilers name all modules uppercase
# MODLINK = @NAME=$(dir $*)`echo $(notdir $*) | tr '[:lower:]' '[:upper:]'`_M.mod; test ! -e $$NAME || ln -sf $(PWD)/$$NAME $(PWD)/$(basename $<)_m.mod
# endif

#default_goal: default
# GNU make 3.80 and earlier don't have .DEFAULT_GOAL

# # all objects to be made from Common directory must appear here
# ALL_COMOBJ = scalapack.o inversion.o write_matrix.o vcoul_generator.o check_inversion.o  \
#       sort.o blas.o scalapack.o lapack.o misc.o input_utils.o \
#       symmetries.o hdf5_io.o wfn_io_hdf5.o epsread_hdf5.o epswrite_hdf5.o io_utils.o so32su2.o

# ALL_COMMON_OBJ = $(addprefix ,$(ALL_COMOBJ))

# # these are involved in the modules which all routines must use
# GLOBOBJ = global.o typedefs.o nrtype.o push_pop.o message.o peinfo.o timing.o intrinsics.o scalapack_aux.o
# GLOBALOBJS = $(addprefix ,$(GLOBOBJ))
# GLOBALMODS = $(GLOBALOBJS:.o=_m.mod)

# If you add any module in the Common directory, add the corresponding .o and _m.mod in a rule here, with some dependency
# e.g. example.o example_m.mod : global_m.mod
# f_defs.h : 
# $(ALL_COMMON_OBJ) : $(GLOBALMODS) f_defs.h compiler.h
# global.o global_m.mod: scalapack_aux_m.mod nrtype_m.mod timing_m.mod typedefs_m.mod \
# peinfo_m.mod push_pop_m.mod message_m.mod f_defs.h compiler.h
# misc.o misc_m.mod : blas_m.mod global_m.mod f_defs.h scalapack_m.mod
# peinfo.o peinfo_m.mod : nrtype_m.mod f_defs.h intrinsics_m.mod
# scalapack_aux.o scalapack_aux_m.mod : push_pop_m.mod f_defs.h
# timing.o timing_m.mod : push_pop_m.mod nrtype_m.mod f_defs.h intrinsics_m.mod peinfo_m.mod
# input_utils.o input_utils_m.mod : blas_m.mod global_m.mod f_defs.h
# blas.o blas_m.mod : global_m.mod f_defs.h
# intrinsics.o intrinsics_m.mod : f_defs.h compiler.h
# lapack.o lapack_m.mod : global_m.mod f_defs.h
# scalapack.o scalapack_m.mod : global_m.mod f_defs.h
# inversion.o inversion_m.mod: lapack_m.mod scalapack_m.mod undef.h
# sort.o sort_m.mod : global_m.mod f_defs.h
# typedefs.o typedefs_m.mod : nrtype_m.mod f_defs.h
# nrtype.o nrtype_m.mod : f_defs.h
# push_pop.o push_pop_m.mod : peinfo_m.mod message_m.mod nrtype_m.mod f_defs.h
# message.o message_m.mod : peinfo_m.mod nrtype_m.mod f_defs.h
# check_inversion.o check_inversion_m.mod : global_m.mod
# io_utils.o io_utils_m.mod : global_m.mod
# write_matrix.o : scalapack_m.mod
# write_matrix.o write_matrix_m.mod : global_m.mod scalapack_m.mod io_utils_m.mod f_defs.h hdf5_io_m.mod
# symmetries.o symmetries_m.mod : global_m.mod misc_m.mod sort_m.mod
# hdf5_io.o hdf5_io_m.mod : global_m.mod
# wfn_io_hdf5.o wfn_io_hdf5_m.mod : hdf5_io_m.mod global_m.mod
# epsread_hdf5.o epsread_hdf5_m.mod : global_m.mod hdf5_io_m.mod
# epswrite_hdf5.o epswrite_hdf5_m.mod : global_m.mod hdf5_io_m.mod wfn_io_hdf5_m.mod
# so32su2.o so32su2_m.mod : so32su2.f90 global_m.mod f_defs.h

#ifeq ($(COMPFLAG),-DOPEN64)
# other compilers provide this module, so only Open64 needs its own here
# To use it, copy in the file http://www.open64.net/doc/db/d6b/omp__lib_8f-source.html to this directory
# and uncomment the dependencies below,  needed to deal with the fact that the module name doesn't end on _m.
# intrinsics_m.mod : OMP_LIB.mod 
# OMP_LIB.mod : omp_lib.o 
#endif

#common: $(ALL_COMMON_OBJ) $(GLOBALOBJS) 

# spglib: $(SPGLIB)/libsymspg.a
# clean-spglib:
# 	cd $(SPGLIB) && $(MAKE) clean
# cleanall-spglib:
# 	cd $(SPGLIB) && $(MAKE) cleanall

# SPGLIB_SRC = $(addprefix $(SPGLIB)/, cell.c debug.c hall_symbol.c kpoint.c lattice.c mathfunc.c pointgroup.c primitive.c \
# refinement.c site_symmetry.c sitesym_database.c spacegroup.c spg_database.c \
# spglib.c symmetry.c spglib_f.c spglib_f_meta.c)
# SPGLIB_OBJ = $(SPGLIB_SRC:.c=.o)

# SPGLIB_HEADERS = $(addprefix $(SPGLIB)/, cell.h debug.h hall_symbol.h kpoint.h lattice.h mathfunc.h \
# pointgroup.h primitive.h refinement.h site_symmetry.h sitesym_database.h \
# spacegroup.h spg_database.h spglib.h symmetry.h)

# SPGLIB_OBJ : $(SPGLIB_HEADERS)

# ifndef AR
# AR = /usr/bin/ar
# endif
# $(SPGLIB)/libsymspg.a : $(SPGLIB_OBJ)
# 	$(AR) ru $@ $^

#Rules

#.SUFFIXES:
# remove all implicit suffix rules

# #GNU and gfortran don't update mod files unless interfaces have changed, frustrating make's dependency-checking
# ifneq ($(findstring -DG,$(COMPFLAG)),)
# RM_MOD_CMD = @$(REMOVE) $(basename $<)_m.mod
# endif

# ifeq ($(findstring -P,$(FCPP)),)
# # if not running with cpp -P, keep .p.f
# ######
# RM_P_F_CMD = @$(REMOVE) $(basename $<).p.f
# ######
# endif

# # both files are made at once by this rule
# %.o %_m.mod : %.f90
# 	$(RM_MOD_CMD)
# 	$(f90_CPP)
# 	$(F90_CMD)
# 	$(RM_P_F_CMD)
# 	$(MODLINK)

# # Files including _inc.f90 by cpp have their own rules to keep .p.f, and not insert line markers,
# # so you can tell which of the various preprocessed versions caused any runtime error.
# %.o %_m.mod : %.F90 %_inc.f90
# 	$(RM_MOD_CMD)
# 	$(F90_CPP)
# 	$(F90_CMD)
# 	$(MODLINK)

# # FHJ: Use this rule if you want to generate an assembly file.
# # You should only do this if you care about fine optimization.
# %.s : %.p.f
# 	$(F90_ASM_CMD)

# rules to make .p.f by hand for debugging purposes
# %.p.f : %.f90
# 	$(f90_CPP)
# %.p.f : %.F90 %_inc.f90
# 	$(F90_CPP)

# %.o : %.cpp
# 	$(CC_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@

# %.o : %.cc
# 	$(CC_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@

# Alas, it seems that if you compile a .c with gcc or .cpp with g++ you get a leading underscore in the symbol
# table, but it you compile a .c with g++ you do not. Therefore we must make a distinction between C and C++ (for spglib).
# %.o : %.c
# 	$(C_COMP) $(INCLUDE) $(C_CPPOPT) -c $(C_OPTS) $< -o $@

# clean:
# 	$(MAKE_CLEAN) *.mod *.x


# all targets which are not the name of a file should be listed here
#.PHONY: default all all-j
