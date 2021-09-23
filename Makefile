## Makefile (D. Prendergast Jul 2008)
##

PREFIX=.
include $(PREFIX)/Common/common-rules.mk

pre:
	$(MAKE) donkey && $(MAKE) common && $(MAKE) spglib

default_: epsilon sigma bse
default:
	$(MAKE) pre && $(MAKE) default_

list:
	@echo
	@echo "BerkeleyGW Makefile"
	@echo
	@echo "Available make targets:"
	@echo
	@echo "make default"
	@echo "make all"
	@echo "make all-flavors"
	@echo "make install"
	@echo "make check"
	@echo "make check-save"
ifdef TESTSCRIPT
	@echo "make check-jobscript"
	@echo "make check-jobscript-save"
endif
	@echo "make core-tools"
	@echo "make other-tools"
	@echo "make common"
	@echo "make library"
	@echo "make epsilon"
	@echo "make epsilon-all"
	@echo "make epsilon-tools"
	@echo "make sigma"
	@echo "make sigma-all"
	@echo "make sigma-tools"
	@echo "make bse"
	@echo "make bse-all"
	@echo "make bse-tools"
	@echo "make plotxct"
	@echo "make meanfield"
	@echo "make spglib"
	@echo "make sapo"
#BEGIN_INTERNAL_ONLY
	@echo "make nonlinearoptics"
	@echo "make cumexp"
#END_INTERNAL_ONLY
	@echo "make epm"
	@echo "make siesta2bgw"
	@echo "make bgw2para"
	@echo "make abi2bgw"
	@echo "make kgrid"
	@echo "make icm"
	@echo "make utilities"
	@echo "make surface"
	@echo "make qhull"
	@echo
	@echo "make clean"
	@echo "make clean-flavored"
	@echo "make clean-keepmod"
	@echo "make cleanall"
	@echo "make cleanall-bin"
	@echo "make clean-cpp"
	@echo

help: list

# All
#
all_: epsilon-all sigma-all bse-all other-tools library manual
all:
	$(MAKE) pre && $(MAKE) all_
all-j:
	@echo
	@echo 'WARNING: "all-j" target is deprecated.'
	@echo
	@sleep 2 && $(MAKE) all

all-flavors:
	$(MAKE) clean
	cp flavor_real.mk flavor.mk
	$(MAKE) all
	$(MAKE) clean-flavored
	cp flavor_cplx.mk flavor.mk
	$(MAKE) all

# Deliberately make this always be executed, by not mentioning the file manual.html
# It depends on so many files it would be tedious to list them all here.
manual: printsvninfo
	bin/assemble_manual.sh > manual.html

install: all
ifdef INSTDIR
	mkdir -p $(INSTDIR)/bin
	install bin/*.x $(INSTDIR)/bin/
	install bin/*.sh $(INSTDIR)/bin/
	install bin/*.py $(INSTDIR)/bin/
	install bin/*.pl $(INSTDIR)/bin/

	mkdir -p $(INSTDIR)/lib
	install library/*.a $(INSTDIR)/lib/
	mkdir -p $(INSTDIR)/include
	install library/*.mod $(INSTDIR)/include/

	mkdir -p $(INSTDIR)/share
	mkdir -p $(INSTDIR)/share/BerkeleyGW
# install cannot work on a whole directory
	cp -rf examples  $(INSTDIR)/share/BerkeleyGW/
	cp -rf testsuite $(INSTDIR)/share/BerkeleyGW/
	install manual.html $(INSTDIR)/share/BerkeleyGW/
else
	$(error Error: Please define installation prefix INSTDIR via 'make install INSTDIR='.)
endif

check: all
	cd testsuite && $(MAKE) check

check-save: all
	cd testsuite && $(MAKE) check-save

ifdef TESTSCRIPT
check-jobscript: all
	cd testsuite && $(MAKE) check-jobscript
check-jobscript-save: all
	cd testsuite && $(MAKE) check-jobscript-save
endif

core-tools: epsilon-tools sigma-tools bse-tools

other-tools: plotxct meanfield surface cumexp

# We use the shell AND operator (&&) so that if the cd command fails,
# the script will fail without trying to invoke the make command in
# the wrong directory, which could cause problems.

# Common
#
common:
	cd $(COMMON) && $(MAKE)
# Library
#
library:
	cd library && $(MAKE)

# Epsilon
#
epsilon:
	cd Epsilon && $(MAKE)
#
epsilon-all:
	cd Epsilon && $(MAKE) all
#
epsilon-tools:
	cd Epsilon && $(MAKE) tools

# Sigma
#
sigma:
	cd Sigma && $(MAKE)
#
sigma-all:
	cd Sigma && $(MAKE) all
#
sigma-tools:
	cd Sigma && $(MAKE) tools

# BSE
#
bse:
	cd BSE && $(MAKE)
#
bse-all:
	cd BSE && $(MAKE) all
#
bse-tools:
	cd BSE && $(MAKE) tools

# PlotXct
#
plotxct:
	cd PlotXct && $(MAKE)

# MeanField
#
meanfield:
	cd MeanField && $(MAKE) all

sapo:
	cd MeanField/SAPO && $(MAKE)

#BEGIN_INTERNAL_ONLY
other-tools: nonlinearoptics cumexp
nonlinearoptics:
	cd NonLinearOptics && $(MAKE)

# GW+C
#
cumexp:
	cd GW+C && $(MAKE)

# ABINIT2BGW
#
abi2bgw:
	cd MeanField/ABINIT && $(MAKE)
#END_INTERNAL_ONLY

# EPM
#
epm:
	cd MeanField/EPM && $(MAKE) all

# SIESTA2BGW
#
siesta2bgw:
	cd MeanField/SIESTA && $(MAKE)

# BGW2PARA
#
bgw2para:
	cd MeanField/PARATEC && $(MAKE)

# kgrid
#
kgrid:
	cd MeanField/ESPRESSO && $(MAKE)

# ICM
#
icm:
	cd MeanField/ICM && $(MAKE)

# Utilities
#
utilities:
	cd MeanField/Utilities && $(MAKE)

# Surface
#
surface:
	cd Visual && $(MAKE)

# Qhull targets
#
# defined in common-rules.mk

# Clean
#
clean: clean-common clean-epsilon clean-sigma clean-bse clean-plotxct clean-meanfield \
       clean-surface clean-qhull clean-manual
#
# FIXME: most of common is actually flavorless.
clean-flavored: clean-common clean-epsilon clean-sigma clean-bse-flavored clean-plotxct clean-meanfield-flavored
#
clean-common:
	cd $(COMMON) && $(MAKE) clean
#
clean-library:
	cd library && $(MAKE) clean
#
clean-epsilon:
	cd Epsilon && $(MAKE) clean
#
clean-sigma:
	cd Sigma && $(MAKE) clean
#
clean-bse:
	cd BSE && $(MAKE) clean
#
clean-bse-flavored:
	cd BSE && $(MAKE) clean-flavored
#
clean-plotxct:
	cd PlotXct && $(MAKE) clean
#
clean-meanfield:
	cd MeanField && $(MAKE) clean
#
clean-meanfield-flavored:
	cd MeanField && $(MAKE) clean-flavored
#
clean-manual:
	-$(REMOVE) manual.html
#
clean-sapo:
	cd MeanField/SAPO && $(MAKE) clean
#
#BEGIN_INTERNAL_ONLY
clean-flavored: clean-nonlinearoptics
clean: clean-cumexp
clean-nonlinearoptics:
	cd NonLinearOptics && $(MAKE) clean
clean-cumexp:
	cd GW+C && $(MAKE) clean
clean-abi2bgw:
	cd MeanField/ABINIT && $(MAKE) clean
#END_INTERNAL_ONLY
#
clean-epm:
	cd MeanField/EPM && $(MAKE) clean
#
clean-siesta2bgw:
	cd MeanField/SIESTA && $(MAKE) clean
#
clean-bgw2para:
	cd MeanField/PARATEC && $(MAKE) clean
#
clean-kgrid:
	cd MeanField/ESPRESSO && $(MAKE) clean
#
clean-icm:
	cd MeanField/ICM && $(MAKE) clean
#
clean-utilities:
	cd MeanField/Utilities && $(MAKE) clean
#
clean-surface:
	cd Visual && $(MAKE) clean
#
clean-cpp: clean-surface clean-icm clean-wfn_utils clean-qhull clean-spglib

# Cleanall
#
cleanall: cleanall-bin clean-library clean
#
# Get rid of the symbolic links and the executables they point to
cleanall-bin:
	-cd bin && $(MAKE) cleanall-bin
#
cleanall-cpp: cleanall-surface cleanall-icm clean-wfn_utils cleanall-qhull cleanall-spglib
#
uninstall: cleanall-bin
#
cleanall-common:
	cd $(COMMON) && $(MAKE) cleanall
#
cleanall-epsilon:
	cd Epsilon && $(MAKE) cleanall
#
cleanall-sigma:
	cd Sigma && $(MAKE) cleanall
#
cleanall-bse:
	cd BSE && $(MAKE) cleanall
#
cleanall-plotxct:
	cd PlotXct && $(MAKE) cleanall
#
cleanall-meanfield:
	cd MeanField && $(MAKE) cleanall
#
cleanall-sapo:
	cd MeanField/SAPO && $(MAKE) cleanall
#
#BEGIN_INTERNAL_ONLY
cleanall-nonlinearoptics:
	cd NonLinearOptics && $(MAKE) cleanall
cleanall-cumexp:
	cd GW+C && $(MAKE) cleanall
cleanall-abi2bgw:
	cd MeanField/ABINIT && $(MAKE) cleanall
#END_INTERNAL_ONLY
#
cleanall-epm:
	cd MeanField/EPM && $(MAKE) cleanall
#
cleanall-siesta2bgw:
	cd MeanField/SIESTA && $(MAKE) cleanall
#
cleanall-bgw2para:
	cd MeanField/PARATEC && $(MAKE) cleanall
#
cleanall-kgrid:
	cd MeanField/ESPRESSO && $(MAKE) cleanall
#
cleanall-icm:
	cd MeanField/ICM && $(MAKE) cleanall
#
cleanall-utilities:
	cd MeanField/Utilities && $(MAKE) cleanall
#
cleanall-surface:
	cd Visual && $(MAKE) cleanall

# Clean-keepmod: retains Fortran .mod files. Do not define for C/C++ only targets.
clean-keepmod: clean-keepmod-common clean-keepmod-epsilon clean-keepmod-sigma clean-keepmod-bse clean-keepmod-plotxct \
	clean-keepmod-meanfield
#
clean-keepmod-common:
	cd $(COMMON) && $(MAKE) clean-keepmod
#
clean-keepmod-epsilon:
	cd Epsilon && $(MAKE) clean-keepmod
#
clean-keepmod-sigma:
	cd Sigma && $(MAKE) clean-keepmod
#
clean-keepmod-bse:
	cd BSE && $(MAKE) clean-keepmod
#
clean-keepmod-plotxct:
	cd PlotXct && $(MAKE) clean-keepmod
#
clean-keepmod-meanfield:
	cd MeanField && $(MAKE) clean-keepmod
#
clean-keepmod-sapo:
	cd MeanField/SAPO && $(MAKE) clean-keepmod
#
#BEGIN_INTERNAL_ONLY
cleanall-keepmod: clean-keepmod-nonlinearoptics clean-keepmod-cumexp
clean-keepmod-nonlinearoptics:
	cd NonLinearOptics && $(MAKE) clean-keepmod
clean-keepmod-cumexp:
	cd GW+C && $(MAKE) clean-keepmod
clean-keepmod-abi2bgw:
	cd MeanField/ABINIT && $(MAKE) clean-keepmod
#END_INTERNAL_ONLY
clean-keepmod-epm:
	cd MeanField/EPM && $(MAKE) clean-keepmod
#
clean-keepmod-siesta2bgw:
	cd MeanField/SIESTA && $(MAKE) clean-keepmod
#
clean-keepmod-bgw2para:
	cd MeanField/PARATEC && $(MAKE) clean-keepmod
#
clean-keepmod-kgrid:
	cd MeanField/ESPRESSO && $(MAKE) clean-keepmod
#
clean-keepmod-utilities:
	cd MeanField/Utilities && $(MAKE) clean-keepmod
