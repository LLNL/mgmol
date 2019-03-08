# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory.
# Written by J.-L. Fattebert, D. Osei-Kuffuor and I.S. Dunn.
# LLNL-CODE-743438
# All rights reserved.
# This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
# Please also read this link https://github.com/llnl/mgmol/LICENSE
# 
CURRENT_DIR = $(PWD)
MAIN_DIR = $(CURRENT_DIR)/..
INSTALL_DIR = $(MAIN_DIR)/install

.PHONY: DEFAULT COMMON_TARGETS makedirs clean distclean depend

include $(INSTALL_DIR)/Makefile.in

GITHASH := $(shell git show -s --format="%h" )
 
BINDIR = ../bin

MODULES = tools pb DistMatrix sparse_linear_algebra \
          linear_algebra numerical_kernels local_matrices \
          radial
 
SRC =   \
 GrassmanCGFactory.cc \
 DMStrategyFactory.cc \
 manage_memory.cc \
 SubCell.cc \
 SinCosOps.cc \
 PolakRibiereSolver.cc \
 EnergySpreadPenalty.cc \
 SpreadPenalty.cc \
 SpreadPenaltyVolume.cc \
 setup.cc \
 MVP_DMStrategy.cc \
 HamiltonianMVP_DMStrategy.cc \
 MVPSolver.cc \
 HamiltonianMVPSolver.cc \
 OrbitalsPreconditioning.cc \
 DFTsolver.cc \
 NonOrthoDMStrategy.cc \
 FullyOccupiedNonOrthoDMStrategy.cc \
 EigenDMStrategy.cc \
 Masks4Orbitals.cc \
 AOMMprojector.cc \
 hdf_tools.cc \
 DistributedIonicData.cc \
 IonData.cc \
 SubspaceProjector.cc \
 MDfiles.cc \
 OrbitalsExtrapolation.cc \
 OrbitalsExtrapolationOrder2.cc \
 OrbitalsExtrapolationOrder3.cc \
 runfire.cc \
 FIRE.cc \
 IonicAlgorithm.cc \
 FIRE_IonicStepper.cc \
 tools.cc \
 MGmol.cc \
 MGmol_NEB.cc \
 ABPG.cc \
 GrassmanLineMinimization.cc \
 GrassmanCG.cc \
 GrassmanCGSparse.cc \
 LBFGS.cc \
 IonicStepper.cc \
 Energy.cc \
 GramMatrix.cc \
 DensityMatrix.cc \
 DensityMatrixSparse.cc \
 ConstraintSet.cc \
 DistanceConstraint.cc \
 MultiDistanceConstraint.cc \
 getLAeigen.cc \
 LDAFunctional.cc \
 LDAonGrid.cc \
 PBEFunctional.cc \
 PBEonGrid.cc \
 PBEonGridSpin.cc \
 Electrostatic.cc \
 ProjectedMatrices.cc \
 ProjectedMatricesSparse.cc \
 ProjectedMatrices2N.cc \
 ProjectedMatricesMehrstellen.cc \
 ReplicatedWorkSpace.cc \
 Rho.cc \
 Mesh.cc \
 Control.cc \
 Hamiltonian.cc \
 BlockVector.cc \
 KBPsiMatrixInterface.cc \
 KBPsiMatrixSparse.cc \
 KBprojectorSparse.cc \
 MultipoleExpansion.cc \
 SpreadsAndCenters.cc \
 Preconditioning.cc \
 OrbitalsTransform.cc \
 NOLMOTransform.cc \
 LocalizationRegions.cc \
 Hartree.cc \
 ShiftedHartree.cc \
 PBdiel.cc \
 Species.cc \
 HDFrestart.cc \
 MasksSet.cc \
 FunctionsPacking.cc \
 ColoredRegions.cc \
 LocGridOrbitals.cc \
 ExtendedGridOrbitals.cc \
 mcstep.cc \
 MD_IonicStepper.cc \
 LBFGS_IonicStepper.cc \
 lbfgsrlx.cc \
 AndersonMix.cc \
 MLWFTransform.cc \
 Ion.cc \
 GridMask.cc \
 GridMaskMult.cc \
 GridMaskMax.cc \
 Ions.cc \
 restart.cc \
 md.cc \
 get_vnlpsi.cc \
 quench.cc \
 mlwf.cc \
 readInput.cc \
 Forces.cc \
 computeHij.cc \
 Potentials.cc \
 BasicDataDistributors.cc \
 ShortSightedInverse.cc \
 PCGSolver_Diel.cc \
 PBdiel_CG.cc \
 PCGSolver.cc \
 Hartree_CG.cc \
 ClusterOrbitals.cc
           
MGMOL_SRC = \
       main.cc

ALL_SRC = $(SRC)
ALL_SRC += $(MGMOL_SRC)

include $(patsubst %,%/module.mk, $(MODULES))

BUILD_SUFFIX := $(PLAT)
OBJECTS = $(SRC:%.cc=$(OBJDIR)/%.o)
OBJECTS += $(MAIN_SRC:%.cc=$(OBJDIR)/%.o)

OBJDIR := objs-$(PLAT)
EXENAME = mgmol

DEFS += -std=c++11 $(MACHINE_FLAGS) -DGITHASH="$(GITHASH)" -DUSE_MPI -DUSE_DIS_MAT -DSCALAPACK

INCLUDE += $(patsubst %,-I%,$(MODULES))

LD = $(CXX)

TARGET_LIB =



DEFAULT: opt

COMMON_TARGETS: makedirs

opt:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/mgmol-$(PLAT) \
	BUILD_SUFFIX=$(PLAT) OBJDIR=objs-$(PLAT) \
	CXX="$(MPCC)" \
	CLIBS="$(NAT_COMPLIB)" \
	CXXFLAGS="$(DEFS) $(INCLUDE) $(OPT_FLAGS)" \
	MAIN_SRC=$(MGMOL_SRC) EXENAME=mgmol-$(PLAT)

gnu:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/mgmol-$(PLAT)-gnu \
	BUILD_SUFFIX=$(PLAT) OBJDIR=objs-$(PLAT) \
	CXX="$(MPGCC)" \
	CLIBS="$(GCC_COMPLIB)" \
	CXXFLAGS="$(DEFS) $(INCLUDE) $(GCC_OPT_FLAGS)" \
	MAIN_SRC=$(MGMOL_SRC) EXENAME=mgmol-$(PLAT)-gnu

mp:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/mgmol-$(PLAT)-mp \
	BUILD_SUFFIX=$(PLAT) OBJDIR=objs-$(PLAT) \
	CXX="$(MPCC)" \
	CLIBS="$(NAT_COMPLIB)" \
	TARGET_LIB="$(MPBLASLIB)" \
	CXXFLAGS="$(DEFS) -DUSE_MP $(INCLUDE) $(OPT_FLAGS)" \
	MAIN_SRC=$(MGMOL_SRC) EXENAME=mgmol-$(PLAT)-mp

debug:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/mgmol-$(PLAT)-dbg \
	BUILD_SUFFIX=$(PLAT) OBJDIR=objs-$(PLAT) \
	CXX="$(MPCC)" \
	CLIBS="$(NAT_COMPLIB)" \
	CXXFLAGS="$(DEFS) $(INCLUDE) $(DEBUG_FLAGS)" \
	MAIN_SRC=$(MGMOL_SRC) EXENAME=mgmol-$(PLAT)-dbg

debug-mp:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/mgmol-$(PLAT)-dbg-mp \
	BUILD_SUFFIX=$(PLAT) OBJDIR=objs-$(PLAT) \
	CXX="$(MPCC)" \
	CLIBS="$(NAT_COMPLIB)" \
	TARGET_LIB="$(MPBLASLIB)" \
	CXXFLAGS="$(DEFS) -DUSE_MP $(INCLUDE) $(DEBUG_FLAGS)" \
	MAIN_SRC=$(MGMOL_SRC) EXENAME=mgmol-$(PLAT)-dbg-mp

gprof:
	@$(MAKE) --no-print-directory COMMON_TARGETS \
	$(BINDIR)/mgmol-$(PLAT) \
	BUILD_SUFFIX=$(PLAT) OBJDIR=objs-$(PLAT) \
	CXX="$(MPCC) -pg" \
	CLIBS="$(NAT_COMPLIB)" \
	CXXFLAGS="$(DEFS) $(INCLUDE) $(OPT_FLAGS) -pg" \
	MAIN_SRC=$(MGMOL_SRC) EXENAME=mgmol-$(PLAT)

.SUFFIXES: .cc

$(OBJDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

LDFLAGS = $(MGMOL_LD_LIBRARY_PATH) $(TARGET_LIB) $(PLIBS) $(LIBS) $(MPI_LIB) $(COMMONLIBS) $(CLIBS)

$(BINDIR)/$(EXENAME): $(OBJECTS)
	$(LD) -o $(BINDIR)/$(EXENAME) $(OBJECTS) $(LDFLAGS)

makedirs:
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR) ; fi
	$(foreach DIR, $(MODULES), @mkdir -p $(OBJDIR)/$(DIR))
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR) ; fi

#------------------------------------------------------------------------------
# generate dependencies in makefile: use -Y to avoid library header files
# that are likely to be different on other platforms.
depend:
	@makedepend -f .depend_$(BUILD_SUFFIX) -p'$(OBJDIR)/' -o '.o' -Y $(ALL_SRC) -I$(patsubst %, %, $(MODULES)) 2> /dev/null

.depend_$(BUILD_SUFFIX): $(ALL_SRC)
	@touch .depend_$(BUILD_SUFFIX)
	@$(MAKE) depend BUILD_SUFFIX=$(BUILD_SUFFIX)
#------------------------------------------------------------------------------
#  Cleanup object files
clean :
	rm -rf $(OBJDIR)/*.o
	rm -rf $(OBJDIR)/*/*.o
	rm -f $(BINDIR)/$(EXENAME)-$(PLAT)*
	rm -f .depend_$(PLAT)*

distclean: clean
	rm -rf objs-*
	rm -f .depend*
	rm -f $(BINDIR)/*

ifdef BUILD_SUFFIX
  -include .depend_$(BUILD_SUFFIX)
endif
