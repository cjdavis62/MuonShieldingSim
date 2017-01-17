# --------------------------------------------------------------
# GNUmakefile for physics list user.  
# JPW. Fri Jul 25 10:39:58 CEST 2003
# --------------------------------------------------------------

.PHONY: all
all: lib bin

name := qshields

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif
ifndef G4PL
  G4PL = LBE
endif

G4VIS_USE := 1 
G4VIS_USE_VRML := 1 
G4VIS_USE_VRMLFILE = 1

include $(G4INSTALL)/config/architecture.gmk

#
# define G4LISTS_BASE, if you have your own physics lists area installed
# point G4LISTS_BASE to the directory, that contains the subdirectory 'lists'.
#
ifndef G4LISTS_BASE
  EXTRALIBS += -L$(G4LIB)/$(G4SYSTEM)
  G4LISTS_BASE = $(G4INSTALL)/hadronic_lists
else
  EXTRALIBS += -L$(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
endif

 EXTRALIBS += -lG4esources
 EXTRALIBS += -l$(G4PL)
#

.PHONY: all
 
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

LDFLAGS  += -L/usr/lib -lc
LOADLIBS += $(G4STATIC) -lG4esources -lc

install: all
	@echo "Installing $(name) to /usr/remote/geant4/bin/$(name)-$(G4VERSION)"
	@install $(G4BIN)/$(G4SYSTEM)/$(name) /usr/remote/geant4/bin/$(name)-$(G4VERSION)
	@scp $(G4BIN)/$(G4SYSTEM)/$(name) crio.mib.infn.it:/usr/remote/geant4/bin/$(name)-$(G4VERSION)
