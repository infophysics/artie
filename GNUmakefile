
name := Artie
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: lib bin

# ROOT support
CPPFLAGS += -I$(shell root-config --incdir) -g
EXTRALIBS = $(shell root-config --glibs)

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
