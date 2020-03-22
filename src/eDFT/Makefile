IDIR =../../include
BDIR =../../bin
ODIR = obj
OODIR = ../IntPak/obj
SDIR =.
FC = gfortran -I$(IDIR) -Wall -g -Wno-unused -Wno-unused-dummy-argument
ifeq ($(DEBUG),1)
FFLAGS = -Wall -g -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
else
FFLAGS = -O3
endif

ifeq ($(PROFIL),1)
  FC += -p -fno-inline
endif

LIBS = -lc++ ../../lib/*.a
LIBS += -lblas -llapack

SRCF90 = $(wildcard *.f90)

SRC = $(wildcard *.F)

OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90)) $(patsubst %.F,$(ODIR)/%.o,$(SRC))

$(BDIR)/eDFT: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

numgrid.mod: numgrid.f90 
	$(FC) -c -o $@ $< $(FFLAGS) 

$(ODIR)/%.o: %.f90 numgrid.mod
	$(FC) -c -o $@ $< $(FFLAGS) 

$(ODIR)/%.o: %.F numgrid.mod
	$(FC) -c -o $@ $< $(FFLAGS) 

debug: 
	DEBUG=1 make $(BDIR)/eDFT

profil: 
	PROFIL=1 make $(BDIR)/eDFT

clean:
	rm -f $(ODIR)/*.o $(BDIR)/eDFT $(BDIR)/debug numgrid.mod
