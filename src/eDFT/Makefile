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

LIBS = ~/Dropbox/quack/lib/*.a
#LIBS = -lblas -llapack

SRCF90 = $(wildcard *.f90)

SRC = $(wildcard *.f)

OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90)) $(patsubst %.f,$(ODIR)/%.o,$(SRC))

$(ODIR)/%.o: %.f90 
	$(FC) -c -o $@ $< $(FFLAGS) 

$(ODIR)/%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS) 

$(BDIR)/eDFT: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

debug: 
	DEBUG=1 make $(BDIR)/eDFT

profil: 
	PROFIL=1 make $(BDIR)/eDFT

clean:
	rm -f $(ODIR)/*.o $(BDIR)/eDFT $(BDIR)/debug
