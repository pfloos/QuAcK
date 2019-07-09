IDIR =../../include
BDIR =../../bin
LDIR =../../lib
ODIR = obj
OODIR = ../IntPak/obj
SDIR =.
FC = gfortran -I$(IDIR) 
AR = libtool
ifeq ($(DEBUG),1)
FFLAGS = -Wall -g -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
else
FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O3
endif

LIBS = ../../lib/*.a
#LIBS = -lblas -llapack

SRCF90 = $(wildcard *.f90)

SRC = $(wildcard *.f)

OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90)) $(patsubst %.f,$(ODIR)/%.o,$(SRC)) 

$(ODIR)/%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) 

$(ODIR)/%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS) 

$(LDIR)/utils.a: $(OBJ)
	$(AR) -static -o $@ $^ 

debug: 
	DEBUG=1 make $(LDIR)/utils.a

clean:
	rm -f $(ODIR)/*.o $(LDIR)/utils.a
