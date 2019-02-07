IDIR =../../include
BDIR =../../bin
LDIR =../../lib
ODIR = obj
SDIR =.
FC = gfortran -g -I$(IDIR) 
ifeq ($(DEBUG),1)
FFLAGS = -Wall -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
else
FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O2 
endif

LIBS = $(LDIR)/*.a $(LDIR)/slatec/src/static/libslatec.a 

SRC = $(wildcard *.f90)

OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(SRC))


$(ODIR)/%.o: %.f90 
	$(FC) -c -o $@ $< $(FFLAGS) 2>&1

$(BDIR)/MCQC: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)
#	$(FC) -o $(BDIR)/$@ $^ $(FFLAGS) $(LIBS)

debug: 
	DEBUG=1 make clean $(BDIR)/MCQC

clean:
	rm -f $(ODIR)/*.o $(BDIR)/MCQC $(BDIR)/debug
