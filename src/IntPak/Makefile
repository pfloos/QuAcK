IDIR =../../include
LDIR =../../lib
BDIR =../../bin
ODIR = obj
SDIR =.
FC = gfortran
FFLAGS = -Wall -O3 -I$(IDIR) 
DFLAGS = -g -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant


LIBS = $(LDIR)/*.a $(LDIR)/slatec/src/static/libslatec.a


SRC = $(wildcard *.f90)

OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(SRC))


$(ODIR)/%.o: %.f90 
	$(FC) -c -o $@ $< $(FFLAGS) 

$(BDIR)/IntPak: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS)

debug: $(OBJ)
	$(FC) -o $(BDIR)/$@ $^ $(FFLAGS) $(LIBS) $(DFLAGS)

clean:
	rm -f $(ODIR)/*.o $(BDIR)/IntPak $(BDIR)/debug
