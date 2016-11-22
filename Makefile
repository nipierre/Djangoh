CC=gcc
CX=g++
CXX=gfortran

CCFLAGS = -g -O1 -W -Wall #-pedantic -fPIC
ROOTFLAGS = `root-config --cflags --glibs`
ROOTVERSION = -D ROOT5
PFMFLAG = -DCCAGE

PROGRAM=djangoh
DJANGOH=/sps/compass/npierre/djangoh
LHAPDF=/sps/compass/npierre/lhapdf5
LHAPDF_LIB=$(LHAPDF)/lib
LIBS=-L$(LHAPDF_LIB) -lLHAPDF -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -L/afs/in2p3.fr/cernlib/@sys/pro/lib -Wl,-static -lmathlib -lpacklib -lkernlib -Wl,-dy -llapack -lm -lnsl -lcrypt -ldl -lg2c
#LIBS=-L$(LHAPDF_LIB) -lLHAPDF $(shell cernlib mathlib kernlib packlib)

CFFLAGS = -Wall -pedantic -fno-automatic

OBJS = $(addsuffix .o, $(basename $(SRCS)))

SRCS = djangoh_h.f djangoh_l.f djangoh_u.f djangoh_t.f \
	sophia.f pythia-6.4.24.f jetset7409.f polpdf.f

%.o: %.f
	$(CXX) $(CFFLAGS) -o $@ -c $<

lhaglue.o: lhaglue.f
	$(CXX) $(CFFLAGS) -ffree-form -c lhaglue.f -o lhaglue.o

gmc_random.o: gmc_random.f
	$(CXX) $(CFFLAGS) -ffree-form -c gmc_random.f -o gmc_random.o
	
djangoh: $(OBJS) lhaglue.o gmc_random.o
	$(CXX) $(CFFLAGS) -o $@ $(OBJS) lhaglue.o gmc_random.o $(LIBS)
	if [ -a compass-nc-test_out.dat ]; then rm compass-nc-test_*; fi;
	if [ -a luevents.dat ]; then rm luevents.dat; fi;
	./djangoh < compass-nc-test.in

pdist: pdist.cc
	$(CX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	./pdist

mult: multiplicity.cc multiplicity.h
	$(CX) $(CCFLAGS) $(PFMFLAG) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< 
	./mult

mult_d: multiplicity.cc multiplicity.h
	$(CX) $(CCFLAGS) $(PFMFLAG) $(ROOTFLAGS) $(ROOTVERSION) -DDEBUG -o $@ $<
	./mult_d

clean:
	rm -f *.o $(PROGRAM)
	rm -f lhaglue.o
	if [ -a mult ]; then rm mult; fi;
	if [ -a mult_d ]; then rm mult_d; fi;
	if [ -a pdist ]; then rm pdist; fi;

cleanfile:
	if [ -a compass-nc-test_out.dat ]; then rm compass-nc-test_*; fi;
	if [ -a luevents.dat ]; then rm luevents.dat; fi;
