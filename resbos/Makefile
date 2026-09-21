# Makefile to compile ResBos and link to CERN Fortran or C++ (ROOT)
# libraries
# Written by P. Nadolsky in April 2008
# To tun it, for example,
# Make resbos_root_WZ_cms_zp
#========================================================================
# User-defined parameters

#TARGET=resbos_root

#FC = g77 #Fortran compiler
FC = gfortran #Fortran compiler
CXX=g++  #C++ compiler

# Compiler parameters with optimization or debugging
# CXXFLAGS0 =  -O3 -g
CXXFLAGS0 =  -g
FFLAGS =  -O3 -g -fno-automatic

#CXXFLAGS0 = -g   
#FFLAGS = -g -fno-automatic


#Uncomment the following two lines to link to Fortran CERN libraries 
# to write PAW HBOOKS and NTUPles
#DMYPAW= -DMYPAW
# PAWLIBS = -Lcernlib -lpacklib
#PAWLIBS = -L/usr/lib/cernlib/2006/lib -lpacklib
#PAWLIBS = -L/usr/lib64/cernlib/2006/lib -lpacklib
#PAWLIBS = -L/usr/lib64/cernlib/2006-g77/lib -lpacklib
#PAWLIBS = -Lcernlib -lpacklib -lkernlib -lnsl
#PAWLIBS = -L/usr/lib64/cernlib/2006-g77/lib -lpacklib -lstdc++

# Uncomment the next four lines to link to CERN C++ libraries to
# enable writing to ROOT ntuples
#
# ROOTLIBS and ROOTINCLUDE are locations of ROOT libraries and header files.
# Find the ROOT directory automatically by running root-config 
# or specify it manually by editing the variable ROOTDIR 
ROOTDIR= $(shell root-config --prefix)
DMYROOT= -DMYROOT
ROOTLIBS     := $(shell root-config  --libs) -lstdc++
ROOTINCLUDE     := -I $(shell root-config --incdir)

# End of user-defined parameters ========================================

LDFLAGS = $(FFLAGS) 
CXXFLAGS=$(CXXFLAGS0) $(DMYROOT) -std=c++17

%.co: %.c
	$(CXX) $(CXXFLAGS) $(DMYROOT) $(ROOTINCLUDE) -c -o $@ $<

# Main executable code
resbos_root: resbos_root.o resbosm.o resbosma.o pentagon.o hvvresph.o froot.co
	$(FC) $(LDFLAGS) -o resbos_root resbos_root.o resbosm.o resbosma.o \
	pentagon.o hvvresph.o froot.co $(PAWLIBS) $(ROOTLIBS) 

# Main executable code
resbos_root_WZ: resbos_root_WZ.o resbosm.o resbosma.o pentagon.o hvvresph.o froot.co
	$(FC) $(LDFLAGS) -o resbos_root_WZ resbos_root_WZ.o resbosm.o resbosma.o \
	pentagon.o hvvresph.o froot.co $(PAWLIBS) $(ROOTLIBS) 

# Main executable code
resbos_root_WZ_cms_zp: resbos_root_WZ_cms_zp.o resbosm.o resbosma.o pentagon.o hvvresph.o froot.co
	$(FC) $(LDFLAGS) -o resbos_root_WZ_cms_zp resbos_root_WZ_cms_zp.o resbosm.o resbosma.o \
	pentagon.o hvvresph.o froot.co $(PAWLIBS) $(ROOTLIBS) 

# A utility to add the hbook files (not needed for the main calculation)
hbook+: hbook+.o 
	$(FC) $(LDFLAGS) -o  hbook+ hbook+.o $(PAWLIBS)

# A utility to add the rootnt files (not needed for the main calculation)
rootnt+: rootnt+.co 
	$(CXX) $(CXXFLAGS) -o  rootnt+ rootnt+.co $(ROOTLIBS)

# A utility to add the gbook files (not needed for the main calculation)
gbook+: gbook+.c
	$(CXX) -c $(CXXFLAGS) -o $@ $<

# Main executable code
tmp: tmp.o froot.co
	$(FC) $(LDFLAGS) tmp.o froot.co $(ROOTLIBS) 


%.o: %.f
	$(FC) -c $(FFLAGS) $<

%.co: %.c
	$(CXX) -c $(CXXFLAGS) $(DMYROOT) $(DMYPAW) $(ROOTINCLUDE) -o $@ $<

clean: 
	rm $ *.o *.co
