#Makefile for Legacy
#Author: Pavel Nadolsky, July 8, 1999
#Use tabs whenever you see them

#COMP=g77296
COMP=gfortran
CXX=g++

# For cteq4 --  ct4m.ini
#PRZPAC=PrzPac701
#EVLPAC=EvlPac721
#QCDPAC=QcdPac66
#UTLPAC=UtlPac6
#
# For cteq5 --  ct5m.ini
#PRZPAC=PrzPac732
#EVLPAC=EvlPac732
#QCDPAC=QcdPac70
#UTLPAC=UtlPac601
#
# For cteq6
#PRZPAC=PrzPac801
#EVLPAC=EvlPac801
#QCDPAC=QcdPac70
#UTLPAC=UtlPac601
#
# For cteq6.1
#PRZPAC=PrzPac92a
#PRZPAC=PrzPac01e
#EVLPAC=EvlPac92an
#QCDPAC=QcdPac92a
#UTLPAC=UtlPac92a
#
# For cteq6.1 and cteq6.6 (via link to 6.5)
# PRZPAC=PrzPac01e_66
# EVLPAC=EvlPac92a
# QCDPAC=QcdPac92a
# UTLPAC=UtlPac92a

# Jan 31, 2009 
PRZPAC=PrzPac02b
EVLPAC=EvlPac02b
QCDPAC=QcdPac02b
UTLPAC=UtlPac02b
EWKPAC=EwkPac02b
#LHAPAC=LhaPac02b
#
TARGET=w_asym
#
#COMPFLAGS=-c -g 
#COMPFLAGS=-c -O 
#
COMPFLAGS=-c -O2 -g -w -fno-automatic 
#COMPFLAGS=-c -O -fno-automatic
#COMPFLAGS=-c -O -fvxt -fno-automatic
#COMPFLAGS=-c -g -w -fno-automatic
#
LHAPDFINCLUDE = $(shell lhapdf-config --incdir)
LHAPDFLIBS = $(shell lhapdf-config --ldflags) -lLHAPDF
LHAPDFCPP = $(shell lhapdf-config --cppflags) --std=c++11
# 
LINKFLAGS= -O2 -fno-automatic -o $(TARGET)
$(TARGET): w_asym.o asym.o gauss.o util.o CT14Pdf.o $(EVLPAC).o $(PRZPAC).o $(QCDPAC).o\
                          $(UTLPAC).o $(EWKPAC).o lhapdf.co 
	$(COMP) $(LINKFLAGS) w_asym.o asym.o gauss.o util.o CT14Pdf.o $(EVLPAC).o \
	     $(PRZPAC).o $(QCDPAC).o $(UTLPAC).o $(EWKPAC).o  lhapdf.co $(LHAPDFLIBS) -lstdc++

%.co: %.cpp
	$(CXX) $(CXXFLAGS) $(LHAPDFCPP) -c -o $@ $<

w_asym.o: w_asym.f 
	$(COMP) $(COMPFLAGS)  w_asym.f

asym.o: asym.f 
	$(COMP) $(COMPFLAGS)  asym.f 

gauss.o: gauss.f 
	$(COMP) $(COMPFLAGS)  gauss.f

util.o: util.f 
	$(COMP) $(COMPFLAGS)  util.f

CT14Pdf.o: CT14Pdf.f 
	$(COMP) $(COMPFLAGS)  CT14Pdf.f

$(EVLPAC).o: $(EVLPAC).for
	$(COMP) $(COMPFLAGS)  $(EVLPAC).for

$(PRZPAC).o: $(PRZPAC).for
	$(COMP) $(COMPFLAGS)  $(PRZPAC).for
	
$(QCDPAC).o: $(QCDPAC).for
	$(COMP) $(COMPFLAGS)  $(QCDPAC).for
	
$(UTLPAC).o: $(UTLPAC).for
	$(COMP) $(COMPFLAGS)  $(UTLPAC).for
	
$(EWKPAC).o: $(EWKPAC).for
	$(COMP) $(COMPFLAGS)  $(EWKPAC).for
	
#$(LHAPAC).o: $(LHAPAC).for
#	$(COMP) $(COMPFLAGS)  $(LHAPAC).for
	
clean:;		@rm -f *.o *.co core w_asym

#end 
