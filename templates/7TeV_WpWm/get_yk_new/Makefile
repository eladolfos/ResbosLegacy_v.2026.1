
COMPFLAGS = -O3 -fno-automatic
LINKFLAGS = -O3  -fno-automatic

#COMPFLAGS = -g -fvxt -fno-automatic
#LINKFLAGS = -fvxt -g  -fno-automatic

#All :   get_yk_new   get_yk_new_16line  get_yk_new_h0
All :   get_yk_new

get_yk_new: get_yk_new.o 
	gfortran $(LINKFLAGS) -o  get_yk_new get_yk_new.o 
        
get_yk_new.o: get_yk_new.f 
	gfortran $(COMPFLAGS)  -c get_yk_new.f

get_yk_new_16line: get_yk_new_16line.o 
	gfortran $(LINKFLAGS) -o  get_yk_new_16line get_yk_new_16line.o 
        
get_yk_new_16line.o: get_yk_new_16line.f 
	gfortran $(COMPFLAGS)  -c get_yk_new_16line.f

get_yk_new_h0: get_yk_new_h0.o 
	gfortran $(LINKFLAGS) -o  get_yk_new_h0 get_yk_new_h0.o 
        
get_yk_new_h0.o: get_yk_new_h0.f 
	gfortran $(COMPFLAGS)  -c get_yk_new_h0.f

clean:	
	rm *.o get_yk_new

#end 
