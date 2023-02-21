all: libFluidProp.a rbc-fluid #libFluidProp.so 

#testC: src/testC.c obj/libFluidProp.o
#	cc src/testC.c -L/home/falke/weis_sp/.local/lib obj/libFluidProp.o -ldl -lCoolProp -o testC -lm -lstdc++
	

rbc-fluid: obj/rbc-fluid.o obj/libFluidProp.o
	cc obj/rbc-fluid.o -lFluidPropC -ldl -lCoolProp -o rbc-fluid -lm -lstdc++
	cp rbc-fluid ~/bin/

obj/rbc-fluid.o: src/rbc-fluid.c
	cc -c src/rbc-fluid.c
	mkdir -p obj
	mv rbc-fluid.o obj/

#libFluidProp.so: src/libFluidProp.cpp obj/libFluidProp.o
#	g++ -c -Wall -Werror -fPIC src/libFluidProp.cpp -o obj/libFluidProp.o -DCOOLPROP_LIB -I../include
#	g++ -shared -lm -o libFluidProp.so obj/libFluidProp.o
#	cp libFluidProp.so ~/.local/lib/

libFluidProp.a: obj/libFluidProp.o 
	ar rcs libFluidProp.a obj/libFluidProp.o
	cp libFluidProp.a ~/.local/lib/
	cp src/libFluidProp.h ~/.local/include/

obj/libFluidProp.o: src/libFluidProp.cpp
	g++ -c src/libFluidProp.cpp -o obj/libFluidProp.o -DCOOLPROP_LIB -I../include 


install:
	cp libFluidProp.a ~/.local/lib/
	cp src/libFluidProp.h ~/.local/include/
	cp libFluidProp.a /usr/local/lib/
	cp src/libFluidProp.h /usr/local/include/
	
