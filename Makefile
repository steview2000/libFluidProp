all: testC

testC: testC.o libFluidPropC.a
	cc testC.o -lFluidPropC -ldl -lCoolProp -o testC

testC.o: testC.c
	cc -c testC.c

libFluidPropC.a: libFluidPropC.cpp
	g++ -c libFluidPropC.cpp -o libFluidPropC.o -DCOOLPROP_LIB -I../include 
	ar rcs libFluidPropC.a libFluidPropC.o

install:
	cp libFluidPropC.a /usr/local/lib/
	
