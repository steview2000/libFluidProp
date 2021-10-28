all: fluid-prop testC

testC: src/testC.c obj/libFluidPropC.o
	cc src/testC.c -L./ obj/libFluidPropC.o -ldl -lCoolProp -o testC -lm -lstdc++
	

fluid-prop: obj/fluid-prop.o libFluidPropC.a
	cc obj/fluid-prop.o -L./ -lFluidPropC -ldl -lCoolProp -o fluid-prop -lm -lstdc++

obj/fluid-prop.o: src/fluid-prop.c
	cc -c src/fluid-prop.c
	mv fluid-prop.o obj/

libFluidPropC.a: src/libFluidPropC.cpp
	g++ -c src/libFluidPropC.cpp -o obj/libFluidPropC.o -DCOOLPROP_LIB -I../include 
	ar rcs libFluidPropC.a obj/libFluidPropC.o

install:
	cp libFluidPropC.a /usr/local/lib/
	cp src/libFluidPropC.h /usr/local/include/
	
