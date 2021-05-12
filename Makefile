all: fluid-prop


fluid-prop: obj/fluid-prop.o libFluidPropC.a
	cc obj/fluid-prop.o -lFluidPropC -ldl -lCoolProp -o fluid-prop

obj/fluid-prop.o: src/fluid-prop.c
	cc -c src/fluid-prop.c
	mv fluid-prop.o obj/

libFluidPropC.a: src/libFluidPropC.cpp
	g++ -c src/libFluidPropC.cpp -o obj/libFluidPropC.o -DCOOLPROP_LIB -I../include 
	ar rcs libFluidPropC.a obj/libFluidPropC.o

install:
	cp libFluidPropC.a /usr/local/lib/
	
