all: libFluidProp.a fluid-prop #libFluidProp.so 

#testC: src/testC.c obj/libFluidProp.o
#	cc src/testC.c -L/home/falke/weis_sp/.local/lib obj/libFluidProp.o -ldl -lCoolProp -o testC -lm -lstdc++
	

fluid-prop: obj/fluid-prop.o obj/libFluidProp.o
	cc obj/fluid-prop.o -lFluidPropC -ldl -lCoolProp -o fluid-prop -lm -lstdc++
	cp fluid-prop ~/bin/

obj/fluid-prop.o: src/fluid-prop.c
	cc -c src/fluid-prop.c
	mkdir -p obj
	mv fluid-prop.o obj/

#libFluidProp.so: src/libFluidProp.cpp obj/libFluidProp.o
#	g++ -c -Wall -Werror -fPIC src/libFluidProp.cpp -o obj/libFluidProp.o -DCOOLPROP_LIB -I../include
#	g++ -shared -lm -o libFluidProp.so obj/libFluidProp.o
#	cp libFluidProp.so ~/.local/lib/

libFluidProp.a: obj/libFluidProp.o
	ar rcs libFluidProp.a obj/libFluidProp.o
	cp libFluidProp.a ~/.local/lib/

obj/libFluidProp.o: src/libFluidProp.cpp
	g++ -c src/libFluidProp.cpp -o obj/libFluidProp.o -DCOOLPROP_LIB -I../include 


install:
	cp libFluidProp.a /usr/local/lib/
	cp src/libFluidProp.h /usr/local/include/
	
