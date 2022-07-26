ll: libFluidPropC.so libFluidPropC.a fluid-prop

#testC: src/testC.c obj/libFluidPropC.o
#	cc src/testC.c -L/home/falke/weis_sp/.local/lib obj/libFluidPropC.o -ldl -lCoolProp -o testC -lm -lstdc++
	

fluid-prop: obj/fluid-prop.o obj/libFluidPropC.o
	cc obj/fluid-prop.o -lFluidPropC -ldl -lCoolProp -o fluid-prop -lm -lstdc++
	cp fluid-prop ~/bin/

obj/fluid-prop.o: src/fluid-prop.c
	cc -c src/fluid-prop.c
	mkdir -p obj
	mv fluid-prop.o obj/

libFluidPropC.so: src/libFluidPropC.cpp obj/libFluidPropC.o
	g++ -c -Wall -Werror -fPIC src/libFluidPropC.cpp -o obj/libFluidPropC.o -DCOOLPROP_LIB -I../include
	g++ -shared -lm -o libFluidPropC.so obj/libFluidPropC.o
	cp libFluidPropC.so ~/.local/lib/

libFluidPropC.a: obj/libFluidPropC.o
	ar rcs libFluidPropC.a obj/libFluidPropC.o
	cp libFluidPropC.a ~/.local/lib/

obj/libFluidPropC.o: src/libFluidPropC.cpp
	g++ -c src/libFluidPropC.cpp -o obj/libFluidPropC.o -DCOOLPROP_LIB -I../include 


install:
	cp libFluidPropC.a /usr/local/lib/
	cp src/libFluidPropC.h /usr/local/include/
	
