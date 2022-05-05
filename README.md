# This is a C library to calculate fluid properties. It is a backend to CoolProp but can also use
  the REFPROP database (NIST) 

C library to calculate fluid properties - based heavily on CoolProp
# Requirements: CoolProp from CoolProp.org
  CoolProp needs to be installed and a dynamic linked library needs to be compliled out of the
  source. see: http://www.coolprop.org/dev/coolprop/wrappers/SharedLibrary/index.html
 
# Alternatively Refprop can be used:
  For this, first build libRefprop.so (see https://github.com/jowr/librefprop.so):
 	1. Get a copy:
  		git clone --recursive https://github.com/jowr/librefprop.so.git
	2.  Copy the REFPROP Fortran code to the *fortran* directory.
	3.  Put the *fluids* and *mixtures* folders from REFPROP into the *files* folder.
	4.  Call `make` to prepare the files. 
	5.  Either you use `sudo make install` to copy the files to `/usr/local/lib`, `/usr/local/include` and `/opt/refprop` **or** you run `make install` as normal user to copy the files to `$(HOME)/.refprop/lib`, `$(HOME)/.refprop/include` and `$(HOME)/.refprop`.
		
	
## For linux:
	# Check out the sources for CoolProp
	git clone https://github.com/CoolProp/CoolProp --recursive
	# Move into the folder you just created
	cd CoolProp
	# Make a build folder
	mkdir build && cd build
	# Generate builder (defaults to 64-bit on 64-bit machine)
	cmake .. -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release
	# Build
	cmake --build .
	
	Copy by hand the libraries:
	
	# Change "64" to match your system bitness
	sudo cp libCoolProp.so /usr/local/lib/libCoolProp.so.64.:version:
	pushd /usr/local/lib
	sudo ln -sf libCoolProp.so.64.:version: libCoolProp.so.5
	sudo ln -sf libCoolProp.so.5 libCoolProp.so
	popd
	
	sudo ldconfig

## Compiling and install as user:
	make
	
	To prevent linking errors later on, make sure that:
	LIBRARY_PATH  is set correctly

## install as root:
	sudo make install
