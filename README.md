# This is a C library to calculate fluid properties. It is heavily based on CoolProp

C library to calculate fluid properties - based heavily on CoolProp
# Requirements: CoolProp from CoolProp.org
  CoolProp needs to be installed and a dynamic linked library needs to be compliled out of the
  source. see: http://www.coolprop.org/dev/coolprop/wrappers/SharedLibrary/index.html
  
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
	
	# Change "32" to match your system bitness
	sudo cp libCoolProp.so /usr/local/lib/libCoolProp.so.32.:version:
	pushd /usr/local/lib
	sudo ln -sf libCoolProp.so.32.:version: libCoolProp.so.5
	sudo ln -sf libCoolProp.so.5 libCoolProp.so
	popd
	
