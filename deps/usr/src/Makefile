### LIBINT2 ###
# only tested with g++, version as follows though current prerelease tag is 2.1.0-beta2
VERSION='2.2.0-beta1'

# LIBINT2 BUILD CONFIGURATION
OPT_AM=4 # Optimize maximally for up to angular momentum N (N <= max-am)
MAX_AM=7 # Support Gaussians of angular momentum up to N. Can specify values for each derivative levels as a list N0,N1,N2... 
OPT_LEVEL=3

OPT_AM=1
MAX_AM=2
OPT_LEVEL=0


# DERIVATIVE VARIABLES (do not modify directly)
CXXFLAGS        =  -std=c++11 -mavx -O$(OPT_LEVEL) 	
# CONFIGURE FLAGS FOR THE COMPILATION OF THE COMPILER
CONFIGURE_GEN_FLAGS =  --with-opt-am=$(OPT_AM) 
CONFIGURE_GEN_FLAGS += --with-max-am=$(MAX_AM) 
CONFIGURE_GEN_FLAGS += --enable-shared 			# build shared libraries
CONFIGURE_GEN_FLAGS += --disable-unrolling  		# Do not unroll shell sets into integrals
CONFIGURE_GEN_FLAGS += --with-real-type=libint2::simd::VectorAVXDouble
CONFIGURE_GEN_FLAGS += --prefix=$(CURDIR)
CONFIGURE_GEN_FLAGS += --with-incdirs="-I$(CURDIR)/libint2_source/include"
CONFIGURE_GEN_FLAGS += CXXGENFLAGS="$(CXXFLAGS)"
CONFIGURE_GEN_FLAGS += CXXFLAGS="$(CXXFLAGS)"
# CONFIGURE FLAGS FOR THE ACTUAL COMPILATION
CONFIGURE_FLAGS =  --enable-shared 			# build shared libraries
CONFIGURE_FLAGS += --with-real-type=libint2::simd::VectorAVXDouble
CONFIGURE_FLAGS += --prefix=$(CURDIR)/libint2_binaries
CONFIGURE_FLAGS += CXXFLAGS="$(CXXFLAGS)"




# TARGETS
libint2_configure_the_compilation_of_the_compiler=libint2_build/Makefile
$(libint2_configure_the_compilation_of_the_compiler):
	@echo CONFIGURE_THE_COMPILATION_OF_THE_COMPILER
	rm -rf libint2_source libint2_build
	git clone https://github.com/evaleev/libint.git libint2_source
	cd libint2_source && ./autogen.sh
	mkdir libint2_build
	cd libint2_build && ../libint2_source/configure $(CONFIGURE_GEN_FLAGS)

libint2_compile_the_compiler=libint2_build/libint-$(VERSION).tgz
$(libint2_compile_the_compiler): $(libint2_configure_the_compilation_of_the_compiler)
	@echo COMPILE_THE_COMPILER
	cd libint2_build && $(MAKE)
	cd libint2_build && $(MAKE) export

libint2_extract_and_configure_the_compiler=libint2_build/libint-2.2.0-beta1/Makefile
$(libint2_extract_and_configure_the_compiler): $(libint2_compile_the_compiler)
	@echo EXTRACT_AND_CONFIGURE_THE_COMPILER
	cd libint2_build && tar -xvzf libint-$(VERSION).tgz
	cd libint2_build/libint-$(VERSION) && ./configure $(CONFIGURE_FLAGS)

libint2_compile=../lib/libint2-$(VERSION).so
$(libint2_compile): $(libint2_extract_and_configure_the_compiler)
	@echo COMPILE_LIBINT2
	cd libint2_build/libint-$(VERSION) && $(MAKE)
	cd libint2_build/libint-$(VERSION) && $(MAKE) install
	cp libint2_build/libint-$(VERSION)/lib/.libs/libint2.so ../lib/libint2-$(VERSION).so

.PHONY: libint2_compile
libint2_compile: $(libint2_compile)

.PHONY: libint2_extract_and_configure_the_compiler
libint2_extract_and_configure_the_compiler: $(libint2_extract_and_configure_the_compiler)





### LIBINT2JL ###
INCLUDE = -Ilibint2_build/include -Ilibint2_source/include -Ilibint2_build/include/libint2 -Ilibint2_source/include/libint2 -Ilibint2_binaries/include/

libint2jl: ../lib/libint2-$(VERSION).so
	@echo COMPILE_LIBINT2JL
	g++ -o ../lib/libint2jl.so -shared -fPIC -g -std=c++11 -I/usr/include/eigen3 $(INCLUDE) libint2jl.cpp ../lib/libint2-$(VERSION).so 
