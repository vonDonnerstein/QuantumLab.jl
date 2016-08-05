### LIBINT2 ###
# only tested with g++, version as follows though current prerelease tag is 2.1.0-beta2
VERSION=2.2.0-beta1

# for OSX (after a fresh install): (PlATFORM == Darwin)
# 1) install homebrew by following the instructions on brew.sh
# 2) brew install autoconf automake gmp eigen boost
# 3) make

# LIBINT2 BUILD CONFIGURATION
OPT_AM=4 # Optimize maximally for up to angular momentum N (N <= max-am)
MAX_AM=7 # Support Gaussians of angular momentum up to N. Can specify values for each derivative levels as a list N0,N1,N2... 
OPT_LEVEL=3

CXX=clang++
SIMD=avx


# DERIVATIVE VARIABLES (do not modify directly)
PLATFORM = $(shell uname)
SO_EXTENSION=so
ifeq ($(PLATFORM),Darwin)
  SO_EXTENSION=dylib
endif 
CXXFLAGS        = -std=c++11 -O$(OPT_LEVEL) 	
ifeq ($(SIMD),avx)
  CXXFLAGS     += -mavx
endif
ifeq ($(SIMD),sse)
  CXXFLAGS     += -msse -msse2 -msse3 -msse4
endif

# CONFIGURE FLAGS FOR THE COMPILATION OF THE COMPILER
CONFIGURE_GEN_FLAGS =  --with-opt-am=$(OPT_AM) 
CONFIGURE_GEN_FLAGS += --with-max-am=$(MAX_AM) 
CONFIGURE_GEN_FLAGS += --enable-shared 			# build shared libraries
CONFIGURE_GEN_FLAGS += --disable-unrolling  		# Do not unroll shell sets into integrals
ifeq ($(SIMD),avx)
  CONFIGURE_GEN_FLAGS += --with-real-type=libint2::simd::VectorAVXDouble
endif
ifeq ($(SIMD),sse)
  CONFIGURE_GEN_FLAGS += --with-real-type=libint2::simd::VectorSSEDouble
endif
CONFIGURE_GEN_FLAGS += --prefix=$(CURDIR)
CONFIGURE_GEN_FLAGS += --with-incdirs="-I$(CURDIR)/libint2_source/include"
CONFIGURE_GEN_FLAGS += CXXGENFLAGS="$(CXXFLAGS)"
CONFIGURE_GEN_FLAGS += CXXFLAGS="$(CXXFLAGS)"

# CONFIGURE FLAGS FOR THE ACTUAL COMPILATION
CONFIGURE_FLAGS =  --enable-shared 			# build shared libraries
# v--- This is not yet supported by libint (as of 2.2.0-beta1) ---v
#ifeq ($(SIMD),avx)
#  CONFIGURE_FLAGS += --with-real-type=libint2::simd::VectorAVXDouble
#endif
# ^--- This is not yet supported by libint (as of 2.2.0-beta1) ---^
CONFIGURE_FLAGS += --prefix=$(CURDIR)/libint2_binaries
CONFIGURE_FLAGS += CXXFLAGS="$(CXXFLAGS)"




# TARGETS
all:	libint2jl

libint2_configure_the_compilation_of_the_compiler=libint2_build/Makefile
$(libint2_configure_the_compilation_of_the_compiler):
	@echo CONFIGURE_THE_COMPILATION_OF_THE_COMPILER
	rm -rf libint2_source libint2_build
	git clone https://github.com/evaleev/libint.git libint2_source
	cd libint2_source && ./autogen.sh
	mkdir libint2_build
	cd libint2_build && ../libint2_source/configure $(CONFIGURE_GEN_FLAGS)
.PHONY: libint2_configure_the_compilation_of_the_compiler
libint2_configure_the_compilation_of_the_compiler: $(libint2_configure_the_compilation_of_the_compiler)

libint2_compile_the_compiler=libint2_build/libint-$(VERSION).tgz
$(libint2_compile_the_compiler): $(libint2_configure_the_compilation_of_the_compiler)
	@echo COMPILE_THE_COMPILER
	cd libint2_build && $(MAKE)
	cd libint2_build && $(MAKE) export
.PHONY: libint2_compile_the_compiler
libint2_compile_the_compiler: $(libint2_compile_the_compiler)

libint2_extract_and_configure_the_compiler=libint2_build/libint-$(VERSION)/Makefile
$(libint2_extract_and_configure_the_compiler): $(libint2_compile_the_compiler)
	@echo EXTRACT_AND_CONFIGURE_THE_COMPILER
	cd libint2_build && tar -xvzf libint-$(VERSION).tgz
	cd libint2_build/libint-$(VERSION) && ./configure $(CONFIGURE_FLAGS)
.PHONY: libint2_extract_and_configure_the_compiler
libint2_extract_and_configure_the_compiler: $(libint2_extract_and_configure_the_compiler)

libint2_compile=../lib/libint2-$(VERSION).$(SO_EXTENSION)
$(libint2_compile): $(libint2_extract_and_configure_the_compiler)
	@echo COMPILE_LIBINT2
	cd libint2_build/libint-$(VERSION) && $(MAKE)
	cd libint2_build/libint-$(VERSION) && $(MAKE) install
	cp libint2_build/libint-$(VERSION)/lib/.libs/libint2.$(SO_EXTENSION) ../lib/libint2-$(VERSION)-$(SIMD).$(SO_EXTENSION)
	@tput setaf 1
	@tput bold
	@echo -e "WARNING: libint2-$(VERSION).so has been (re-)compiled. \n" \
	                "Please make sure the 'const MAX_AM' in src/SubModules/LibInt2Module.jl \n" \
			"is in line with your actual compile options."
	@tput sgr0
ifeq ($(PLATFORM),Darwin)
	install_name_tool -id @rpath/libint2-$(VERSION)-$(SIMD).dylib ../lib/libint2-$(VERSION)-$(SIMD).dylib
endif
.PHONY: libint2_compile
libint2_compile: $(libint2_compile)







### LIBINT2JL ###
INCLUDE = -Ilibint2_binaries/include/ -Ilibint2_binaries/include/libint2 -I/usr/include/eigen3
ifeq ($(PLATFORM),Darwin)
  INCLUDE = -Ilibint2_binaries/include/ -Ilibint2_binaries/include/libint2 -I/usr/local/opt/eigen/include/eigen3
endif

libint2jl: ../lib/libint2-$(VERSION).$(SO_EXTENSION)
	@echo COMPILE_LIBINT2JL
ifneq ($(PLATFORM),Darwin)
	$(CXX) $(CXXFLAGS) -Wno-deprecated -o ../lib/libint2jl-$(VERSION)-$(SIMD).$(SO_EXTENSION) -shared -fPIC -g -std=c++11 $(INCLUDE) libint2jl.cpp ../lib/libint2-$(VERSION)-$(SIMD).$(SO_EXTENSION)
else
	$(CXX) $(CXXFLAGS) -Wno-deprecated -o ../lib/libint2jl-$(VERSION)-$(SIMD).dylib -dynamiclib -fPIC -g -std=c++11 $(INCLUDE) -install_name libint2jl-$(VERSION)-$(SIMD).dylib libint2jl.cpp ../lib/libint2-$(VERSION)-$(SIMD).dylib 
	install_name_tool -add_rpath @loader_path/. ../lib/libint2jl-$(VERSION)-$(SIMD).dylib
endif



### test ###
test.exe:
	$(CXX) $(CXXFLAGS) -g -o test.exe $(INCLUDE) -Ilibint2_build/libint-2.2.0-beta1/tests/eri libint2_build/libint-2.2.0-beta1/tests/eri/test.cc libint2_binaries/lib/libint2.a

.PHONY: clean
clean:
	rm -rf libint2_build/ libint2_source/ libint2_binaries
	rm ../lib/libint2-$(VERSION)-$(SIMD).so
	rm ../lib/libint2jl-$(VERSION)-$(SIMD).so
