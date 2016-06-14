TOPDIR    ?= /home/henry/devel/libint_build
SRCTOPDIR ?= /home/henry/devel/libint

INCLUDE = -I$(TOPDIR)/include -I$(SRCTOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(SRCTOPDIR)/include/libint2

sharedobject:
	g++ -o ../lib/libint2jl.so -shared -fPIC -g -std=c++11 -I/usr/include/eigen3 $(INCLUDE) libint2jl.cpp $(TOPDIR)/lib/.libs/libint2.so
