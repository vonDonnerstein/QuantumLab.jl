TOPDIR = /home/henry/devel/libint_build

SRCTOPDIR = /home/henry/devel/libint

#INCLUDE = -I$(TOPDIR)/src/bin -I$(TOPDIR)/src/bin/libint -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(SRCTOPDIR)/include -I$(SRCTOPDIR)/src/bin \
#	  -I$(SRCTOPDIR)/src/lib/libint -I$(SRCTOPDIR)/src/bin/libint
INCLUDE = -I$(TOPDIR)/include -I$(SRCTOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(SRCTOPDIR)/include/libint2

all:	regular shared
regular:	
	g++ -o libinttest.exe -g -std=c++11 -I/usr/include/eigen3 $(INCLUDE) libinttest.cpp $(TOPDIR)/lib/.libs/libint2.a

shared_for_windows:
	x86_64-w64-mingw32-gcc -v -o libinttest.dll -shared -fPIC -g -std=c++11 -I/usr/include/eigen3 $(INCLUDE) -DWIN64 libinttest.cpp -Wl,--out-implib -Wl,--export-all -Wl,--enable-auto-import $(TOPDIR)/lib/.libs/libint2.a -lm -lgmpxx -lgmp

shared:
	g++ -o libinttest.so -shared -fPIC -g -std=c++11 -I/usr/include/eigen3 $(INCLUDE) libinttest.cpp $(TOPDIR)/lib/.libs/libint2.so
	cp libinttest.so /home/henry/devel/QuantumLab/test/
#$(TOPDIR)/src/bin/libint/libINT.a -DHAVE_CONFIG_H  -lm -lgmpxx -lgmp
