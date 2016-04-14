all: binaries.c binaries.h
	mpicc -fPIC -shared binaries.c -o libbinaries.so -lm -lsort -lkdtree -lhdf5
	mv -f libbinaries.so /usr/local/lib/
local: binaries.c binaries.h
	mpicc -fPIC -shared binaries.c -o libbinaries.so -lm -lsort -lkdtree -lhdf5
