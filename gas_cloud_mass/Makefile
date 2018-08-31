include ../QuickFlash-1.0.0/settings.make

CCFLAGS += -I../QuickFlash-1.0.0/include -O0 -g
LDFLAGS += ../QuickFlash-1.0.0/lib/libquickflash.a



BIN_FILES = gas_cloud_mass


all : $(BIN_FILES)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o

gas_cloud_mass : gas_cloud_mass.o
	$(LD) -o gas_cloud_mass gas_cloud_mass.o $(LDFLAGS) $(LDFLAGS_HDF5)

gas_cloud_mass.o : gas_cloud_mass.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o gas_cloud_mass.o gas_cloud_mass.cpp



