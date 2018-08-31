include ../QuickFlash-1.0.0/settings.make

CCFLAGS += -I../QuickFlash-1.0.0/include
LDFLAGS += ../QuickFlash-1.0.0/lib/libquickflash.a



BIN_FILES = halo


all : $(BIN_FILES)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o

halo : halo.o
	$(LD) -o halo halo.o $(LDFLAGS) $(LDFLAGS_HDF5)

halo.o : halo.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o halo.o halo.cpp
