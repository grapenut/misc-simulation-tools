include ../QuickFlash-1.0.0/settings.make

CCFLAGS += -I../QuickFlash-1.0.0/include
LDFLAGS += ../QuickFlash-1.0.0/lib/libquickflash.a



BIN_FILES = phase hist entropy


all : $(BIN_FILES)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o

phase : phase.o
	$(LD) -o phase phase.o $(LDFLAGS) $(LDFLAGS_HDF5)

phase.o : phase.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o phase.o phase.cpp

hist : hist.o
	$(LD) -o hist hist.o $(LDFLAGS) $(LDFLAGS_HDF5)

hist.o : hist.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o hist.o hist.cpp

entropy : entropy.o
	$(LD) -o entropy entropy.o $(LDFLAGS) $(LDFLAGS_HDF5)

entropy.o : entropy.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o entropy.o entropy.cpp

