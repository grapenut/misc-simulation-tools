include ../QuickFlash-1.0.0/settings.make

CCFLAGS += -I../QuickFlash-1.0.0/include
LDFLAGS += ../QuickFlash-1.0.0/lib/libquickflash.a



BIN_FILES = parts phase hist


all : $(BIN_FILES)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o

parts : parts.o
	$(LD) -o parts parts.o $(LDFLAGS) $(LDFLAGS_HDF5)

parts.o : parts.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o parts.o parts.cpp

phase : phase.o
	$(LD) -o phase phase.o $(LDFLAGS) $(LDFLAGS_HDF5)

phase.o : phase.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o phase.o phase.cpp

hist : hist.o
	$(LD) -o hist hist.o $(LDFLAGS) $(LDFLAGS_HDF5)

hist.o : hist.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o hist.o hist.cpp

