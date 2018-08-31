include ../QuickFlash-1.0.0/settings.make

CCFLAGS += -I../QuickFlash-1.0.0/include -O0 -g
LDFLAGS += ../QuickFlash-1.0.0/lib/libquickflash.a



BIN_FILES = extract_dens_temp


all : $(BIN_FILES)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o

extract_dens_temp : extract_dens_temp.o
	$(LD) -o extract_dens_temp extract_dens_temp.o $(LDFLAGS) $(LDFLAGS_HDF5)

extract_dens_temp.o : extract_dens_temp.cpp kdtree.h pointbox.h
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o extract_dens_temp.o extract_dens_temp.cpp



