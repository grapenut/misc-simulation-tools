include ../QuickFlash-1.0.0/settings.make

CCFLAGS += -I../QuickFlash-1.0.0/include -O0 -g
LDFLAGS += ../QuickFlash-1.0.0/lib/libquickflash.a



BIN_FILES = radial_profile


all : $(BIN_FILES)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o

radial_profile : radial_profile.o
	$(LD) -o radial_profile radial_profile.o $(LDFLAGS) $(LDFLAGS_HDF5)

radial_profile.o : radial_profile.cpp
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o radial_profile.o radial_profile.cpp



