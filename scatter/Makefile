include ../QuickFlash-1.0.0/settings.make

CCFLAGS += -I../QuickFlash-1.0.0/include -O0 -g
LDFLAGS += ../QuickFlash-1.0.0/lib/libquickflash.a



BIN_FILES = scatter


all : $(BIN_FILES)

clean : clean-obj
	rm -f $(BIN_FILES)

clean-obj :
	rm -f *.o

scatter : scatter.o
	$(LD) -o scatter scatter.o $(LDFLAGS) $(LDFLAGS_HDF5)

scatter.o : scatter.cpp
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o scatter.o scatter.cpp



