CXX = g++
INCL = -I/data1/r900-1/milos/gsl-1.15/include
CXXFLAGS = -Wno-deprecated -O3 -Wall  $(INCL)
LIBS = -static -L/data1/r900-1/milos/gsl-1.15/lib -lgsl -lgslcblas
LDFLAGS = $(LIBS)
RESULT = radii

SRC = \
  radii.cc

OBJ = $(SRC:.cc=.o)

$(RESULT): $(OBJ) 
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(RESULT) $(LDFLAGS)

clean: 
	rm -f $(OBJ) $(RESULT) core core.*

radii.o: radii.cc function.hh