CXX = g++
INCL = -I/data1/r900-1/milos/gsl-1.15/include -I.
CXXFLAGS = -g -O2 -Wall $(INCL)
LIBS = -L/data1/r900-1/milos/gsl-1.15/lib -lgsl -lgslcblas
LDFLAGS = $(LIBS)
RESULT = main

SRC = \
  main.cc

OBJ = $(SRC:.cc=.o)

$(RESULT): $(OBJ) 
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(RESULT) $(LDFLAGS)

clean: 
	rm -f $(OBJ) $(RESULT)
