
CC = g++
SRC = ../../c13.05/source

CFLAGS = -Wall -I$(SRC) -I$(SRC)/sys_gcc
LFLAGS = -L$(SRC) -L$(SRC)/sys_gcc -lcloudy
DEBUG = -g

C_FILES := $(shell ls *.cpp)
O_FILES = $(C_FILES:.cpp=.o)

.SUFFIX:
.SUFFIX: .cpp .h

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

all: photoeq


photoeq: photoeq.o
	$(CC) $(O_FILES) -o photoeq $(LFLAGS)

debug: 
	$(CC) $(O_FILES) -o debug $(LFLAGS) $(DEBUG)

clean:
	rm $(O_FILES) photoeq debug *~


photoeq.o: logspace.o cloudy.o data.h logspace.h cloudy.h
logspace.o:
cloudy.o:

data.h:
logspace.h:
cloudy.h:

