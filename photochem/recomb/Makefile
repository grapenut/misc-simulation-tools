
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
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $<

all: recomb


recomb: $(O_FILES)
	$(CC) $(O_FILES) -o recomb $(LFLAGS)

debug: $(O_FILES)
	$(CC) $(O_FILES) -o debug $(LFLAGS) $(DEBUG)

clean:
	rm $(O_FILES) recomb debug *~ core.*


print.h: data.h
data.h:
logspace.h:
cloudy.h: data.h

cloudy.o: cloudy.h
logspace.o: logspace.h data.h
recomb.o: data.h logspace.h cloudy.h print.h
print.o: print.h
