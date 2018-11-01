# We require the C11 standard for testing.
CCFLAGS=-std=c++11 -O2
LDFLAGS=-lm
# list all of object files
OBJ=myfun.o mmio.o
EXE=pr3
# compiler
CC=g++

all: $(EXE)

$(EXE): $(OBJ) main.o
	$(CC) $(CCFLAGS) -o $@ $^ $(LDFLAGS)

# source files and headerfiles that should be included
%.o: %.cpp myfun.h mmio.h
	$(CC) $(CCFLAGS) -c $<

clean:
	rm edit $(OBJ) main.o
