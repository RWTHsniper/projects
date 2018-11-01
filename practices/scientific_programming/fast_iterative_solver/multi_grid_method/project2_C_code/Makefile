# We require the C99 standard for testing.
CCFLAGS=-std=c99 -O2
LDFLAGS=-lm
# list all of object files
OBJ= myfun.o
EXE= pr2
# compiler
CC= gcc

all: $(EXE)

$(EXE): $(OBJ) main.o
	$(CC) $(CCFLAGS) -o $@ $^ $(LDFLAGS)

# headerfile that should be included
%.o: %.c myfun.h
	$(CC) $(CCFLAGS) -c $<

clean:
	rm edit $(OBJ) main.o
