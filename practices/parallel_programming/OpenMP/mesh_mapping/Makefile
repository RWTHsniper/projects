CCFLAGS= -O3 -fopenmp
LDFLAGS=
OBJ=io.o projection.o definitions.o tests.o
EXE=pr2

all: $(EXE)

$(EXE): $(OBJ) main.o
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.c interfaces.h
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -rf $(OBJ) $(EXE) main.o
