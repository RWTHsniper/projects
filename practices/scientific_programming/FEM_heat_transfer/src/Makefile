CC = icpc
SOURCE = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCE))
EXECUTABLE = 2d_Unsteady_Diffusion
CFLAGS =-O0 -g  -Wno-deprecated -Wall -I/usr/include/vtk
#CFLAGS =-O3 -Wno-deprecated -Wall -I/usr/include/vtk
LDFLAGS = -L/usr/lib64/vtk  -L$(TVLIB) -ltvheap_64 -Wl,-rpath,$(TVLIB)
LIBS = -lvtkCommon -lvtkFiltering -lvtkGraphics -lvtkIO -lvtkRendering -lvtkWidgets

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(LIBS) $(OBJECTS) -o $(EXECUTABLE)
	@echo DONE!

-include $(OBJECTS:.o=.d)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $*.cpp -o $*.o
	@$(CC) -MM $(CFLAGS) $*.cpp > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -rf *.o *.d $(EXECUTABLE)
	@echo ALL CLEANED UP!

rebuild:
	make clean
	make
