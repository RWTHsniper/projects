all : obs.exe f.exe fd.exe
	./obs.exe 10 30
	./f.exe 10 30
	./fd.exe 10 30

%.exe : %.cpp f.hpp
	g++ $< -o$@

clean :
	rm -fr *.exe O.txt

.PHONY: clean
