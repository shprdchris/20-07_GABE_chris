COMPILER = g++
OFLAG = -O3
FLAGS = -m64 -g -fopenmp $(OFLAG)
LINKS= -L/opt/local/lib -lfftw3l_threads -lfftw3l_omp -lfftw3l  

do: fclean dirmake oclean compile 

xcmake : do
	
compile: g2header.h g2parameters.h g2model.o g2functions.o g2spectra.o g2output.o g2main.o g2init.o 
	$(COMPILER) $(FLAGS) g2model.o g2functions.o g2spectra.o g2output.o g2init.o g2main.o $(LINKS) -o gabe

clean : fclean oclean

g2main.o: g2header.h g2parameters.h g2main.cpp
	$(COMPILER) -c $(FLAGS)  g2main.cpp

g2spectra.o: g2header.h g2parameters.h g2spectra.cpp
	$(COMPILER) -c $(FLAGS)  g2spectra.cpp

g2model.o: g2header.h g2parameters.h g2model.cpp
	$(COMPILER) -c $(FLAGS)  g2model.cpp

g2functions.o: g2header.h g2parameters.h g2functions.cpp
	$(COMPILER) -c $(FLAGS)  g2functions.cpp

g2output.o: g2header.h g2parameters.h g2output.cpp
	$(COMPILER) -c $(FLAGS)  g2output.cpp

g2init.o: g2header.h g2parameters.h g2init.cpp
	$(COMPILER) -c $(FLAGS)  g2init.cpp

oclean: 
	rm -f *.o 

fclean:
#	rm -r gabe2_0.dSYM
	rm -f gabe2_0 *.dat *.txt
	rm -f ./slices/*.dat
	
dirmake:
	@ if test -d slices; then echo directory "'slices'" exists; else mkdir slices; echo made directory "'slices'"; fi	


run: clean compile
	./gabe
