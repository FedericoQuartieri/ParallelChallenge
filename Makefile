FLAGS = -I ${mkEigenInc} -O3 -fopenmp # -I ${mkLisInc} -L${mkLisLib} -llis
GCC = g++
MTX_ARGS = 83456566 3

all: execlean main.exe
	./main.exe ${MTX_ARGS}

%.exe: %.cpp
	@$(GCC) $(FLAGS) $< -o $@
	

execlean:
	@rm -f **.exe
