FLAGS = -O3 -fopenmp # -I ${mkLisInc} -L${mkLisLib} -llis
GCC = g++
MTX_ARGS = 12345678 10

all: execlean main.exe
	./main.exe ${MTX_ARGS}

%.exe: %.cpp
	@$(GCC) $(FLAGS) $< -o $@
	

execlean:
	@rm -f **.exe
