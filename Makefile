FLAGS = -O3 -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp# -I ${mkLisInc} -L${mkLisLib} -llis
GCC = clang++
#######max_size;max_cutoff;passo;diveder for avg
MTX_ARGS = 1000000 15 10 8

all: allclean main.exe
	./main.exe ${MTX_ARGS}

%.exe: %.cpp
	@$(GCC) $(FLAGS) $< -o $@
	

allclean: pngclean csvclean execlean
	@rm -f plot_commands.gp

pngclean:
	@rm -f **.png


csvclean:
	@rm -f **.csv

execlean:
	@rm -f **.exe
