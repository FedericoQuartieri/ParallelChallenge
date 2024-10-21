FLAGS = -O3 -fopenmp # -I ${mkLisInc} -L${mkLisLib} -llis
FLAGS2 = -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp
GCC = g++
#######max_size;max_cutoff;passo;diveder for avg
MTX_ARGS = 1009000 12 13 8

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
