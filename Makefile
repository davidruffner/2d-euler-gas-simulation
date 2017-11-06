

default:
	gcc -Wall -o eulerSolve2D eulerSolve2D.c config.c myconfig.c euler.c -lm
	gcc -Wall -o eulerSolve2DnoGrav eulerSolve2DnoGrav.c config.c myconfig.c euler.c -lm
clean:
	rm -r eulerSolve2D eulerSolve2DnoGrav
	rm *~
