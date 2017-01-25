CFLAGS = `pkg-config --cflags opencv`
LIBS = `pkg-config --libs opencv`
OP = -Ofast

Proc: Proc.cpp Sobel.o Median.o 
	g++ $(CFLAGS) Proc.cpp Sobel.o Median.o $(LIBS) -o Proc $(OP) -lpthread

Sobel.o: Sobel.c
	g++ $(CFLAGS) -c Sobel.c $(OP)

Median.o: Median.c
	g++ $(CFLAGS) -c Median.c $(OP) -lpthread

clean:
	rm -f *~ *.o Proc
