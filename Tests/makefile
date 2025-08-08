CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6`

test1: test1.o
	$(CC) -o test1 test1.o $(LDFLAGS)

test1.o: test1.cpp
	$(CC) -c -o test1.o test1.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o test1
